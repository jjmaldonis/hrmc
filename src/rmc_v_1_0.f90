
! Reverse Monte Carlo structural simulation
!
! Performs RMC refinement against reduced density
! function data from electron, x-ray, or neutron
! diffraction data and against fluctuation electron
! microscopy V(k) data.
!
! Voyles research group: Jinwoo Hwang, Feng, Yi, and
! Paul Voyles, begun 12/16/08
! 
! merged with version of fy, 12/20/08, pmv
! more program flow added 12/20/08 pmv
! add new variable: electron scattering coefficients for calculating
! which are a_e, b_e, c_e and d_e
! add new module scattering_factors, add USE without :: following them for
! visual fortran debug by fy on 12/23
! need to not define gr_sim etc. arrays
! Re-written by Jinwoo Hwang - 03/20/2009
! Hard sphere cutoff applied.
! Debugged and tested by Jinwoo Hwang - 03/25/2009
! Vk allocated before initial Vk - JWH 04/02/2009
! hutch_move_atom added in reject part - JWH 04/15/2009.
! gr_e_err is added - jwh 04/25/2009
! See commit history for changes by jjm on Github: (check all branches)
!   https://github.com/refreshx2/OdieCode/tree/master/060911_rmc_eam_gr_vk_t4

program rmc

#ifdef USE_LMP
    use LAMMPS
#endif
    use omp_lib
    use rmc_global
    use readinputs
    use model_mod
    use gr_mod
    use fem_mod
    use rmc_functions
    use eam_mod
    implicit none
    include 'mpif.h'
    ! LAMMPS objects
#ifdef USE_LMP
    type (C_ptr) :: lmp
    real (C_double), pointer :: te1 => NULL()
    real (C_double), pointer :: te2 => NULL()
    character (len=512) :: lmp_cmd_str
#endif
    ! RMC / Femsim objects
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: param_filename  
    character (len=256) :: outbase
    character (len=256) :: jobID, c, step_str
    character (len=512) :: comment
    character (len=256) :: time_elapsed, vki_fn, vku_fn, vkf_fn, output_model_fn, energy_fn, final_model_fn, gri_fn, chi_squared_file, acceptance_rate_fn, beta_fn
    logical, dimension(4) :: used_data_sets
    real :: temperature
    real :: max_move
    real :: Q, res
    real :: pixel_distance
    real, dimension(4) :: weights
    real, pointer, dimension(:) :: gr_e,r_e,gr_e_err
    real, pointer, dimension(:) :: gr_n,r_n
    real, pointer, dimension(:) :: gr_x,r_x
    real, pointer, dimension(:) :: vk, vk_exp, k, vk_exp_err, v_background, vk_as !vk_as is for autoslice/multislice (whatever it's called)
    real, pointer, dimension(:,:) :: cutoff_r 
    real, pointer, dimension(:,:) :: scatfact_e
    real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: chi2_old, chi2_new, scale_fac, scale_fac_initial, del_chi, beta, chi2_gr, chi2_vk, chi2_no_energy
    real :: rmin_e, rmax_e
    real :: rmin_n, rmax_n
    real :: rmin_x, rmax_x
    real :: R
    integer :: i, j
    integer :: w
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: fem_algorithm
    integer :: istat, status2, length
    integer :: total_steps
    integer :: iseed2
    real :: randnum
#ifndef USE_LMP
    real :: te1, te2 ! For eam, not lammps version of energy
#endif
    logical :: square_pixel, accepted, use_rmc, use_multislice
    integer :: ipvd, nthr
    doubleprecision :: t0, t1, t2 !timers
    real :: x! This is the parameter we will use to fit vsim to vas.
    integer, dimension(1000) :: acceptance_array
    real :: avg_acceptance = 1.0

    !------------------- Program setup. -----------------!

#ifdef TIMING
    write(*,*) "Using timing version of code!"
#endif

    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr) !http://www.open-mpi.org/doc/v1.5/man3/MPI_Init_thread.3.php
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

    nthr = omp_get_max_threads()
    if(myid.eq.0)then
        write(*,*)
        write(*,*) "Using", numprocs, "processors."
        write(*,*) "OMP found a max number of threads of", nthr
        write(*,*)
    endif

#ifdef USE_LMP
    if(myid.eq.0) then
        call lammps_open ('lmp -log log.simple -screen none', mpi_comm_world, lmp)
        call lammps_file (lmp, 'lmp_energy.in')
    endif
#endif

if(myid .eq. 0) then
    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        jobID = '_temp'
    end if
    call get_command_argument(2, c, length, istat)
    !if (istat == 0) then
    !    model_filename = trim(c)
    !else
    !    stop "istat for get_command_arg was nonzero", model_filename
    !end if
    !model_filename = trim(model_filename)

    ! Set output filenames.
    outbase = ""
    write(time_elapsed, "(A12)") "time_elapsed"
    time_elapsed = trim(trim(time_elapsed)//jobID)//".txt"
    write(vki_fn, "(A10)") "vk_initial"
    vki_fn = trim(trim(vki_fn)//jobID)//".txt"
    write(vku_fn, "(A9)") "vk_update"
    vku_fn = trim(trim(vku_fn)//jobID)//".txt"
    write(vkf_fn, "(A8)") "vk_final"
    vkf_fn = trim(trim(vkf_fn)//jobID)//".txt"
    write(output_model_fn, "(A12)") "model_update"
    output_model_fn = trim(trim(output_model_fn)//jobID)//".txt"
    write(final_model_fn, "(A11)") "model_final"
    final_model_fn = trim(trim(final_model_fn)//jobID)//".txt"
    write(energy_fn, "(A6)") "energy"
    energy_fn = trim(trim(energy_fn)//jobID)//".txt"
    write(gri_fn, "(A10)") "gr_initial"
    gri_fn = trim(trim(gri_fn)//jobID)//".txt"
    write(chi_squared_file, "(A11)") "chi_squared"
    chi_squared_file = trim(trim(chi_squared_file)//jobID)//".txt"
    write(acceptance_rate_fn, "(A15)") "acceptance_rate"
    acceptance_rate_fn = trim(trim(acceptance_rate_fn)//jobID)//".txt"
    write(beta_fn, "(A14)") "beta_weighting"
    beta_fn = trim(trim(beta_fn)//jobID)//".txt"
endif

    !------------------- Read inputs and initialize. -----------------!

    ! Set input filenames.
    param_filename = 'param_file.in'
    !call execute_command_line("cat param_file.in")
    
    ! Start timer.
    t0 = omp_get_wtime()

    ! Read input model
    call read_model(param_filename, model_filename, comment, m, istat)
    call check_model(m, istat)
    call recenter_model(0.0, 0.0, 0.0, m)

    ! Read input parameters
    allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
    call read_inputs(param_filename,model_filename,temperature, max_move, cutoff_r, &
        used_data_sets, weights, gr_e, r_e, gr_e_err, gr_n, r_n, gr_x, &
        r_x, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, &
        scale_fac_initial, Q, fem_algorithm, pixel_distance, total_steps, &
        rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

    if(myid .eq. 0) then
    write(*,*) "Model filename: ", trim(model_filename)
    write(*,*)
    endif

    scale_fac = scale_fac_initial
    res = 0.61/Q
    nk = size(k)
    beta=1./((8.6171e-05)*temperature)
    iseed2 = 104756
    if(myid.eq.0) write(*,*) "random number generator seed =", iseed2
    square_pixel = .TRUE. ! RMC uses square pixels, not round.
    use_rmc = .TRUE.
    use_multislice = .FALSE.

#ifdef USE_LMP
    call lammps_command (lmp, 'run 0')
    call lammps_extract_compute (te1, lmp, 'pot', 0, 0)
#else
    call read_eam(m)
    call eam_initial(m,te1)
#endif
    if(myid .eq. 0) write(*,*) "Energy = ", te1

    call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat,  square_pixel)
    allocate(vk(size(vk_exp)), vk_as(size(vk_exp)))
    vk = 0.0; vk_as = 0.0

    ! Print warning message if we are using too many cores.
    if(myid.eq.0) then
        call print_sampled_map(m, res, square_pixel)
        if(pa%npix /= 1) then
            if(numprocs > 3*nrot) write(*,*) "WARNING: You are using too many cores!"
        else
            if(numprocs > nrot) write(*,*) "WARNING: You are using too many cores!"
        endif
    endif

    !------------------- Call femsim. -----------------!

    ! Fem updates vk based on the intensity calculations and v_background.
    call fem(m, res, k, vk, vk_as, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)

    if(myid.eq.0)then
        write(*,*) "Femsim took", t1-t0, "seconds on processor", myid

        ! Write initial vk 
        open(unit=52,file=trim(vki_fn),form='formatted',status='unknown')
!            write(*,*) "Writing V(k)."
            do i=1, nk
!                write (*,*) k(i),vk(i)
                write(52,*) k(i),vk(i)
            enddo
        close(52)
    endif

    !------------------- Start RMC. -----------------!

    call mpi_barrier(mpi_comm_world, mpierr)

    if(use_rmc) then ! End here if we only want femsim. Set the variable above.

        ! Calculate initial chi2
        chi2_no_energy = chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, vk_exp, vk_exp_err, &
            gr_e_sim_cur, gr_n_sim_cur, gr_x_sim_cur, vk, scale_fac,&
            rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)

        chi2_old = chi2_no_energy + te1
#ifndef USE_LMP
        e2 = e1 ! eam
#endif

        if(myid.eq.0)then
            write(*,*)
            write(*,*) "Initialization complete. Starting Monte Carlo."
            write(*,*) "Initial Conditions:"
            write(*,*) "   Step =      ", i
            write(*,*) "   Energy =     ", te1
            write(*,*) "   Temperature =", temperature
            write(*,*)
            ! Reset time_elapsed, energy_function, chi_squared_file
            open(35,file=trim(time_elapsed),form='formatted',status='unknown')
                t1 = omp_get_wtime()
                write(35,*) numprocs, "processors are being used."
                write(35,*) "Step, Time elapsed, Avg time per step, This step's time"
            close(35)
            open(36,file=trim(chi_squared_file),form='formatted',status='unknown')
                write(36,*) "step, chi2, energy"
                write(36,*) i, chi2_no_energy, te1
            close(36)
            open(37,file=trim(acceptance_rate_fn),form='formatted',status='unknown',access='append')
                write(37,*) "step, acceptance rate averaged over last 1000 steps"
            close(37)
        endif


        t0 = omp_get_wtime()
        i=0
        ! RMC loop begins. The loop never stops.
        do while (i .ge. 0)
            i=i+1
            t2 = omp_get_wtime()

            if(myid .eq. 0) write(*,*) "Starting step", i

#ifdef TIMING
            if( i > 100) then
                if(myid .eq. 0) write(*,*) "STOPPING MC AFTER 100 STEPS"
                call mpi_finalize(mpierr)
                stop ! Stop after 100 steps for timing runs.
            endif
#endif

            call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            ! check_curoffs returns false if the new atom placement is too close to
            ! another atom. Returns true if the move is okay. (hard shere cutoff)
            do while( .not. check_cutoffs(m,cutoff_r,w) )
                ! Check_cutoffs returned false so reset positions and try again.
                m%xx%ind(w) = xx_cur
                m%yy%ind(w) = yy_cur
                m%zz%ind(w) = zz_cur
                call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            end do
            ! Update hutches, data for chi2, and chi2/del_chi
            call hutch_move_atom(m,w,xx_new, yy_new, zz_new)
    
#ifdef USE_LMP
            write(lmp_cmd_str, "(A9, I4, A3, F, A3, F, A3, F)") "set atom ", w, " x ", xx_new, " y ", yy_new, " z ", zz_new
            call lammps_command(lmp, trim(lmp_cmd_str))
            call lammps_command (lmp, 'run 0')
            call lammps_extract_compute (te2, lmp, 'pot', 0, 0)
#else
            call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
#endif
            if(myid .eq. 0) write(*,*) "Energy = ", te2

            ! Calculate a randnum for accept/reject
            randnum = ran2(iseed2)
            accepted = .true.
            if(log(1.0-randnum) > -(te2-chi2_old)*beta) then
                accepted = .false.
                call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
#ifndef USE_LMP
                e2 = e1 ! eam
#else
                write(lmp_cmd_str, "(A9, I4, A3, F, A3, F, A3, F)") "set atom ", w, " x ", xx_cur, " y ", yy_cur, " z ", zz_cur
                call lammps_command(lmp, trim(lmp_cmd_str))
#endif
                if(myid.eq.0) write(*,*) "MC move rejected solely due to energy."
            endif

            if(accepted) then
            ! Use multislice every 10k steps if specified.
            if(use_multislice .and. mod(i,10000) .eq. 0) then
                call fem_update(m, w, res, k, vk, vk_as, v_background, scatfact_e, mpi_comm_world, istat, square_pixel, .true.)
                call update_scale_factor(scale_fac, scale_fac_initial, vk, vk_as)
            else
                call fem_update(m, w, res, k, vk, vk_as, v_background, scatfact_e, mpi_comm_world, istat, square_pixel, .false.)
                !write(*,*) "I am core", myid, "and I have exited from fem_update into the main rmc block."
            endif

            chi2_no_energy = chi_square(used_data_sets,weights,gr_e, gr_e_err, &
                gr_n, gr_x, vk_exp, vk_exp_err, gr_e_sim_new, gr_n_sim_new, &
                gr_x_sim_new, vk, scale_fac, rmin_e, rmax_e, rmin_n, rmax_n, &
                rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)

            chi2_new = chi2_no_energy + te2
            del_chi = chi2_new - chi2_old
            call mpi_barrier(mpi_comm_world, mpierr) ! Pretty sure this is unnecessary

            ! Test if the move should be accepted or rejected based on del_chi
            if(del_chi <0.0)then
                ! Accept the move
                call fem_accept_move(mpi_comm_world)
#ifndef USE_LMP
                e1 = e2 ! eam
#endif
                chi2_old = chi2_new
                accepted = .true.
                if(myid.eq.0) write(*,*) "MC move accepted outright."
            else
                ! Based on the random number above, even if del_chi is negative, decide
                ! whether to move or not (statistically).
                if(log(1.-randnum)<-del_chi*beta)then
                    ! Accept move
                    call fem_accept_move(mpi_comm_world)
#ifndef USE_LMP
                    e1 = e2 ! eam
#endif
                    chi2_old = chi2_new
                    accepted = .true.
                    if(myid.eq.0) write(*,*) "MC move accepted due to probability. del_chi*beta = ", del_chi*beta
                else
                    ! Reject move
                    call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                    call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
                    call fem_reject_move(m, w, xx_cur, yy_cur, zz_cur, mpi_comm_world)
#ifndef USE_LMP
                    e2 = e1 ! eam
#else
                    write(lmp_cmd_str, "(A9, I4, A3, F, A3, F, A3, F)") "set atom ", w, " x ", xx_cur, " y ", yy_cur, " z ", zz_cur
                    call lammps_command(lmp, trim(lmp_cmd_str))
#endif
                    accepted = .false.
                    if(myid.eq.0) write(*,*) "MC move rejected."
                endif
            endif
            endif ! if(accepted) from above

            if(accepted) then
                acceptance_array(mod(i,1000)+1) = 1
            else
                acceptance_array(mod(i,1000)+1) = 0
            endif
            if(i .gt. 1000) avg_acceptance = sum(acceptance_array)/1000.0
            ! Writing to 0 is stderr
            if(i .gt. 1000 .and. avg_acceptance .le. 0.05) write(0,*) "WARNING!  Acceptance rate is low:", avg_acceptance

            ! Periodically save data.
            if(myid .eq. 0) then
            if(mod(i,1000)==0)then
                ! Write to vk_update
                write(vku_fn, "(A9)") "vk_update"
                write(step_str,*) i
                vku_fn = trim(trim(trim(trim(vku_fn)//jobID)//"_")//step_str)//".txt"
                open(32,file=trim(vku_fn),form='formatted',status='unknown')
                    do j=1, nk
                        write(32,*)k(j),vk(j)
                    enddo
                close(32)
                ! Write to model_update
                ! This takes a bit of time.
                write(output_model_fn, "(A12)") "model_update"
                write(step_str,*) i
                output_model_fn = trim(trim(trim(trim(output_model_fn)//jobID)//"_")//step_str)//".xyz"
                open(33,file=trim(output_model_fn),form='formatted',status='unknown')
                    write(33,*)"updated model"
                    write(33,*)m%lx,m%ly,m%lz
                    do j=1,m%natoms
                        write(33,*)m%znum%ind(j), m%xx%ind(j), m%yy%ind(j), m%zz%ind(j)
                    enddo
                    write(33,*)"-1"
                close(33)
            endif
            if(mod(i,1)==0)then
                if(accepted) then
                    ! Write chi2 info
                    open(36,file=trim(chi_squared_file),form='formatted',status='unknown',access='append')
                        write(36,*) i, chi2_no_energy, te2
                    close(36)
                endif
            endif
            if(mod(i,1)==0)then
                ! Write to time_elapsed
                open(35,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
                    t1 = omp_get_wtime()
                    write (35,*) i, t1-t0, (t1-t0)/i, t1-t2
                close(35)
            endif
            if(mod(i,1000)==0 .and. i .gt. 100)then
                ! Write to acceptance rate
                open(40,file=trim(acceptance_rate_fn),form='formatted',status='unknown',access='append')
                    write(40,*) i, avg_acceptance
                close(40)
            endif
            endif ! myid == 0

            ! Every 50,000 steps lower the temp, max_move, and reset beta.
            if(mod(i,50000)==0)then
                temperature = temperature * sqrt(0.7)
                if(myid.eq.0) write(*,*) "Lowering temp to", temperature, "at step", i
                max_move = max_move * sqrt(0.94)
                beta=1./((8.6171e-05)*temperature)
            endif
            
        enddo !RMC do loop
        write(*,*) "Monte Carlo Finished!"

        ! The rmc loop finished. Write final data.
        if(myid.eq.0)then
            ! Write final vk
            open(unit=54,file=trim(vkf_fn),form='formatted',status='unknown')
            do i=1, nk
                write(54,*)k(i),vk(i)
            enddo
            close(54)
            ! Write final model
            open(unit=55,file=trim(final_model_fn),form='formatted',status='unknown')
            write(55,*)"final model"
            write(55,*)m%lx,m%ly,m%lz
            do i=1,m%natoms
                write(55,*)m%znum%ind(i), m%xx%ind(i), m%yy%ind(i), m%zz%ind(i)
            enddo
            write(55,*)"-1"; close(55)
            ! Write final energy.
            open(56,file=trim(energy_fn),form='formatted', status='unknown',access='append')
            write(56,*) i, te2
            close(56)
            ! Write final time spent.
            open(57,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
            t1 = omp_get_wtime()
            write (57,*) i, t1-t0
            write(57,*) "Finshed.", numprocs, "processors."
            close(57)
        endif
    endif ! Use RMC

#ifdef USE_LMP
    if(myid.eq.0) then
    call lammps_close (lmp)
    endif
#endif
    call mpi_finalize(mpierr)

end program rmc
