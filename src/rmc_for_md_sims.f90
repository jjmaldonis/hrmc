
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


program rmc

    use rmc_global
    use readinputs
    use model_mod
    use gr_mod
    use fem_mod
    use rmc_functions
    use eam_mod
    implicit none
    include 'mpif.h'
    type(model) :: m
    character (len=256) :: model_filename, c
    integer :: len
    character (len=256) :: vk_filename
    character (len=256) :: param_filename  
    character (len=256) :: outbase
    character (len=512) :: comment
    logical, dimension(4) :: used_data_sets
    !logical :: not_too_short
    real :: temperature
    real :: max_move
    real :: Q, res
    real :: pixel_distance
    real, dimension(4) :: weights
    real, pointer, dimension(:) :: gr_e,r_e,gr_e_err
    real, pointer, dimension(:) :: gr_n,r_n
    real, pointer, dimension(:) :: gr_x,r_x
    real, pointer, dimension(:) :: vk, vk_exp, k, vk_exp_err, v_background
    real, pointer, dimension(:,:) :: cutoff_r 
    real, pointer, dimension(:,:) :: scatfact_e
    real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: chi2_old, chi2_new, scale_fac, del_chi, beta, chi2_gr, chi2_vk
    real :: rmin_e, rmax_e
    real :: rmin_n, rmax_n
    real :: rmin_x, rmax_x
    real :: R
    integer :: i, j
    integer :: w
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: fem_algorithm
    integer :: istat, status2
    integer :: total_steps
    integer :: iseed2
    real :: randnum
    real :: te1, te2
    logical :: square_pixel, use_femsim
    doubleprecision :: t0, t1 !timers

    write(*,*)
    write(*,*) "This is the dev version of rmc!"
    write(*,*)

    call mpi_init(mpierr)
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

    call get_command_argument(1, c, len, istat)
    if (istat == 0) then
        model_filename = trim(c)
    else
        write (*,*) 'Cannot read model file name.  Exiting.'
        stop
    end if
    call get_command_argument(2, c, len, istat)
    if (istat == 0) then
        param_filename = trim(c)
    else
        write (*,*) 'Cannot read param file name.  Exiting.'
        stop
    end if
    call get_command_argument(3, c, len, istat)
    if (istat == 0) then
        vk_filename = trim(c)
    else
        write (*,*) 'Cannot read vk file name.  Exiting.'
        stop
    end if
    outbase = ""

    ! Read input model
    call read_model(model_filename, comment, m, istat)
    call check_model(m, istat)
    ! recenter_model is called in read_model so it doesnt need to be called
    ! agian I dont think. Its small though.
    call recenter_model(0.0, 0.0, 0.0, m)

    ! Read input parameters
    allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
    call read_inputs(param_filename,temperature, max_move, cutoff_r, used_data_sets, weights, gr_e, r_e, gr_e_err, gr_n, r_n, &
    gr_x, r_x, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, scale_fac, Q, fem_algorithm, pixel_distance, total_steps, &
    rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

    res = 0.61/Q
    nk = size(k)
    beta=1./((8.6171e-05)*temperature)

    iseed2 = 104756

    square_pixel = .TRUE. ! RMC uses square pixels, not round.
    use_femsim = .FALSE.

    !!!call read_eam(m)   
    !!!call eam_initial(m, te1)

    !!!write(*,*)"te1=", te1

    ! Initialize and calculate initial gr
    call scatt_power(m,used_data_sets,istat)
    !write (*,*) 'Znum reduced, scatt pow initialized.'
    !!!call gr_initialize(m,r_e,gr_e,r_n,gr_n,r_x,gr_x,used_data_sets,istat)

    ! call gr_no_hutch(m,used_data_sets) ! Jason 20130725 bc not in rmc in Jinwoos 040511c_t1

    ! Initialize and calculate initial vk
    call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat,  square_pixel)

    !call print_sampled_map(m, res, square_pixel)

    allocate(vk(size(vk_exp)))

    call fem(m, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)

    if(myid.eq.0)then
        ! Write initial gr
        !open(unit=51,file=outbase//"_gr_initial.txt",form='formatted',status='unknown')
        !    do i=1, mbin_x
        !        R = del_r_x*(i)-del_r_x
        !        write(51,*)R, gr_x_sim_cur(i)
        !    enddo
        !close(51)
        ! Write initial vk 
        open(unit=52,file=vk_filename,form='formatted',status='unknown')
            do i=1, nk
                !write (*,*) k(i),vk(i)
                write(52,*)k(i),vk(i)
            enddo
        close(52)
    endif

stop

    ! Initial chi2
    chi2_old = chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, vk_exp, vk_exp_err, &
        gr_e_sim_cur, gr_n_sim_cur, gr_x_sim_cur, vk, scale_fac,&
        rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)
    chi2_old = chi2_old + te1 

    e2 = e1

    if(myid.eq.0)then
       write(*,*)"initial"
       write(*,*)chi2_old, chi2_gr, chi2_vk, te1, temperature
       write(*,*)"  i, chi2_gr, chi2_vk, te1, temperature"
    endif

    write(*,*)
    write(*,*)"Initialization complete. Starting Monte Carlo."
    write(*,*)
    t0 = mpi_wtime()
    ! RMC step begins

    !*********update here*******************
    i=1315708
    !***************************************

    ! Fortran's max int is 2147483647, and 2147483647 + 1 = -2147483648, so this
    ! loop will stop at i= one after i=2147483647.
    ! The loops also stops if the temperature goes below a certain temp (30.0).
    do while (i >0)
        i=i+1

        call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
        ! check_curoffs returns false if the new atom placement is too close to
        ! another atom. Returns true if the move is okay. (hard shere cutoff)
        do while( .not. check_cutoffs(m,cutoff_r,w) )
            ! Check_cutoffs returned false so reset positions and try again.
            m%xx(w) = xx_cur
            m%yy(w) = yy_cur
            m%zz(w) = zz_cur
            call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
        end do

        ! Update hutches, data for chi2, and chi2/del_chi
        call hutch_move_atom(m,w,xx_new, yy_new, zz_new)
        call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
        !call gr_hutch_mc(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new,used_data_sets,istat) ! Jason 20130725 bc not in rmc in Jinwoos 040511c_t1
        call fem_update(m, w, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat, square_pixel)
        write(*,*) "Finished updating eam, gr, and fem data."
        
        chi2_new = chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, vk_exp, vk_exp_err,&
            gr_e_sim_new, gr_n_sim_new, gr_x_sim_new, vk, scale_fac,&
            rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)

        chi2_new = chi2_new + te2
        del_chi = chi2_new - chi2_old

        call mpi_bcast(del_chi, 1, mpi_real, 0, mpi_comm_world, mpierr)
           
        ! TODO Jason - Check if this is being done as well as it could be.
        ! the accept and reject functions shouldn't both have the same do loops
        ! in them.
        !write(*,*) "Accepting or rejecting move."
        randnum = ran2(iseed2)
        ! Test if the move should be accepted or rejected based on del_chi
        if(del_chi <0.0)then
            ! Accept the move
            e1 = e2
            !call accept_gr(m, used_data_sets) ! Jason 20130725 bc not in rmc in Jinwoos 040511c_t1
            call fem_accept_move(mpi_comm_world)
            if(myid.eq.0)then
                write(*,*)i, chi2_gr, chi2_vk, te2, temperature
            endif
            chi2_old = chi2_new
            write(*,*) "MC move accepted outright."
        else
            ! Based on the random number above, even if del_chi is negative, decide
            ! whether to move or not (statistically).
            if(log(1.-randnum)<-del_chi*beta)then
                ! Accept move
                e1 = e2
                !call accept_gr(m, used_data_sets) ! Jason 20130725 bc not in rmc in Jinwoos 040511c_t1
                call fem_accept_move(mpi_comm_world)
                if(myid.eq.0)then
                    write(*,*)i, chi2_gr, chi2_vk, te2, temperature
                endif
                chi2_old = chi2_new
                write(*,*) "MC move accepted due to probability. del_chi*beta = ", del_chi*beta
            else
                ! Reject move
                e2 = e1
                call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
                !call reject_gr(m,used_data_sets) ! Jason 20130725 bc not in rmc in Jinwoos 040511c_t1
                call fem_reject_move(m, mpi_comm_world)
                write(*,*) "MC move rejected."
            endif
        endif

        ! Every 50,000 steps lower the temp, max_move, and reset beta.
        if(mod(i,50000)==0)then
            temperature = temperature * sqrt(0.7)
            max_move = max_move * sqrt(0.94)
            beta=1./((8.6171e-05)*temperature)
            if(myid.eq.0)then
                if(temperature.lt.30.0)then
                    stop
                endif
            endif
        endif

        ! Periodically save data. Are we saving the energy here??? We should!
        ! The energy will tell us if there is a steady decline, or if it
        ! statistically jumped way up. We could also write del_chi (or instead).
        if(mod(i,1000)==0)then
            if(myid.eq.0)then
                open(31,file=outbase//'_gr_update.txt',form='formatted',status='unknown')
                open(32,file=outbase//'_vk_update.txt',form='formatted',status='unknown')
                open(33,file=outbase//'_model_update.txt',form='formatted',status='unknown')
                open(34,file=outbase//'_energy_function.txt',form='formatted',status='unknown')
                open(35,file=outbase//'_time_elapsed.txt',form='formatted',status='unknown')
                ! Write to gr_update        
                !do j=1, mbin_x
                !    R = del_r_x*(j)-del_r_x
                !    write(31,*)R, gr_x_sim_new(j)
                !enddo
                ! Write to vk_update
                do j=1, nk
                    write(32,*)k(j),vk(j)
                enddo
                ! Write to model_update
                write(33,*)"updated model"
                write(33,*)m%lx,m%ly,m%lz
                do j=1,m%natoms
                    write(33,*)m%znum(j), m%xx(j), m%yy(j), m%zz(j)
                enddo
                write(33,*)"-1"
                ! Write to energy_function
                write(34,*) i, te1
                ! Write to time_elapsed
                t1 = mpi_wtime()
                write (34,*) i, t1-t0
                ! Close files
                close(31)
                close(32)
                close(33)
                close(34)
                close(35)
            endif
        endif
    ENDDO

    ! The rmc loop finished. Write final data.
    if(myid.eq.0)then
        t1 = mpi_wtime()
        write(*,*)"time=", t1-t0, "sec"
        write(*,*)t1, t0

        ! Write final gr
        !open(unit=53,file=outbase//"_gr_update_final.txt",form='formatted',status='unknown')
        !do i=1, mbin_e
        !    R = del_r_e*(i)-del_r_e
        !    write(53,*)R, gr_e_sim_new(i)
        !enddo
        !close(53)
        
        ! Write final vk
        open(unit=54,file=outbase//"_vk_update_final.txt",form='formatted',status='unknown')
        do i=1, nk
            write(54,*)k(i),vk(i)
        enddo
        close(54)
        
        ! Write final model
        open(unit=55,file=outbase//"_model_update_final.txt",form='formatted',status='unknown')
        write(55,*)"updated model"
        write(55,*)m%lx,m%ly,m%lz
        do i=1,m%natoms
            write(55,*)m%znum(i), m%xx(i), m%yy(i), m%zz(i)
        enddo
        write(55,*)"-1"
        close(55)
        ! Write final energy.
        open(56,file=outbase//'_energy_function.txt',form='formatted', status='unknown')
        write(56,*) te1, i
        close(56)
        ! Write final time spent.
        open(57,file=outbase//'_time_elapsed.txt',form='formatted',status='unknown')
        t1 = mpi_wtime()
        write (34,*) i, t1-t0
        close(57)
    endif
    call mpi_finalize(mpierr)

end program rmc
