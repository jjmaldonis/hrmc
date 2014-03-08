!
!This module loads the experiment data.
!It was written by Jason Maldonis from the old read_inputs.f90 file on
!06/28/2013
!
module ReadInputs

contains

    subroutine read_inputs(param_filename,model_fn, temperature, max_move, cutoff_r, used_data_sets, weights, gr_e, r_e, gr_e_err, &
    gr_n, r_n, gr_x, r_x, V, k, V_err, V_background, ntheta, nphi, npsi, scale_fac, Q, fem_algorithm, pixel_distance, total_steps, &
    rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)
    !It reads input file names, temperature, maximum atom movement distance
    !cutoff distance, used_data_sets, weight factors calculating chi square
    !structure factor pairs derived from electron scattering, neutron scattering
    !and X-ray scattering
    !FEM data V(k), k and V_err triplets
    !Rotation number and fem algorithm
    !Pixel distance
    !status

    !param_filename=input file name containing initilizing parameters,such as temperature
    !temperature=beginning temperature for RMC
    !max_move=maximum allowable movement distance
    !cutoff_r=cut off distance for different pairs
    !used_data_sets=A logical array determine whether electron, neutron, X-ray or FEM data are available or not
    !weights=
    !gr_e and r_e=electron G(r) data
    !gr_e_err=
    !gr_n and r_n=neutron G(r) data
    !gr_x and r_x=X-ray G(r) data
    !V=FEM intensity variance
    !k=scattering vector
    !V_err=FEM measurement variance error
    !V_background=
    !ntheta, nphi, npsi=rotation number for calculating simulated V
    !scale_fac=
    !Q=resolution for FEM
    !fem_algorithm=which algorithm is used to do variance calculation
    !pixel_distance=pixel distance
    !total_steps=
    !rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x=
    !status2=whether param_filename is opened with success or not

        implicit none
        character (len=*), intent(in) :: param_filename  !Assume size array, be careful
        character (len=80), intent(out) :: model_fn
        real, intent(out) :: temperature
        real, intent(out) :: max_move
        real, intent(out), dimension(:,:) :: cutoff_r
        logical, intent(out), dimension(4) :: used_data_sets
        real, intent(out), dimension(4) :: weights
        real, pointer, dimension(:) :: gr_e,r_e,gr_e_err
        real, pointer, dimension(:) :: gr_n,r_n
        real, pointer, dimension(:) :: gr_x,r_x
        real, pointer, dimension(:) :: v, k, v_err, v_background
        real, intent(out) :: rmin_e, rmax_e
        real, intent(out) :: rmin_n, rmax_n
        real, intent(out) :: rmin_x, rmax_x
        integer, intent(out) :: ntheta, nphi, npsi
        real, intent(out) :: q
        real, intent(out) :: scale_fac
        integer, intent(out) ::fem_algorithm
        real, intent(out) :: pixel_distance
        integer, intent(out) ::total_steps
        integer, intent(out) :: status2

        !Variables declared in this subroutine, local variables
        !comment1=comment in param_filename
        character (len=80) comment1  
        character (len=80) scatteringfile ! The scattering file name, at most 20 characters
        character (len=80) femfile ! The fem data file name, at most 20 characters
            ! Note: electron, neutron, and x-ray data all use this name, but in order.
        integer filenamelength !The file name length in scatteringfile or femfile
        integer i
        real indicator_end !Indicator_end=-1 means the end of file reaches
        integer status1 !Indicate the status of opening file in this subroutine
        logical file_end !Indicate fileend has reached
            ! Files for electron, neutron, x-ray scattering and fem data
        integer num_line !Number of lines in each data file except comment line
        real, pointer, dimension(:) :: tempdata !Temperature data
        integer stat_allocate1, stat_allocate2, stat_allocate3,stat_allocate4 !Allocate status, 0 means success

        open(20, file=param_filename,iostat=status2, status='old')
        if(status2 .ne. 0) then !open fails
            print *, 'cannot open file with name: ', param_filename
            return
        endif
        read(20, '(a80)') comment1 !read comment and it is neglected in the input
        read(20, '(a80)') model_fn; model_fn = adjustl(model_fn)
        read(20, *) total_steps !WRONG. I changed this to starting step #
        read(20, *) temperature
        read(20, *) max_move
        read(20, * ) !commented line
        read(20, *) cutoff_r ! This is a nelements by nelements matrix

        ! Read each number in the file in column order
        ! i stands for which file we are using
        i=0

        !********************************************************************
        !*******************READ electron scattering file**********************
        !********************************************************************
        i=i+1
        num_line=0

        read(20, '(a)') scatteringfile
        scatteringfile = adjustl(scatteringfile)
        read(20, *) weights(i)
        read(20, *) rmin_e, rmax_e

        if(scatteringfile(1:1) .ne. '#') then !Electron scattering file exists
            used_data_sets(i) = .TRUE.
            filenamelength=len_trim(scatteringfile)
            file_end = .false.
            open(30,file=scatteringfile(1:filenamelength),iostat=status1,status='old') 
            if(status1 .eq. 0) then !Open succeeds
                read(30, '(a80)') comment1
                do while( .not. file_end) ! Count how many data pairs are in the file. -1 dentoes EOF
                    read(30, *) indicator_end
                    if(abs(indicator_end+1.0) .le. 1e-6) then !File end is reached
                        exit !Exit this loop
                    else
                        num_line = num_line + 1
                    endif
                enddo
                rewind(30) !Go to the beginning of the file
                !Read r_e and gr_e data
                !read r_e first, then gr_e
                read(30, '(a80)') comment1
                allocate(tempdata(3*num_line),stat=stat_allocate1)
                if(stat_allocate1 .eq. 0) then
                    read(30, *) tempdata
                    allocate(r_e(num_line), stat=stat_allocate2)
                    allocate(gr_e(num_line), stat=stat_allocate3)
                    allocate(gr_e_err(num_line), stat=stat_allocate4)   !JWH 04/25/2009
                    if ((stat_allocate2 .eq. 0) .and. (stat_allocate3 .eq. 0) .and. (stat_allocate4 .eq. 0)  ) then
                        r_e=tempdata(1:3*num_line:3)
                        gr_e=tempdata(2:3*num_line:3)
                        gr_e_err=tempdata(3:3*num_line:3)
                    else
                        print *, 'electron scattering part 2 or 3 fails!'
                        return
                    endif
                else
                        print *, 'electron scattering part allocation fail'
                        return
                endif
                deallocate(tempdata)
            else
                print *, 'Open electron scattering file fails.'
            endif
            close(30)
        else 
            used_data_sets(i) = .false.
        endif !Electron scattering file exists

        !*******************************************************************
        !*****************Read neutron scattering file**********************
        !*******************************************************************

        i=i+1
        num_line=0

        read(20, '(a)') scatteringfile
        scatteringfile = adjustl(scatteringfile)
        read(20, *) weights(i)
        read(20, *) rmin_n, rmax_n

        if(scatteringfile(1:1) .ne. '#') then !electron scattering file exists
            used_data_sets(i) = .true.
            filenamelength=len_trim(scatteringfile)
            file_end=.false.
            open(30,file=scatteringfile(1:filenamelength),iostat=status1,status='old') 
            if(status1 .eq. 0) then !open succeeds
                read(30, '(a80)') comment1
                do while( .not. file_end) ! Count how many data pairs are in the file. -1 denotes EOF
                    read(30, *) indicator_end
                    if(abs(indicator_end+1.0) .le. 1e-6) then !file end is reached
                        !print *, 'neutron scattering file end is reached!'
                        exit !exit this loop
                    else
                        num_line = num_line + 1
                        !read(30, *) 
                    endif
                enddo
                rewind(30) !go to the beginning of the file
                !read r_n and gr_n data
                !read r_n first, then gr_n
                read(30, '(a80)') comment1
                allocate(tempdata(2*num_line),stat=stat_allocate1)
                if(stat_allocate1 .eq. 0) then
                    read(30, *) tempdata
                    allocate(r_n(num_line), stat=stat_allocate2)
                    allocate(gr_n(num_line), stat=stat_allocate3)
                    if ((stat_allocate2 .eq. 0) .and. (stat_allocate3 .eq. 0)) then
                        r_n=tempdata(1:2*num_line:2)
                        gr_n=tempdata(2:2*num_line:2)
                    else
                        print *, 'neutron scattering part 2 or 3 fails!'
                        return
                    endif
                else
                    print *, 'neutron scattering part allocation fail'
                    return
                endif
                deallocate(tempdata)
            else
                print *, 'open neutron scattering file fails'
            endif
            close(30)
        else 
            used_data_sets(i) = .false.
        endif !neutron scattering file exists

        !*******************************************************************
        !*****************read x-ray scattering file**********************
        !*******************************************************************

        i=i+1
        num_line=0

        read(20, '(a)') scatteringfile
        scatteringfile = adjustl(scatteringfile)
        read(20, *) weights(i)
        read(20, *) rmin_x, rmax_x

        if(scatteringfile(1:1) .ne. '#') then !x-ray scattering file exists 
            used_data_sets(i) = .true.
            filenamelength=len_trim(scatteringfile)
            file_end=.false.
            open(30, file="xrd_zr50_kramer_reduced.txt", iostat=status1,status='old')
            if(status1 .eq. 0) then !open succeeds
                read(30, '(a80)') comment1
                do while( .not. file_end) ! Count how many data pairs are in the file. -1 denotes EOF
                    read(30, *) indicator_end
                    if(abs(indicator_end+1.0) .le. 1e-6) then !file end is reached
                        !print *, 'x-ray scattering file end is reached!'
                        exit !exit this loop
                    else
                        num_line = num_line + 1
                        !read(30, *) 
                    endif
                enddo                          
                rewind(30) !go to the beginning of the file
                !read r_x and gr_x data
                !read r_x first, then gr_x
                read(30, '(a80)') comment1
                allocate(tempdata(2*num_line),stat=stat_allocate1)
                if(stat_allocate1 .eq. 0) then
                    read(30, *) tempdata
                    allocate(r_x(num_line), stat=stat_allocate2)
                    allocate(gr_x(num_line), stat=stat_allocate3)
                    if ((stat_allocate2 .eq. 0) .and. (stat_allocate3 .eq. 0)) then
                        r_x=tempdata(1:2*num_line:2)
                        gr_x=tempdata(2:2*num_line:2)
                    else
                        print *, 'x-ray scattering part 2 or 3 fails!'
                        return
                    endif
                else
                    print *, 'x-ray scattering part allocation fail'
                    return
                endif
                deallocate(tempdata)
            else
                print *, 'open x-ray scattering file fails'
            endif
            close(30)
        else 
            used_data_sets(i) = .false.
        endif !X-ray scattering file exists

        !*******************************************************************
        !******************Read FEM file************************************
        !*******************************************************************

        i=i+1
        num_line=0

        read(20, '(a)') femfile
        femfile = adjustl(femfile)
        read(20, *) weights(i)
        read(20, *) scale_fac
        read(20, *) nphi, npsi, ntheta
        read(20, *) q
        read(20, *) fem_algorithm
        read(20, *) pixel_distance

        if(femfile(1:1) .ne. '#') then !fem file exists
            used_data_sets(i) = .true.
            filenamelength=len_trim(femfile)
            file_end=.false.
            open(30,file=femfile(1:filenamelength),iostat=status1,status='old') 
            if(status1 .eq. 0) then !open succeeds
                read(30, '(a80)') comment1 !first line is comment
                ! Count how many data pairs are in the file. -1 denotes EOF
                do while( .not. file_end)
                    read(30, *) indicator_end
                    if(abs(indicator_end+1.0) .le. 1e-6) then !file end is reached
                        exit !exit this loop
                    else
                        num_line = num_line + 1
                    endif
                enddo
                rewind(30) !go to the beginning of the file
                read(30, '(a80)') comment1
                allocate(tempdata(4*num_line),stat=stat_allocate1)
                !read k ,v, and v_err data
                !read k first, then v, last v_err
                if(stat_allocate1 .eq. 0) then
                    read(30, *) tempdata
                    allocate(k(num_line), stat=stat_allocate2)
                    allocate(v(num_line), stat=stat_allocate3)
                    allocate(v_err(num_line), stat=stat_allocate4)
                    allocate(v_background(num_line))

                    if ((stat_allocate2 .eq. 0) .and. (stat_allocate3 .eq. 0) .and. (stat_allocate4 .eq. 0)) then
                        k=tempdata(1:4*num_line:4)
                        v=tempdata(2:4*num_line:4)
                        v_err=tempdata(3:4*num_line:4)
                        v_background=tempdata(4:4*num_line:4)
                    else
                        print *, 'fem part 2 or 3, or 4 fails!'
                        return
                    endif !allocate2, 3 and 4
                else
                    print *, 'fem part allocation fail'
                    return
                endif
            deallocate(tempdata)
            else
                print *, 'open fem file fails', femfile(1:filenamelength)
            endif
            close(30)
        else 
          used_data_sets(i) = .false.
        endif !FEM file exists
        close(20)
    end subroutine read_inputs

end module readinputs

