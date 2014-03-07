!************************************************************
!******A module is defined here to load experiment data******
!*******This subroutine is written by FY on 12/18/2008*******
!************************************************************
MODULE ReadInputs

CONTAINS

!It reads input file names, temperature, maximum atom movement distance
!cutoff distance, used_data_sets, weight factors calculating chi square
!structure factor pairs derived from electron scattering, neutron scattering
!and X-ray scattering
!FEM data V(k), k and V_err triplets
!Rotation number and fem algorithm
!Pixel distance
!status

!   change log
!   gr_e_err is added - jwh 04/25/2009
!   modified to print what's read from the fem file right away - JWH 05/08/2009
!   The order of reading FEM angles are changed to nphi, npsi, ntheta = JWH 09/04/09




!If this subroutine used a module or modules, do not include implicit none line

SUBROUTINE read_inputs(param_filename,temperature, max_move, cutoff_r, used_data_sets, weights, gr_e, r_e, gr_e_err, &
gr_n, r_n, gr_x, r_x, V, k, V_err, V_background, ntheta, nphi, npsi, scale_fac, Q, fem_algorithm, pixel_distance, total_steps, &
rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, status2)

!param_filename=input file name containing initilizing parameters,such as temperature
!temperature=beginning temperature for RMC
!max_move=maximum allowable movement distance
!cutoff_r=cut off distance for different pairs
!used_data_sets=A logical array determine whether electron, neutron, X-ray or FEM data are available or not
!gr_e and r_e=electron G(r) data
!gr_n and r_n=neutron G(r) data
!gr_x and r_x=X-ray G(r) data
!V=FEM intensity variance
!k=scattering vector
!V_err=FEM measurement variance error
!ntheta, nphi, npsi=rotation number for calculating simulated V
!fem_algorithm=which algorithm is used to do variance calculation
!Q=resolution for FEM
!pixel_distance=pixel distance
!status=whether param_filename is opened with success or not

IMPLICIT NONE
!INTEGER, INTENT(IN) :: name_length
CHARACTER (LEN=*), INTENT(IN) :: param_filename  !Assume size array, be careful

REAL, INTENT(OUT) :: temperature
REAL, INTENT(OUT) :: max_move
REAL, INTENT(OUT), DIMENSION(:,:) :: cutoff_r !
LOGICAL, INTENT(OUT), DIMENSION(4) :: used_data_sets
REAL, INTENT(OUT), DIMENSION(4) :: weights
REAL, POINTER, DIMENSION(:) :: gr_e,r_e,gr_e_err
REAL, POINTER, DIMENSION(:) :: gr_n,r_n
REAL, POINTER, DIMENSION(:) :: gr_x,r_x
REAL, POINTER, DIMENSION(:) :: V, k, V_err, V_background
REAL, INTENT(OUT) :: rmin_e, rmax_e
REAL, INTENT(OUT) :: rmin_n, rmax_n
REAL, INTENT(OUT) :: rmin_x, rmax_x
INTEGER, INTENT(OUT) :: ntheta, nphi, npsi
REAL, INTENT(OUT) :: Q
REAL, INTENT(OUT) :: scale_fac
INTEGER, INTENT(OUT) ::fem_algorithm
REAL, INTENT(OUT) :: pixel_distance
INTEGER, INTENT(OUT) ::total_steps
INTEGER, INTENT(OUT) :: status2

!Variable declared in this subroutine, local variables
!comment1=comment in param_filename
CHARACTER (LEN=80) comment1  
CHARACTER (LEN=20) ScatteringFile ! The scattering file name, at most 20 characters
    ! Note: Electron, neutron, and X-ray data all use this name, but in order.
CHARACTER (LEN=20) FemFile ! The FEM data file name, at most 20 characters
INTEGER FileNameLength !The file name length in ScatteringFile or FemFile
INTEGER i, j
real Indicator_End !Indicator_End=-1 means the end of file reaches
INTEGER status1 !Indicate the status of opening file in this subroutine
    !Files for electron, neutron, X-ray scattering and FEM data
LOGICAL File_End !Indicate FileEnd has reached
INTEGER Num_Line !Number of lines in each data file except comment line
REAL, POINTER, DIMENSION(:) :: TempData !Temperature data
INTEGER Stat_Allocate1, Stat_Allocate2, Stat_Allocate3,Stat_Allocate4 !Allocate status, 0 means success

OPEN(20, FILE=param_filename,IOSTAT=status2, STATUS='OLD')
IF(status2 .NE. 0) THEN !Open fails
  PRINT *, 'Cannot open file with name: ', param_filename
  RETURN
ENDIF
READ(20, '(A80)') comment1 !Read comment and it is neglected in the input
!PRINT *, comment1
read(20, *) total_steps
read(20, *) temperature
read(20, *) max_move
read(20, * )
read(20, *) cutoff_r !This is a nelements by nelements matrix

!Fortran read each number in the file in column order

i=0

!********************************************************************
!*******************READ electron scattering file**********************
!********************************************************************
i=i+1
num_line=0

read(20, '(a)') scatteringfile
scatteringfile = adjustl(scatteringfile)
!write(*,*)scatteringfile
read(20, *) weights(i)
read(20, *) rmin_e, rmax_e



if(scatteringfile(1:1) .ne. '#') then !Electron scattering file exists
  used_data_sets(i) = .TRUE.

  FileNameLength=LEN_TRIM(ScatteringFile)
  File_End=.FALSE.

  !Read r_e and gr_e data
  !read r_e first, then gr_e
  !First line is comment
  !At first, count how many data pairs are in the file
  !-1 is chosen as the last point
  !print *, 'scatteringfile(1:filenamelength) 1 is:  ', scatteringfile(1:filenamelength) !debug
  open(30,file=scatteringfile(1:filenamelength),iostat=status1,status='old') 
   if(status1 .eq. 0) then !Open succeeds
     read(30, '(a80)') comment1
     do while( .not. file_end)
        read(30, *) indicator_end
        if(abs(indicator_end+1.0) .le. 1e-6) then !File end is reached
          !PRINT *, 'Indicator_End is: ', Indicator_End
          !PRINT *, 'Electron scattering file end is reached'
          exit !Exit this loop
        else
          num_line = num_line + 1
          !read(30, *) 
        endif
     enddo

     rewind(30) !Go to the beginning of the file

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
       endif !Allocate2 and 3
     else
       print *, 'electron scattering part allocation fail'
       return
     endif
     deallocate(tempdata)
   else
    print *, 'open electron scattering file fails'
   endif
  close(30)
else 
  used_data_sets(i) = .false.
endif !Electron scattering file exists


!*******************READ electron scattering file**********************

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

  !read r_n and gr_n data
  !read r_n first, then gr_n
  !first line is comment
  !at first, count how many data pairs are in the file
  !-1 is chosen as the last point
  open(30,file=scatteringfile(1:filenamelength),iostat=status1,status='old') 
   if(status1 .eq. 0) then !open succeeds
     read(30, '(a80)') comment1
     do while( .not. file_end)
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
       endif !allocate2 and 3
       
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

!*******************read neutron scattering file**********************

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
  !read r_x and gr_x data
  !read r_x first, then gr_x
  !first line is comment
  !at first, count how many data pairs are in the file
  !-1 is chosen as the last point
  !open(30,file=scatteringfile(1:filenamelength),iostat=status1,status='old') - this sucks. it never works. total junk. god damn wasting time. - jwh 052210
  open(30, file="xrd_zr50_kramer_reduced.txt", iostat=status1,status='old')
   if(status1 .eq. 0) then !open succeeds
     read(30, '(a80)') comment1
     do while( .not. file_end)
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
       endif !allocate2 and 3
       
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

!*******************READ X-ray scattering file**********************


!*******************************************************************
!******************Read FEM file************************************
!*******************************************************************

i=i+1
num_line=0

read(20, '(a)') femfile
femfile = adjustl(femfile)
!write(*,*)femfile
read(20, *) weights(i)
read(20, *) scale_fac
read(20, *) nphi, npsi, ntheta
read(20, *) q
!write(*,*)q
read(20, *) fem_algorithm
!write(*,*)fem_algorithm
read(20, *) pixel_distance
!write(*,*)pixel_distance

if(femfile(1:1) .ne. '#') then !fem file exists
  used_data_sets(i) = .true.
  filenamelength=len_trim(femfile)
  file_end=.false.

  !read k ,v, and v_err data
  !read k first, then v, last v_err
  !first line is comment
  !at first, count how many data pairs are in the file
  !-1 is chosen as the last point
  open(30,file=femfile(1:filenamelength),iostat=status1,status='old') 
   if(status1 .eq. 0) then !open succeeds
     read(30, '(a80)') comment1
     do while( .not. file_end)
        read(30, *) indicator_end
        if(abs(indicator_end+1.0) .le. 1e-6) then !file end is reached
          !print *, 'fem file end is reached'
          exit !exit this loop
        else
          num_line = num_line + 1
          !read(30, *) 
        endif
     enddo

     rewind(30) !go to the beginning of the file

     read(30, '(a80)') comment1
     allocate(tempdata(4*num_line),stat=stat_allocate1)
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
    print *, 'open fem file fails'
   endif
  close(30)
else 
  used_data_sets(i) = .false.
endif !FEM file exists

!*******************READ FEM file**********************

CLOSE(20)





!CHARACTER, ALLOCATABLE :: filename


END SUBROUTINE read_inputs

END MODULE ReadInputs

