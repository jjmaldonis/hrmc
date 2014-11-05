!*********************************************************
!**********This module is written by FY on 12/19/2008****
!**********Include three functions***********************
!*********************************************************

!    function chi_square, check_cutoffs and generate_move are included
!

! Changelog
!    unused arguments (r_e etc) removed from chi_sqaure, 12/20/08 pmv
!    subroutine random_move updated. 01/12/08 Jinwoo Hwang
!    function check_cufoff is completed by Feng Yi on 01/16/2009, may need a test program
!    check_cutoffs is modified by Feng Yi on 01/21/2009
!    removed "use ReadInputs" which doesn't seem to be needed; pmv 1/23/09
!    added interface block for function ran2 to enable clean compilation pmv 03/12/09
!    In check_cutoffs, pbc is added to the pair calculation after the hutch is calling. jwh 04/14/2009
!    Added deallocate(atoms) in check_cutoffs, pmv 04/17/09
!    change if condition when deallocate(atoms) in check_cutoffs, FY 04/17/2009
!    gr_e_err is added - jwh 04/25/2009

MODULE rmc_functions

  use model_mod
  implicit none

  interface
     function ran2(idum)
       real :: ran2
       integer :: idum
     end function ran2
  end interface

CONTAINS

    !***************************************************************
    !This function is used to calculate chi square
    FUNCTION chi_square(used_data_sets,weights,gr_e, gr_e_err, gr_n, gr_x, V, V_err,&
            gr_e_sim, gr_n_sim, gr_x_sim, V_sim, scale_fac,&
            rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x, del_r_e, del_r_n, del_r_x, nk, chi2_gr, chi2_vk)

        !used_data_sets=A logical array determine whether electron, neutron, X-ray or FEM data are available
        !gr_e electron G(r) data
        !gr_n neutron G(r) data
        !gr_x X-ray G(r) data
        !V=FEM intensity variance from input file
        !V_err=FEM measurement variance error
        !gr_e_sim=simulated g(r) for electron scattering
        !gr_n_sim=simulated g(r) for neutron scattering
        !gr_x_sim=simulated g(r) for X-ray scattering
        !V_sim=simulated V(k) for variance data
        logical, dimension(4) :: used_data_sets
        real,  dimension(4) :: weights
        real, pointer, dimension(:) :: gr_e, gr_e_err
        real, pointer, dimension(:) :: gr_n
        real, pointer, dimension(:) :: gr_x
        real, pointer, dimension(:) :: v, v_err
        real, pointer, dimension(:) :: gr_e_sim,gr_n_sim,gr_x_sim
        real, pointer, dimension(:) :: v_sim
        real, intent(in) :: scale_fac, rmin_e, rmax_e, rmin_n, rmax_n, rmin_x, rmax_x 
        real, intent (in) :: del_r_e, del_r_n, del_r_x  
        integer, intent(in) ::nk
        real, intent(out) :: chi2_gr, chi2_vk
        integer i, j
        integer nf   !normalization factor   - jwh 04/25/2009
        real chi_square
        real,dimension(4) :: sum1 !record summation of each contribution from diffraction data and fem data
        sum1=0.0
        chi_square=1.0
        
        !Electron diffraction data
        i=1
        if(used_data_sets(i)) then
           !sum1(i)=sum((gr_e_sim-gr_e)**2)*weights(i)
            nf = 0
            do j=int(rmin_e/del_r_e)+1, int(rmax_e/del_r_e)+1    !jwh -032409
                nf = nf + 1
                sum1(i) = sum1(i)+weights(i)*((gr_e_sim(j)-gr_e(j))/gr_e_err(j))**2
            enddo
            sum1(i) = sum1(i)/nf
        endif
        
        !Neutron diffraction data
        i=i+1
        if(used_data_sets(i)) then
           !sum1(i)=sum((gr_n_sim-gr_n)**2)*weights(i)
            nf = 0
            do j=int(rmin_n/del_r_n)+1, int(rmax_n/del_r_n)+1    !jwh -032409
                nf = nf +1
                sum1(i) = sum1(i)+weights(i)*(gr_n_sim(j)-gr_n(j))**2
            enddo
            sum1(i) = sum1(i)/nf
        endif
        
        !X-ray diffraction data
        i=i+1
        if(used_data_sets(i)) then
           !sum1(i)=sum((gr_x_sim-gr_x)**2)*weights(i)
            nf = 0
            do j=int(rmin_x/del_r_x)+1, int(rmax_x/del_r_x)+1    !jwh -032409
                nf = nf + 1
                sum1(i) = sum1(i)+weights(i)*(gr_x_sim(j)-gr_x(j))**2
            enddo
            sum1(i) = sum1(i)/nf
        endif
        
        !FEM data
        i=i+1
        !write(*,*) "DEBUG: sum1(i) starts at:", sum1(i)
        if(used_data_sets(i)) then
           !sum1(i)=sum(((V_sim-V*scale_fac)/(scale_fac*V_err))**2)*weights(i)
            nf = 0
            do j=1,nk
                nf = nf + 1
                !write(*,*) "DEBUG: weights(4)=", weights(4)
                !write(*,*) "DEBUG: v(j)*scale_fac-v_sim(j)=", v(j)*scale_fac-v_sim(j)
                !write(*,*) "DEBUG: scale_fac*V_err(j)=", scale_fac*V_err(j)
                !write(*,*) "DEBUG: Total addition is:", weights(4)*((v(j)*scale_fac-v_sim(j))/(scale_fac*V_err(j)))**2
                sum1(i) = sum1(i)+weights(4)*((v(j)*scale_fac-v_sim(j))/(scale_fac*V_err(j)))**2
                !write(*,*) "DEBUG: New sum1(i)=", sum1(i)
            enddo
            sum1(i) = sum1(i)/nf
            !write(*,*) "DEBUG: Final sum1(i) = sum1(i)/nf=", sum1(i)
            !write(*,*) "DEBUG: where nf=", nf
        endif

        chi2_gr = sum1(3)
        chi2_vk = sum1(4)
        chi_square=sum(sum1)
    end function chi_square


    subroutine update_scale_factor(sf, sf_initial, vsim, vas)
    ! Compute chi2(vsim*x, vas) for variable x to minimize chi2.
    ! (Note, I mean "variable" as in "changing", not the noun.)
    ! Set sf = sf/x.
    ! This isn't going to work because calling this multiple times will
    ! eventually monotonically change sf because vsim vs. vas will always be a
    ! little smaller or a little bigger (likely) and therefore sf will get worse
    ! and worse after teh first call. I need another variable I think.
        real, intent(out) :: sf
        real, intent(in) :: sf_initial
        ! vsim is the simulated vk as computed by the fast intensity algorithm. 
        ! vas is the simulated vk as computed by autoslice.
        real, intent(in), pointer, dimension(:) :: vsim, vas 
        real :: x! This is the parameter we will use to fit vsim to vas.
        integer :: i

        ! Below is the numerical solved solution to mimize
        ! chi2 = summation( (vsim(i)*x - vas(i))**2/vas(i), i=1, i=size(vsim) )
        ! by varying x. x will then be the weighting factor that gets vsim as
        ! close to vas as possible. I.e. vsim*x ~= vas.
        x = 0
        do i=1, size(vsim)
            x = x + vas(i)/vsim(i)
        enddo

        sf = sf_initial/x ! We need to divide because vk_exp ~= vsim/sf
        write(*,*) "Scale factor updated. Old,New =", sf_initial, sf
        write(*,*) "Conversion factor x =", x
    end subroutine update_scale_factor



!SUBROUTINE Calc_Type_Order(m, Type_Order) must be called first, then 
!this function can be called
!The check_cufoffs function return .TRUE. if 
!****no two atoms are within cutoff distance
!Otherwise return .FALSE.

    function check_cutoffs(m,cutoff_r,moved_atom)
        logical check_cutoffs
        real,dimension(:,:) ::cutoff_r
        integer moved_atom
        type(model) m
        integer, dimension(:), pointer :: atoms
        integer  nlist
        integer istat
        real radius, temp_x, temp_y, temp_z
        integer i,j
        integer num1, num2
        real dist_pair

        !find the maximum cut-off distance
        radius=maxval(maxval(cutoff_r, 1),1)

        ! Find all atoms within radius radius of moved_atom and put them 
        ! in the list atoms. Also sets nlist == size(atoms)+1.
        call hutch_list_3D(m, m%xx%ind(moved_atom),m%yy%ind(moved_atom),m%zz%ind(moved_atom), radius, atoms, istat, nlist)
        if (istat .eq. 1) then
            print *, 'memory allocation fails!'
            return
        else if(istat .eq. -1) then
            print *, 'no atom is found'
            check_cutoffs = .true. 
            return
        endif

        !begin to calculate pair distance
        !note, there are (nlist-1) atoms in total in 'atoms'

        !first, determine the type of moved_atom 
        do i=1, m%nelements
            if(m%znum%ind(moved_atom) .eq. m%atom_type(i)) then
                num1 = i
                exit
            endif
        enddo
        
        do i=1, (nlist-1)
            !do not count the pair distance to itself
            ! The below if statement should be irrelevant. moved_atom isn't put
            ! in atoms by hutch_list_3d as far as i know.
            if (atoms(i) .ne. moved_atom) then
                do j=1, m%nelements
                    if(m%atom_type(j) .eq. m%znum%ind(atoms(i))) then
                        num2 = j
                        exit
                    endif
                enddo !j=1, m%nelements
              
                !endif  !032409 - jwh
               
                !calculate the atomic distance
                !compare with cutoff_r
                temp_x = abs(m%xx%ind(moved_atom) - m%xx%ind(atoms(i)))  !pbc added - jwh 04/14/2009
                temp_y = abs(m%yy%ind(moved_atom) - m%yy%ind(atoms(i)))
                temp_z = abs(m%zz%ind(moved_atom) - m%zz%ind(atoms(i)))

                !write(*,*)"abs=", temp_x, temp_y, temp_z
                
                temp_x = temp_x - m%lx*anint(temp_x/m%lx)
                temp_y = temp_y - m%ly*anint(temp_y/m%ly)
                temp_z = temp_z - m%lz*anint(temp_z/m%lz)
                  
                dist_pair = temp_x**2 + temp_y**2 + temp_z**2
                dist_pair = sqrt(dist_pair)
                   
                if (dist_pair  .lt. cutoff_r(num1, num2)) then
                    !write(*,*)dist_pair  , cutoff_r(num1, num2) !debug - jwh 032409
                    check_cutoffs=.false.
                    exit
                else
                    check_cutoffs = .true.  !jwh 032409
                endif
            endif !032409 - jwh
        enddo !i=1, (nlist-1)

        !if (associated(atoms)) deallocate(atoms) ! pmv 4/17/09
        if (nlist .gt. 1) deallocate(atoms) ! fy 4/17/09
    
    end function check_cutoffs
  


  !subroutine random_move---------------------------------------------------------------
  !Generates a random move of one atom.
  !Written for testing purpose here.
  !The mc move will be done some place else later, but it should still update the variables xx_cur, and xx_new (and etc).
  subroutine random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, alpha)
    use mpi  
    type(model), intent(inout) :: m
    integer iseed, w
    real alpha, aa, bb, cc
    real, intent(out) :: xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new    !Cur and new positions of atoms
    real :: rand1, rand2, rand3, rand4


    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    rand3 = ran2(iseed)
    rand4 = ran2(iseed)

    w = int(m%natoms*rand1)+1

    !write(*,*)myid, rand1, rand2, rand3, rand4

    xx_cur = m%xx%ind(w)            !Cur positions of the atom before random move
    yy_cur = m%yy%ind(w)
    zz_cur = m%zz%ind(w)
    
    aa = alpha*(rand2 - 0.5)
    bb = alpha*(rand3 - 0.5)
    cc = alpha*(rand4 - 0.5)
    
    m%xx%ind(w) = m%xx%ind(w) + aa 
    m%yy%ind(w) = m%yy%ind(w) + bb 
    m%zz%ind(w) = m%zz%ind(w) + cc
    
    if(m%xx%ind(w)>m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)-m%lx       !pbc 
    if(m%yy%ind(w)>m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)-m%ly
    if(m%zz%ind(w)>m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)-m%lz
    if(m%xx%ind(w)<-m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)+m%lx
    if(m%yy%ind(w)<-m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)+m%ly
    if(m%zz%ind(w)<-m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)+m%lz
    
    xx_new=m%xx%ind(w)              !new positions of the atom after random move
    yy_new=m%yy%ind(w)
    zz_new=m%zz%ind(w)
  end subroutine random_move


END MODULE rmc_functions
