



module mutations

use model_mod
use rmc_global
use rmc_functions

  !interface
  !   function ran2(idum)
  !     real :: ran2
  !     integer :: idum
  !   end function ran2
  !end interface

contains

subroutine cluster_move(m,atom1,atom2,atoms1,atoms2,cnatoms1,cnatoms2)
    implicit none
    type(model), intent(inout) :: m
    integer :: i,j,k, atom1, atom2, cnatoms1, cnatoms2, istat
    integer, dimension(:), pointer :: atoms1, atoms2
    integer :: iseed
    real :: rand1, rand2, rand3, rand4
    real :: cx1, cy1, cz1, cx2, cy2, cz2 ! Cluster centers of volume
    real :: dx, dy, dz ! Displacement vector
    real :: mind, r, r2, rx1, ry1, rz1, rx2, ry2, rz2

    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    atom1 = int(m%natoms*rand1)+1
    atom2 = int(m%natoms*rand2)+1
    !write(*,*) "Moving clusters:", atom1, atom2

    call hutch_list_3D(m, m%xx%ind(atom1), m%yy%ind(atom1), m%zz%ind(atom1), 3.5, atoms1, istat, cnatoms1)
    cx1 = m%xx%ind(atom1)
    cy1 = m%yy%ind(atom1)
    cz1 = m%zz%ind(atom1)
    do i=1, cnatoms1-1
        cx1 = cx1 + m%xx%ind(atoms1(i))
        cy1 = cy1 + m%yy%ind(atoms1(i))
        cz1 = cz1 + m%zz%ind(atoms1(i))
    enddo
    cx1 = cx1/float(cnatoms1)
    cy1 = cy1/float(cnatoms1)
    cz1 = cz1/float(cnatoms1)

    call hutch_list_3D(m, m%xx%ind(atom2), m%yy%ind(atom2), m%zz%ind(atom2), 3.5, atoms2, istat, cnatoms2)
    cx2 = m%xx%ind(atom2)
    cy2 = m%yy%ind(atom2)
    cz2 = m%zz%ind(atom2)
    do i=1, cnatoms2-1
        cx2 = cx2 + m%xx%ind(atoms2(i))
        cy2 = cy2 + m%yy%ind(atoms2(i))
        cz2 = cz2 + m%zz%ind(atoms2(i))
    enddo
    cx2 = cx2/float(cnatoms2)
    cy2 = cy2/float(cnatoms2)
    cz2 = cz2/float(cnatoms2)

    dx = cx2 - cx1
    dy = cy2 - cy1
    dz = cz2 - cz1
    !write(*,*)  "dx,dy,dz", dx,dy,dz
    !write(*,*) "boxsizes", -m%lx/2.0, m%lx/2

    m%xx%ind(atom1) = m%xx%ind(atom1) + dx
    m%yy%ind(atom1) = m%yy%ind(atom1) + dy
    m%zz%ind(atom1) = m%zz%ind(atom1) + dz
    if(m%xx%ind(atom1) .gt. m%lx*0.5) m%xx%ind(atom1)=m%xx%ind(atom1)-m%lx
    if(m%yy%ind(atom1) .gt. m%ly*0.5) m%yy%ind(atom1)=m%yy%ind(atom1)-m%ly
    if(m%zz%ind(atom1) .gt. m%lz*0.5) m%zz%ind(atom1)=m%zz%ind(atom1)-m%lz
    if(m%xx%ind(atom1) .lt. -m%lx*0.5) m%xx%ind(atom1)=m%xx%ind(atom1)+m%lx
    if(m%yy%ind(atom1) .lt. -m%ly*0.5) m%yy%ind(atom1)=m%yy%ind(atom1)+m%ly
    if(m%zz%ind(atom1) .lt. -m%lz*0.5) m%zz%ind(atom1)=m%zz%ind(atom1)+m%lz
    !write(*,*) m%xx%ind(atom1), m%yy%ind(atom1), m%zz%ind(atom1)
    call hutch_move_atom(m,atom1, m%xx%ind(atom1), m%yy%ind(atom1), m%zz%ind(atom1))
    !write(*,*) m%ha%atom_hutch(atom1,1), m%ha%atom_hutch(atom1,2), m%ha%atom_hutch(atom1,3)
    do i=1, cnatoms1-1
        m%xx%ind(atoms1(i)) = m%xx%ind(atoms1(i)) + dx
        m%yy%ind(atoms1(i)) = m%yy%ind(atoms1(i)) + dy
        m%zz%ind(atoms1(i)) = m%zz%ind(atoms1(i)) + dz
        if(m%xx%ind(atoms1(i)) .gt. m%lx*0.5) m%xx%ind(atoms1(i))=m%xx%ind(atoms1(i))-m%lx
        if(m%yy%ind(atoms1(i)) .gt. m%ly*0.5) m%yy%ind(atoms1(i))=m%yy%ind(atoms1(i))-m%ly
        if(m%zz%ind(atoms1(i)) .gt. m%lz*0.5) m%zz%ind(atoms1(i))=m%zz%ind(atoms1(i))-m%lz
        if(m%xx%ind(atoms1(i)) .lt. -m%lx*0.5) m%xx%ind(atoms1(i))=m%xx%ind(atoms1(i))+m%lx
        if(m%yy%ind(atoms1(i)) .lt. -m%ly*0.5) m%yy%ind(atoms1(i))=m%yy%ind(atoms1(i))+m%ly
        if(m%zz%ind(atoms1(i)) .lt. -m%lz*0.5) m%zz%ind(atoms1(i))=m%zz%ind(atoms1(i))+m%lz
        !write(*,*) m%xx%ind(atoms1(i)), m%yy%ind(atoms1(i)), m%zz%ind(atoms1(i))
        call hutch_move_atom(m,atom1, m%xx%ind(atoms1(i)), m%yy%ind(atoms1(i)), m%zz%ind(atoms1(i)))
        !write(*,*) m%ha%atom_hutch(atoms1(i),1), m%ha%atom_hutch(atoms1(i),2), m%ha%atom_hutch(atoms1(i),3)
    enddo

    m%xx%ind(atom2) = m%xx%ind(atom2) - dx
    m%yy%ind(atom2) = m%yy%ind(atom2) - dy
    m%zz%ind(atom2) = m%zz%ind(atom2) - dz
    if(m%xx%ind(atom2) .gt. m%lx*0.5) m%xx%ind(atom2)=m%xx%ind(atom2)-m%lx
    if(m%yy%ind(atom2) .gt. m%ly*0.5) m%yy%ind(atom2)=m%yy%ind(atom2)-m%ly
    if(m%zz%ind(atom2) .gt. m%lz*0.5) m%zz%ind(atom2)=m%zz%ind(atom2)-m%lz
    if(m%xx%ind(atom2) .lt. -m%lx*0.5) m%xx%ind(atom2)=m%xx%ind(atom2)+m%lx
    if(m%yy%ind(atom2) .lt. -m%ly*0.5) m%yy%ind(atom2)=m%yy%ind(atom2)+m%ly
    if(m%zz%ind(atom2) .lt. -m%lz*0.5) m%zz%ind(atom2)=m%zz%ind(atom2)+m%lz
    !write(*,*) m%xx%ind(atom2), m%yy%ind(atom2), m%zz%ind(atom2)
    call hutch_move_atom(m,atom2, m%xx%ind(atom2), m%yy%ind(atom2), m%zz%ind(atom2))
    !write(*,*) m%ha%atom_hutch(atom2,1), m%ha%atom_hutch(atom2,2), m%ha%atom_hutch(atom2,3)
    do i=1, cnatoms2-1
        m%xx%ind(atoms2(i)) = m%xx%ind(atoms2(i)) - dx
        m%yy%ind(atoms2(i)) = m%yy%ind(atoms2(i)) - dy
        m%zz%ind(atoms2(i)) = m%zz%ind(atoms2(i)) - dz
        if(m%xx%ind(atoms2(i)) .gt. m%lx*0.5) m%xx%ind(atoms2(i))=m%xx%ind(atoms2(i))-m%lx
        if(m%yy%ind(atoms2(i)) .gt. m%ly*0.5) m%yy%ind(atoms2(i))=m%yy%ind(atoms2(i))-m%ly
        if(m%zz%ind(atoms2(i)) .gt. m%lz*0.5) m%zz%ind(atoms2(i))=m%zz%ind(atoms2(i))-m%lz
        if(m%xx%ind(atoms2(i)) .lt. -m%lx*0.5) m%xx%ind(atoms2(i))=m%xx%ind(atoms2(i))+m%lx
        if(m%yy%ind(atoms2(i)) .lt. -m%ly*0.5) m%yy%ind(atoms2(i))=m%yy%ind(atoms2(i))+m%ly
        if(m%zz%ind(atoms2(i)) .lt. -m%lz*0.5) m%zz%ind(atoms2(i))=m%zz%ind(atoms2(i))+m%lz
        !write(*,*) m%xx%ind(atoms2(i)), m%yy%ind(atoms2(i)), m%zz%ind(atoms2(i))
        call hutch_move_atom(m,atom2, m%xx%ind(atoms2(i)), m%yy%ind(atoms2(i)), m%zz%ind(atoms2(i)))
        !write(*,*) m%ha%atom_hutch(atoms2(i),1), m%ha%atom_hutch(atoms2(i),2), m%ha%atom_hutch(atoms2(i),3)
    enddo

    !mind = 10000.0
    !do i=1, m%natoms
    !    do j=1, m%natoms
    !        if(i .ne. j) then
    !            r = sqrt( (m%xx%ind(i)-m%xx%ind(j))**2 + (m%yy%ind(i)-m%yy%ind(j))**2 + (m%zz%ind(i)-m%zz%ind(j))**2)
    !            if(r .lt. mind) mind = r
    !        endif
    !    enddo
    !enddo
    !write(*,*) "MIN DIST BETWEEN ATOMS IS", mind



end subroutine cluster_move

subroutine z_swap(m,atom1,atom2)
    implicit none
    type(model), intent(inout) :: m
    integer :: i,j,k, atom1, atom2
    integer :: iseed, temp
    real :: rand1, rand2, rand3, rand4

    iseed = 791315

    rand1 = ran2(iseed)
    rand2 = ran2(iseed)
    atom1 = int(m%natoms*rand1)+1
    atom2 = int(m%natoms*rand2)+1

    do while ( m%znum%ind(atom1) .eq. m%znum%ind(atom2) )
        write(*,*) "Had to retry", m%znum%ind(atom1), m%znum%ind(atom2)
        rand1 = ran2(iseed)
        rand2 = ran2(iseed)
        atom1 = int(m%natoms*rand1)+1
        atom2 = int(m%natoms*rand2)+1
    enddo

    temp = m%znum%ind(atom1)
    m%znum%ind(atom1) = m%znum%ind(atom2)
    m%znum%ind(atom2) = temp
    
end subroutine z_swap


end module
