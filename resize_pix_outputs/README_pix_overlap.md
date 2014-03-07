In the pixel overlap experiments I increased the pixel overlap to increase the number of times each atom was sampled.

Below I changed lane 334 and added lines 350-351. On line 334 I added the *2 at the end in increased it from 2 to 3 in steps of 1. This makes femsim include the int_i from each pixel and since the pixels are duplicated each atom is being included more. It segfauled on *4 so I stopped there.

 321     subroutine init_pix(m, res, istat, square_pixel)
 322         type(model), intent(in) :: m
 323         real, intent(in) :: res
 324         integer, intent(out) :: istat
 325         logical, intent(in) :: square_pixel
 326         integer :: i, j, k
 327
 328         if(square_pixel) then
 329             pa%phys_diam = res * sqrt(2.0)
 330         else
 331             pa%phys_diam = res
 332         endif
 333         pa%npix_1D = floor( m%lx / pa%phys_diam )
 334         pa%npix = pa%npix_1D**2 *2 !REMOVE 2*! TODO
 335
 336         pa%dr = m%lx/pa%npix_1D - pa%phys_diam
 337
 338         allocate(pa%pix(pa%npix, 2), stat=istat)
 339         if (istat /= 0) then
 340             write (*,*) 'Cannot allocate pixel position array.'
 341             return
 342         endif
 343
 344         k=1
 345         do i=1, pa%npix_1D
 346             do j=1, pa%npix_1D
 347                 pa%pix(k,1) = -m%lx/2.0 + (pa%phys_diam+pa%dr)/2.0 + (pa%phys_diam+pa%dr)*(i-1)
 348                 pa%pix(k,2) = -m%ly/2.0 + (pa%phys_diam+pa%dr)/2.0 + (pa%phys_diam+pa%dr)*(j-1)
 349                 pa%pix(k+pa%npix_1D,1) = -m%lx/2.0 + (pa%phys_diam+pa%dr)/2.0 + (pa%phys_diam+pa%dr)*(i-1) ! REMOVE ME
 350                 pa%pix(k+pa%npix_1D,2) = -m%ly/2.0 + (pa%phys_diam+pa%dr)/2.0 + (pa%phys_diam+pa%dr)*(j-1) ! REMOVE ME
 351                 k = k + 1
 352             enddo
 353         enddo

