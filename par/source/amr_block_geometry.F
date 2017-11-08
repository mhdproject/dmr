!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

         subroutine amr_block_geometry(lb,pe)

!-----------------------------------------------------------
!
! This routine computes cell volumes, areas and edge lengths
! for various grid geometries, for the specified local block lb.

!
!
! Cartesian :
!     coord 1      x
!     coord 2      y
!     coord 3      z
!
! Cylindrical :
!     coord 1      r 
!     coord 2      theta
!     coord 3      z
!
! Spherical :
!     coord 1      r
!     coord 2      theta
!     coord 3      phi           (azimuthal)
!
! Polar (2D) :
!     coord 1      r 
!     coord 2      theta
!     coord 3      z             (has only 1 grid cell in this direction)
!
!
! Written : Peter MacNeice      December 2001

!-----------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit none

#include "amr_shmem.fh"

      integer, intent(in) :: lb,pe

!-----------------------------------------------------------
#ifdef CURVILINEAR

! Local arrays and variables

         real :: cell_face_coord1w(ilw1:iuw1+1)
         real :: cell_face_coord2w(jlw1:juw1+k2d)
         real :: cell_face_coord3w(klw1:kuw1+k3d)
      

         real :: del
         real :: cell_vol_1  ,cell_vol_2  ,cell_vol_3
         real :: cell_area1_1,cell_area1_2,cell_area1_3
         real :: cell_area2_1,cell_area2_2,cell_area2_3
         real :: cell_area3_1,cell_area3_2,cell_area3_3
         real :: cell_leng1_1,cell_leng1_2,cell_leng1_3
         real :: cell_leng2_1,cell_leng2_2,cell_leng2_3
         real :: cell_leng3_1,cell_leng3_2,cell_leng3_3
 
         real,save :: cbnd_box(2,3)

         integer :: ierr_trap
         integer :: mype
         integer,external :: shmem_my_pe

         integer :: lb0,pe0,iloc

         logical :: lfound

!-----------------------

         mype = shmem_my_pe()

         ierr_trap = 0
#ifdef CARTESIAN
         ierr_trap = ierr_trap + 1
#endif /* CARTESIAN */
#ifdef SPHERICAL
         ierr_trap = ierr_trap + 1
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
         ierr_trap = ierr_trap + 1
#endif /* CYLINDRICAL */
#ifdef POLAR
         ierr_trap = ierr_trap + 1
         if(ndim.ne.2) ierr_trap = ierr_trap + 1
#endif /* POLAR */
         if(ierr_trap.gt.1) then
           write(*,*) 'Paramesh ERROR : amr_block_geometry. ',
     .           'Inconsistent choice of curvilinear coord.'
           call amr_abort()
         endif

!-----------------------
!
! If using SHMEM communication library then, if pe is not
! the local pe, get the bounding box info from the remote
! location. If using MPI then you can assume that this information
! has already been stored locally.
#ifdef SHMEM_USED
        call shmem_real_get(cbnd_box(1,1),bnd_box(1,1,lb),6,pe)
        lb0 = lb
        pe0 = pe
#else
        lfound = .false.
        if(pe.eq.mype) then
          if(lb.le.lnblocks) then
            lfound = .true.
            lb0 = lb
            pe0 = mype
          elseif(lb.ge.strt_buffer.and.lb.le.last_buffer) then
            lfound = .true.
            lb0 = lb
            pe0 = mype
          endif
        else
          iloc = strt_buffer
          do while( (iloc.le.last_buffer) .and.
     .              (.not.lfound)  )
            if(laddress(1,iloc).eq.lb.and.
     .         laddress(2,iloc).eq.pe ) lfound = .true.
            if(.not.lfound) iloc = iloc + 1
          enddo
          if(lfound) then
            lb0 = iloc
            pe0 = mype
          endif
        endif

        if(.not.lfound) then
          write(*,*) 'amr_block_geometry ERROR : blk ',
     .        lb,pe,' not found on pe ',mype,
     .        ' strt_buffer:last_buffer ',strt_buffer,last_buffer,
     .        ' laddress ',laddress(:,strt_buffer:last_buffer)
          call amr_abort()
        endif
        cbnd_box(:,:) = bnd_box(:,:,lb0)
#endif /* SHMEM_USED */


!-----------------------
! compute coords of cell interfaces


! for first coordinate direction
         del = (cbnd_box(2,1)-cbnd_box(1,1))/real(nxb)
         do i = il_bnd1,iu_bnd1+1
           cell_face_coord1(i) = cbnd_box(1,1) + del*real(i-1-nguard)
         enddo
         do i = ilw1,iuw1+1
           cell_face_coord1w(i) = cbnd_box(1,1) 
     .                            + del*real(i-1-nguard_work)
         enddo


! for second coordinate direction
         cell_face_coord2 = 0.
         cell_face_coord2w = 0.
         if(ndim.ge.2) then
         del = (cbnd_box(2,2)-cbnd_box(1,2))/real(nyb)
         do j = jl_bnd1,ju_bnd1+1
           cell_face_coord2(j) = cbnd_box(1,2) + del*real(j-1-nguard)
         enddo
         do j = jlw1,juw1+1
           cell_face_coord2w(j) = cbnd_box(1,2) 
     .                            + del*real(j-1-nguard_work)
         enddo
         endif

! for third coordinate direction
         cell_face_coord3 = 0.
         cell_face_coord3w = 0.
         if(ndim.eq.3) then
         del = (cbnd_box(2,3)-cbnd_box(1,3))/real(nzb)
         do k = kl_bnd1,ku_bnd1+1
           cell_face_coord3(k) = cbnd_box(1,3) + del*real(k-1-nguard)
         enddo
         do k = klw1,kuw1+1
           cell_face_coord3w(k) = cbnd_box(1,3) 
     .                            + del*real(k-1-nguard_work)
         enddo
         endif

!-----------------------

! Apply any user specified coordinate transformation
         call user_coord_transfm(lb0,pe0)


!-----------------------
! compute cell volumes

! Note the style used here to compute cell_vol. We
! specify dependence of cell_vol on coord 1 in cell_vol_1,
! specify dependence of cell_vol on coord 2 in cell_vol_2,
! specify dependence of cell_vol on coord 3 in cell_vol_3.
! This style is used throughout this routine.


! first cell volumes for use with UNK data structure
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1


           cell_vol_2 = 1.
           cell_vol_3 = 1.
#ifdef CARTESIAN
           cell_vol_1 =  cell_face_coord1(i+1) - 
     .                   cell_face_coord1(i)
           if(ndim.ge.2)
     .       cell_vol_2 =  cell_face_coord2(j+k2d) -
     .                     cell_face_coord2(j)
           if(ndim.eq.3)
     .       cell_vol_3 =  cell_face_coord3(k+k3d) -
     .                     cell_face_coord3(k)
           if(l2p5d.eq.1) cell_vol_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_vol_1 = ( cell_face_coord1(i+1)**3 - 
     .                    cell_face_coord1(i)**3 ) /3.
           if(ndim.ge.2)
     .       cell_vol_2 = cos( cell_face_coord2(j)   ) -
     .                    cos( cell_face_coord2(j+1) )
           if(ndim.eq.3)
     .        cell_vol_3 = cell_face_coord3(k+k3d) - 
     .                     cell_face_coord3(k)
           if(l2p5d.eq.1) cell_vol_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
     
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_vol_1 = ( cell_face_coord1(i+1)**2 - 
     .                    cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2)
     .       cell_vol_2 =  cell_face_coord2(j+1) -
     .                     cell_face_coord2(j)
           if(ndim.eq.3)
     .       cell_vol_3 =  cell_face_coord2(k+k3d) -
     .                     cell_face_coord3(k)
           if(l2p5d.eq.1) cell_vol_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_vol_1 = ( cell_face_coord1(i+1)**2 - 
     .                    cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2)
           cell_vol_2 =  cell_face_coord2(j+1) -
     .                   cell_face_coord2(j)
#endif /* CYLINDRICAL */

           cell_vol(i,j,k) = cell_vol_1 * cell_vol_2 
     .                                  * cell_vol_3

         enddo
         enddo
         enddo



! now cell volumes for use with WORK data structure
         do k = klw1,kuw1
         do j = jlw1,juw1
         do i = ilw1,iuw1

           cell_vol_2 = 1.
           cell_vol_3 = 1.
#ifdef CARTESIAN
           cell_vol_1 =  cell_face_coord1w(i+1) -
     .                   cell_face_coord1w(i)
           if(ndim.ge.2)
     .       cell_vol_2 =  cell_face_coord2w(j+k2d) -
     .                     cell_face_coord2w(j)
           if(ndim.eq.3)
     .       cell_vol_3 =  cell_face_coord3w(k+k3d) -
     .                     cell_face_coord3w(k)
           if(l2p5d.eq.1) cell_vol_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_vol_1 = ( cell_face_coord1w(i+1)**3 -
     .                    cell_face_coord1w(i)**3 ) /3.
           if(ndim.ge.2)
     .       cell_vol_2 = cos( cell_face_coord2w(j)   ) -
     .                    cos( cell_face_coord2w(j+1) )
           if(ndim.eq.3)
     .       cell_vol_3 = cell_face_coord3w(k+k3d) -
     .                    cell_face_coord3w(k)
           if(l2p5d.eq.1) cell_vol_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_vol_1 = ( cell_face_coord1w(i+1)**2 -
     .                    cell_face_coord1w(i)**2 )*.5
           if(ndim.ge.2)
     .       cell_vol_2 =  cell_face_coord2w(j+1) -
     .                     cell_face_coord2w(j)
           if(ndim.eq.3)
     .       cell_vol_3 =  cell_face_coord2w(k+k3d) -
     .                     cell_face_coord3w(k)
           if(l2p5d.eq.1) cell_vol_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_vol_1 = ( cell_face_coord1w(i+1)**2 -
     .                    cell_face_coord1w(i)**2 )*.5
           if(ndim.ge.2)
     .       cell_vol_2 =  cell_face_coord2w(j+1) -
     .                     cell_face_coord2w(j)
#endif /* POLAR */

           cell_vol_w(i,j,k) = cell_vol_1 * cell_vol_2
     .                                  * cell_vol_3

         enddo
         enddo
         enddo



!-----------------------
! Compute cell face areas


! compute cell area of faces perpendicular to first coord axis
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1+1

           cell_area1_2 = 1.
           cell_area1_3 = 1.
#ifdef CARTESIAN
           cell_area1_1 =  1.
           if(ndim.ge.2)
     .       cell_area1_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
           if(ndim.eq.3)
     .       cell_area1_3 =  cell_face_coord3(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_area1_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_area1_1 = cell_face_coord1(i)**2
           if(ndim.ge.2)
     .       cell_area1_2 = cos( cell_face_coord2(j)   ) -
     .                      cos( cell_face_coord2(j+1) )
           if(ndim.eq.3)
     .       cell_area1_3 = cell_face_coord3(k+k3d) -
     .                      cell_face_coord3(k)
           if(l2p5d.eq.1) cell_area1_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_area1_1 =  cell_face_coord1(i) 
           if(ndim.ge.2)
     .       cell_area1_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
           if(ndim.eq.3)
     .       cell_area1_3 =  cell_face_coord2(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_area1_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_area1_1 =  cell_face_coord1(i) 
           if(ndim.ge.2)
     .       cell_area1_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* POLAR */

           cell_area1(i,j,k) = cell_area1_1 * cell_area1_2 
     .                                      * cell_area1_3

         enddo
         enddo
         enddo



! compute cell area of faces perpendicular to second coord axis
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1+k2d
         do i = il_bnd1,iu_bnd1

           cell_area2_2 = 1.
           cell_area2_3 = 1.
#ifdef CARTESIAN
           cell_area2_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
           if(ndim.eq.3)
     .       cell_area2_3 =  cell_face_coord3(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_area2_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_area2_1 = (cell_face_coord1(i+1)-cell_face_coord1(i))
     .            *(cell_face_coord1(i)+cell_face_coord1(i+1))*.5
           cell_area2_2 = sin( cell_face_coord2(j) )
           if(ndim.eq.3)
     .       cell_area2_3 = cell_face_coord3(k+k3d) -
     .                      cell_face_coord3(k)
           if(l2p5d.eq.1) cell_area2_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_area2_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
           if(ndim.eq.3)
     .       cell_area2_3 =  cell_face_coord2(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_area2_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_area2_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
#endif /* POLAR */

           cell_area2(i,j,k) = cell_area2_1 * cell_area2_2 
     .                                      * cell_area2_3

         enddo
         enddo
         enddo



! compute cell area of faces perpendicular to third coord axis
         do k = kl_bnd1,ku_bnd1+k3d
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1

           cell_area3_2 =  1.
           cell_area3_3 =  1.
#ifdef CARTESIAN
           cell_area3_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
           if(ndim.ge.2)
     .       cell_area3_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_area3_1 = (cell_face_coord1(i+1)-cell_face_coord1(i))
     .            *(cell_face_coord1(i)+cell_face_coord1(i+1))*.5
           if(ndim.ge.2)
     .       cell_area3_2 = cell_face_coord2(j+1) -
     .                      cell_face_coord2(j)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_area3_1 = ( cell_face_coord1(i+1)**2 - 
     .                    cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2)
     .       cell_area3_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_area3_1 = ( cell_face_coord1(i+1)**2 - 
     .                    cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2)
     .       cell_area3_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* POLAR */

           cell_area3(i,j,k) = cell_area3_1 * cell_area3_2 
     .                                      * cell_area3_3

         enddo
         enddo
         enddo

!-----------------------
! Compute cell edge lengths


! compute edge length in direction of first coord axis
         do k = kl_bnd1,ku_bnd1+k3d
         do j = jl_bnd1,ju_bnd1+k2d
         do i = il_bnd1,iu_bnd1

           cell_leng1_2 =  1.
           cell_leng1_3 =  1.
#ifdef CARTESIAN
           cell_leng1_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_leng1_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_leng1_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_leng1_1 =  cell_face_coord1(i+1) -
     .                     cell_face_coord1(i)
#endif /* POLAR */

           cell_leng1(i,j,k) = cell_leng1_1 * cell_leng1_2 
     .                                      * cell_leng1_3

         enddo
         enddo
         enddo


! compute edge length in direction of second coord axis
         do k = kl_bnd1,ku_bnd1+k3d
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1+1

           cell_leng2_1 =  1.
           cell_leng2_2 =  1.
           cell_leng2_3 =  1.
#ifdef CARTESIAN
           if(ndim.ge.2)
     .       cell_leng2_2 =  cell_face_coord2(j+k2d) -
     .                       cell_face_coord2(j)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_leng2_1 =  cell_face_coord1(i)
           if(ndim.ge.2)
     .       cell_leng2_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           cell_leng2_1 =  cell_face_coord1(i)
           if(ndim.ge.2)
     .       cell_leng2_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* CYLINDRICAL */
#ifdef POLAR
           cell_leng2_1 =  cell_face_coord1(i)
           if(ndim.ge.2)
     .       cell_leng2_2 =  cell_face_coord2(j+1) -
     .                       cell_face_coord2(j)
#endif /* POLAR */

           cell_leng2(i,j,k) = cell_leng2_1 * cell_leng2_2
     .                                      * cell_leng2_3

         enddo
         enddo
         enddo


! compute edge length in direction of third coord axis
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1+k2d
         do i = il_bnd1,iu_bnd1+1

           cell_leng3_1 =  1.
           cell_leng3_2 =  1.
           cell_leng3_3 =  1.
#ifdef CARTESIAN
           if(ndim.eq.3)
     .       cell_leng3_3 =  cell_face_coord3(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_leng3_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CARTESIAN */
#ifdef SPHERICAL
           cell_leng3_1 =  cell_face_coord1(i)
           cell_leng3_2 =  sin( cell_face_coord2(j) )
           if(ndim.eq.3)
     .       cell_leng3_3 =  cell_face_coord3(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_leng3_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* SPHERICAL */
#ifdef CYLINDRICAL
           if(ndim.eq.3)
     .       cell_leng3_3 =  cell_face_coord3(k+k3d) -
     .                       cell_face_coord3(k)
           if(l2p5d.eq.1) cell_leng3_3 = 
     .                    cbnd_box(2,3)-cbnd_box(1,3)
#endif /* CYLINDRICAL */
#ifdef POLAR
#endif /* POLAR */

           cell_leng3(i,j,k) = cell_leng3_1 * cell_leng3_2 
     .                                      * cell_leng3_3

         enddo
         enddo
         enddo

!-----------------------
#endif /* CURVILINEAR */

         return
         end subroutine amr_block_geometry
