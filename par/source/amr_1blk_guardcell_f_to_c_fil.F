!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_guardcell_f_to_c_fil(
     .                                  mype,pe,lb,iblock,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit none

!------------------------------------------------------------------------
!
! This routine manages the filling of guard cell information between
! blocks required to fill guard cells on block (pe,lb), assuming that 
! exchange is only required to the current block, from neighbor blocks
! at a higher refinement level. The basis data has been established
! by preceeding calls to guardcell routines, and this routine performs
! the final interpolation required.
! A default, linear interpolation formula is provided as an example.
!
!
! Written :     Peter MacNeice          January 2003
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor number
!      pe             processor address of the selected block
!      lb             local address on proc. pe of the selected block
!      iblock         selects the storage space in data_1blk.fh which is to
!                      be used in this call. If the leaf node is having its
!                      guardcells filled then set this to 1, if its parent
!                      is being filled set it to 2.
!      iopt           a switch to control which data source is to be used
!                      iopt=1 will use 'unk'
!                      iopt=2 will use 'work'
!      nlayers        the number of guard cell layers at each boundary
!      surrblks       the list of addresses of blocks surrounding block lb
!      lcc            a logical switch controlling whether unk or work data
!                      is filled
!      lfc            a logical switch controlling whether facevar data
!                      is filled
!      lec            a logical switch controlling whether unk_e_x(y)(z) data
!                      is filled
!      lnc            a logical switch controlling whether unk_n data
!                      is filled
!      icoord         an integer switch used to select which faces of
!                      the block are to be considered. If icoord=0 all
!                      faces are considered. If icoord=1 only faces perp.
!                      to the y-axis are considered, if icoord=2 only faces
!                      perp. to the x-axis are considered, and if icoord=3
!                      only faces perp. to the z-axis are considered.
!      ldiag          a logical switch which controls whether guardcells
!                      corresponding to neighbor blocks diagonally opposite
!                      block edges and corners are filled.
!
!------------------------------------

#include "amr_shmem.fh"


      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag

!------------------------------------
! local arrays

      real    :: c1(2),c2(2),c3(2),denom(2)
      
      integer :: i,j,k,i0,j0,k0,ii,jj,kk
      integer :: nblk_ind,nblk1,nblk2,n0
      integer :: remote_block,remote_pe
      integer :: il,iu,jl,ju,kl,ku,istr,jstr,kstr
      integer :: select_interpolation 
      integer :: nb0,nb1,nb2,ic,jside,kside

!------------------------------------

        nblk1 = 1
        nblk2 = 27
        do nblk_ind = nblk1,nblk2


         kk = (nblk_ind-1)/9 + 1
         n0 = nblk_ind - (kk-1)*9
         jj = (n0 - 1)/3 + 1
         ii = n0 - (jj - 1)*3
         if(surrblks(3,ii,jj,kk).eq.2) then

           remote_block = surrblks(1,ii,jj,kk)
           remote_pe    = surrblks(2,ii,jj,kk)

!------------------------------------------------
! Apply interpolation operator
!------------------------------------------------
!

        select_interpolation = 2

        if (select_interpolation == 1) then
          call default_restriction
        elseif (select_interpolation == 2) then
          call quadratic_restriction
        elseif (select_interpolation == 3) then
          call gr_restriction
        end if

!------------------------------------------------
        endif                         ! end of surrblks if test
        enddo                         ! end of loop over surrounding neighs


      return


      contains



!-------------------------------------------------------------
!-------------------------------------------------------------
      subroutine default_restriction

      implicit none


! This is the linear restriction which fills guardcells next to
! finer neighbors

! NEEDS MODIFICATION TO SUPPORT CURVILINEAR COORDS


#ifdef CURVILINEAR
! compute geometry variables for the child block (remote_block,remote_pe)
         call amr_block_geometry(remote_block,remote_pe)

         if(iopt.eq.1) then

! Compute volume weighted cell center data for conservative restriction
           if(lcc) then
           do ivar = 1,nvar
             do k = f2c_ind_unk(1,3,nblk_ind),
     .              f2c_ind_unk(2,3,nblk_ind),2
             do j = f2c_ind_unk(1,2,nblk_ind),
     .              f2c_ind_unk(2,2,nblk_ind),2
             do i = f2c_ind_unk(1,1,nblk_ind),
     .              f2c_ind_unk(2,1,nblk_ind),2
               unk1_fl(ivar,i:i+1,j:j+k2d,k:k+k3d) = 
     .                        unk1_fl(ivar,i:i+1,j:j+k2d,k:k+k3d) 
     .                       *cell_vol(i:i+1,j:j+k2d,k:k+k3d)
             enddo
             enddo
             enddo
           enddo
           endif
! Compute area weighted cell face-center data for conservative restriction
           if(lfc) then
           do ivar = 1,nfacevar
             do k = f2c_ind_facex(1,3,nblk_ind),
     .              f2c_ind_facex(2,3,nblk_ind),2
             do j = f2c_ind_facex(1,2,nblk_ind),
     .              f2c_ind_facex(2,2,nblk_ind),2
             do i = f2c_ind_facex(1,1,nblk_ind),
     .              f2c_ind_facex(2,1,nblk_ind)
               if(mod(i,2).eq.1) then
                 facevarx1(ivar,i,j:j+k2d,k:k+k3d,1) = 
     .                      facevarx1(ivar,i,j:j+k2d,k:k+k3d,1) 
     .                     *cell_area1(i,j:j+k2d,k:k+k3d)
               endif
             enddo
             enddo
             enddo

             if(ndim.ge.2) then
               do k = f2c_ind_facey(1,3,nblk_ind),
     .                f2c_ind_facey(2,3,nblk_ind),2
               do j = f2c_ind_facey(1,2,nblk_ind),
     .                f2c_ind_facey(2,2,nblk_ind)
               do i = f2c_ind_facey(1,1,nblk_ind),
     .                f2c_ind_facey(2,1,nblk_ind),2
                 if(mod(j,2).eq.1) then
                   facevary1(ivar,i:i+1,j,k:k+k3d,1) = 
     .                       facevary1(ivar,i:i+1,j,k:k+k3d,1) 
     .                      *cell_area2(i:i+1,j,k:k+k3d)
                 endif
               enddo
               enddo
               enddo
             endif
             if(ndim.eq.3) then
               do k = f2c_ind_facez(1,3,nblk_ind),
     .                f2c_ind_facez(2,3,nblk_ind)
               do j = f2c_ind_facez(1,2,nblk_ind),
     .                f2c_ind_facez(2,2,nblk_ind),2
               do i = f2c_ind_facez(1,1,nblk_ind),
     .                f2c_ind_facez(2,1,nblk_ind),2
                 if(mod(k,2).eq.1) then
                   facevarz1(ivar,i:i+1,j:j+k2d,k,1) = 
     .                        facevarz1(ivar,i:i+1,j:j+k2d,k,1) 
     .                       *cell_area3(i:i+1,j:j+k2d,k)
                 endif
               enddo
               enddo
               enddo
             endif
           enddo
           endif
! Compute distance weighted cell edge-center data for conservative restriction
           if(lec) then
           do ivar = 1,nvaredge
             do k = f2c_ind_unkex(1,3,nblk_ind),
     .              f2c_ind_unkex(2,3,nblk_ind)
             do j = f2c_ind_unkex(1,2,nblk_ind),
     .              f2c_ind_unkex(2,2,nblk_ind)
             do i = f2c_ind_unkex(1,1,nblk_ind),
     .              f2c_ind_unkex(2,1,nblk_ind),2
              if(mod(j,2).eq.1.and.mod(k,2).eq.1) then
                unk_e_x1_fl(ivar,i:i+1,j,k,1) = 
     .                     unk_e_x1_fl(ivar,i:i+1,j,k,1) 
     .                    *cell_leng1(i:i+1,j,k)
              endif
             enddo
             enddo
             enddo

             do k = f2c_ind_unkey(1,3,nblk_ind),
     .              f2c_ind_unkey(2,3,nblk_ind)
             do j = f2c_ind_unkey(1,2,nblk_ind),
     .              f2c_ind_unkey(2,2,nblk_ind),2
             do i = f2c_ind_unkey(1,1,nblk_ind),
     .              f2c_ind_unkey(2,1,nblk_ind)
               if(mod(i,2).eq.1.and.mod(k,2).eq.1) then
                 unk_e_y1_fl(ivar,i,j:j+k2d,k,1) = 
     .                      unk_e_y1_fl(ivar,i,j:j+k2d,k,1) 
     .                     *cell_leng2(i,j:j+k2d,k)
               endif
             enddo
             enddo
             enddo

             if (ndim == 3) then
               do k = f2c_ind_unkez(1,3,nblk_ind),
     .                f2c_ind_unkez(2,3,nblk_ind),2
               do j = f2c_ind_unkez(1,2,nblk_ind),
     .                f2c_ind_unkez(2,2,nblk_ind)
               do i = f2c_ind_unkez(1,1,nblk_ind),
     .                f2c_ind_unkez(2,1,nblk_ind)
                 if(mod(j,2).eq.1.and.mod(i,2).eq.1) then
                   unk_e_z1_fl(ivar,i,j,k:k+k3d,1) =
     .                         unk_e_z1_fl(ivar,i,j,k:k+k3d,1)
     .                        *cell_leng3(i,j,k:k+k3d)
                 endif
               enddo
               enddo
               enddo
             endif

           enddo
           endif

         else

! Compute volume weighted cell center data for conservative restriction
!  of work1.
           ndel = nguard_work - nguard
           do k = f2c_ind_work(1,3,nblk_ind),
     .            f2c_ind_work(2,3,nblk_ind),2
           do j = f2c_ind_work(1,2,nblk_ind),
     .            f2c_ind_work(2,2,nblk_ind),2
           do i = f2c_ind_work(1,1,nblk_ind),
     .            f2c_ind_work(2,1,nblk_ind),2
             work1_fl(i:i+1,j:j+k2d,k:k+k3d,1) =
     .            work1_fl(i:i+1,j:j+k2d,k:k+k3d,1) 
     .           *cell_vol(i-ndel    :i-ndel+1,
     .                     j-ndel*k2d:j-ndel*k2d+k2d,
     .                     k-ndel*k3d:k-ndel*k3d+k3d)
           enddo
           enddo
           enddo

         endif

! now reset geometry factors to appropriate values for the current block lb
         call amr_block_geometry(lb,pe)


#endif /* CURVILINEAR */


! cell centered data
        if(lcc) then

!-------------------------
! Default linear interpolation
!-------------------------
        if(iopt.eq.1) then
!-------------------------

#ifndef GR_TEST

        do k = f2c_ind_unk(1,3,nblk_ind),
     .         f2c_ind_unk(2,3,nblk_ind),2
          k0 = (k+1)/2
        do j = f2c_ind_unk(1,2,nblk_ind),
     .         f2c_ind_unk(2,2,nblk_ind),2
          j0 = (j+1)/2
        do i = f2c_ind_unk(1,1,nblk_ind),
     .         f2c_ind_unk(2,1,nblk_ind),2
          i0 = i/2+1

          if(iblock.eq.1) then
          unk1(1:nvar,i0,j0,k0,iblock) =
     .           ( +unk1_fl(1:nvar,i  ,j    ,k)
     .             +unk1_fl(1:nvar,i+1,j    ,k)
     .             +unk1_fl(1:nvar,i  ,j+k2d,k)
     .             +unk1_fl(1:nvar,i+1,j+k2d,k)
     .             +unk1_fl(1:nvar,i  ,j    ,k+k3d)
     .             +unk1_fl(1:nvar,i+1,j    ,k+k3d)
     .             +unk1_fl(1:nvar,i  ,j+k2d,k+k3d)
     .             +unk1_fl(1:nvar,i+1,j+k2d,k+k3d) )
#ifndef CURVILINEAR
     .           /8.
#else /* CURVILINEAR */
     .           / cell_vol(i0,j0,k0)
#endif /* CURVILINEAR */

           endif

        enddo
        enddo
        enddo

#endif /* GR_TEST */

!-------------------------
        elseif(iopt.gt.1) then
!-------------------------

        do k = f2c_ind_work(1,3,nblk_ind),
     .         f2c_ind_work(2,3,nblk_ind),2
          k0 = (k+1)/2
        do j = f2c_ind_work(1,2,nblk_ind),
     .         f2c_ind_work(2,2,nblk_ind),2
          j0 = (j+1)/2
        do i = f2c_ind_work(1,1,nblk_ind),
     .         f2c_ind_work(2,1,nblk_ind),2
          i0 = i/2+1

          if(iblock.eq.1) then
          work1(i0,j0,k0,iblock) =
     .           ( +work1_fl(i  ,j    ,k)
     .             +work1_fl(i+1,j    ,k)
     .             +work1_fl(i  ,j+k2d,k)
     .             +work1_fl(i+1,j+k2d,k)
     .             +work1_fl(i  ,j    ,k+k3d)
     .             +work1_fl(i+1,j    ,k+k3d)
     .             +work1_fl(i  ,j+k2d,k+k3d)
     .             +work1_fl(i+1,j+k2d,k+k3d) )
#ifndef CURVILINEAR
     .           /8.
#else /* CURVILINEAR */
     .           / cell_vol(i0,j0,k0)
#endif /* CURVILINEAR */
           endif

        enddo
        enddo
        enddo


!-------------------------
        endif                         ! end of iopt if test
!-------------------------

        endif                         ! end of lcc if test

!
! cell-face-centered data
        if(lfc) then

        if(iblock.eq.1) then

          do k = f2c_ind_facex(1,3,nblk_ind),
     .           f2c_ind_facex(2,3,nblk_ind),2
            k0 = (k+1)/2
          do j = f2c_ind_facex(1,2,nblk_ind),
     .           f2c_ind_facex(2,2,nblk_ind),2
            j0 = (j+1)/2
          do i = f2c_ind_facex(1,1,nblk_ind),
     .           f2c_ind_facex(2,1,nblk_ind)
            i0 = i/2+1
            if(mod(i,2).eq.1) then
              facevarx1(1:nbndvar,i0,j0,k0,iblock) = (
     .          +facevarx1_fl(1:nbndvar,i,j,k)
     .          +facevarx1_fl(1:nbndvar,i,j+k2d,k)
     .          +facevarx1_fl(1:nbndvar,i,j,k+k3d)
     .          +facevarx1_fl(1:nbndvar,i,j+k2d,k+k3d) )
#ifndef CURVILINEAR
     .          *.25
#else /* CURVILINEAR */
     .           / cell_area1(i,j,k)
#endif /* CURVILINEAR */

            endif
          enddo
          enddo
          enddo

          if(ndim.ge.2) then
            do k = f2c_ind_facey(1,3,nblk_ind),
     .             f2c_ind_facey(2,3,nblk_ind),2
              k0 = (k+1)/2
            do j = f2c_ind_facey(1,2,nblk_ind),
     .             f2c_ind_facey(2,2,nblk_ind)
              j0 = j/2+1
            do i = f2c_ind_facey(1,1,nblk_ind),
     .             f2c_ind_facey(2,1,nblk_ind),2
              i0 = (i+1)/2
              if(mod(j,2).eq.1) then
                facevary1(1:nbndvar,i0,j0,k0,iblock) = (
     .            +facevary1_fl(1:nbndvar,i,j,k)
     .            +facevary1_fl(1:nbndvar,i+1,j,k)
     .            +facevary1_fl(1:nbndvar,i,j,k+k3d)
     .            +facevary1_fl(1:nbndvar,i+1,j,k+k3d) )
#ifndef CURVILINEAR
     .          *.25
#else /* CURVILINEAR */
     .           / cell_area2(i,j,k)
#endif /* CURVILINEAR */
              endif
            enddo
            enddo
            enddo
          endif

          if(ndim.eq.3) then
            do k = f2c_ind_facez(1,3,nblk_ind),
     .             f2c_ind_facez(2,3,nblk_ind)
              k0 = k/2+1
            do j = f2c_ind_facez(1,2,nblk_ind),
     .             f2c_ind_facez(2,2,nblk_ind),2
              j0 = (j+1)/2
            do i = f2c_ind_facez(1,1,nblk_ind),
     .             f2c_ind_facez(2,1,nblk_ind),2
              i0 = (i+1)/2
              if(mod(k,2).eq.1) then
                facevarz1(1:nbndvar,i0,j0,k0,iblock) = (
     .            +facevarz1_fl(1:nbndvar,i,j,k)
     .            +facevarz1_fl(1:nbndvar,i,j+k2d,k)
     .            +facevarz1_fl(1:nbndvar,i+1,j,k)
     .            +facevarz1_fl(1:nbndvar,i+1,j+k2d,k) )
#ifndef CURVILINEAR
     .          *.25
#else /* CURVILINEAR */
     .           / cell_area3(i,j,k)
#endif /* CURVILINEAR */
              endif
            enddo
            enddo
            enddo
          endif

        endif

        endif                         ! end of lfc if test

!
! cell-edge-centered data
        if(lec) then

          if(iblock.eq.1) then

            do k = f2c_ind_unkex(1,3,nblk_ind),
     .             f2c_ind_unkex(2,3,nblk_ind)
              k0 = 1 + k3d*k/2
            do j = f2c_ind_unkex(1,2,nblk_ind),
     .             f2c_ind_unkex(2,2,nblk_ind)
              j0 = 1 + k2d*j/2
            do i = f2c_ind_unkex(1,1,nblk_ind),
     .             f2c_ind_unkex(2,1,nblk_ind),2
              i0 = (i+1)/2
              if(mod(j,2).eq.1.and.mod(k,2).eq.1) then
                unk_e_x1(1:nbndvare,i0,j0,k0,iblock) = (
     .          + unk_e_x1_fl(1:nbndvare,i,j,k)
     .          + unk_e_x1_fl(1:nbndvare,i+1,j,k) )
#ifndef CURVILINEAR
     .          *.5
#else /* CURVILINEAR */
     .           / cell_leng1(i,j,k)
#endif /* CURVILINEAR */
              endif
            enddo
            enddo
            enddo

            do k = f2c_ind_unkey(1,3,nblk_ind),
     .             f2c_ind_unkey(2,3,nblk_ind)
              k0 = 1 + k3d*k/2
            do j = f2c_ind_unkey(1,2,nblk_ind),
     .             f2c_ind_unkey(2,2,nblk_ind),2
              j0 = (j+1)/2
            do i = f2c_ind_unkey(1,1,nblk_ind),
     .             f2c_ind_unkey(2,1,nblk_ind)
              i0 = 1 + i/2
              if(mod(i,2).eq.1.and.mod(k,2).eq.1) then
                unk_e_y1(1:nbndvare,i0,j0,k0,iblock) = (
     .          + unk_e_y1_fl(1:nbndvare,i,j,k)
     .          + unk_e_y1_fl(1:nbndvare,i,j+k2d,k) )
#ifndef CURVILINEAR
     .          *.5
#else /* CURVILINEAR */
     .           / cell_leng2(i,j,k)
#endif /* CURVILINEAR */
              endif
            enddo
            enddo
            enddo

          if (ndim == 3) then
            do k = f2c_ind_unkez(1,3,nblk_ind),
     .             f2c_ind_unkez(2,3,nblk_ind),2
              k0 = (k+1)/2
            do j = f2c_ind_unkez(1,2,nblk_ind),
     .             f2c_ind_unkez(2,2,nblk_ind)
              j0 = 1 + k2d*j/2
            do i = f2c_ind_unkez(1,1,nblk_ind),
     .             f2c_ind_unkez(2,1,nblk_ind)
              i0 = 1 + i/2
              if(mod(j,2).eq.1.and.mod(i,2).eq.1) then
                unk_e_z1(1:nbndvare,i0,j0,k0,iblock) = (
     .          + unk_e_z1_fl(1:nbndvare,i,j,k)
     .          + unk_e_z1_fl(1:nbndvare,i,j,k+k3d) )
#ifndef CURVILINEAR
     .          *.5
#else /* CURVILINEAR */
     .           / cell_leng3(i,j,k)
#endif /* CURVILINEAR */
              endif
            enddo
            enddo
            enddo
          endif


          endif

        endif                         ! end of lec if test

!
! cell-face-centered data
        if(lnc) then

          if(iblock.eq.1) then
            do k = f2c_ind_unkn(1,3,nblk_ind),
     .             f2c_ind_unkn(2,3,nblk_ind)
              k0 = 1 + k3d*k/2
            do j = f2c_ind_unkn(1,2,nblk_ind),
     .             f2c_ind_unkn(2,2,nblk_ind)
              j0 = 1 + k2d*j/2
            do i = f2c_ind_unkn(1,1,nblk_ind),
     .             f2c_ind_unkn(2,1,nblk_ind)
              i0 = i/2 + 1
              if(mod(i,2).eq.1.and.mod(j,2).eq.1
     .                        .and.mod(k,2).eq.1) then
                unk_n1(1:nvarcorn,i0,j0,k0,iblock)
     .              = unk_n1_fl(1:nvarcorn,i,j,k)
              endif
            enddo
            enddo
            enddo
          endif

        endif                         ! end of lnc if test

      end subroutine default_restriction
!-------------------------------------------------------------
!-------------------------------------------------------------

!-------------------------------------------------------------
!-------------------------------------------------------------
      subroutine quadratic_restriction

      implicit none

      real    :: weight(3,-1:1)
      real    :: xi, xj, www

      integer :: imin,jmin,kmin
      integer :: imax,jmax,kmax
      integer :: iparmin,jparmin,kparmin
      integer :: iparmax,jparmax,kparmax
      integer :: iii,jjj,kkk
      integer :: is,js,ks,iw,jw,kw

! This is the quadratic restriction which fills guardcells next to
! finer neighbors

! NEEDS MODIFICATION TO SUPPORT CURVILINEAR COORDS

      xi = -.5
      do i = -1,1
         weight(1,i) = 1.
         xj = -.5
         do j = -1,1
            if (i .ne. j) then
               weight(1,i) = 
     &              weight(1,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do
      
      xi = -1.5
      do i = -1,1
         weight(2,i) = 1.
         xj = -1.5
         do j = -1,1
            if (i .ne. j) then
               weight(2,i) = 
     &              weight(2,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do


! cell centered data
      if(lcc) then

!-------------------------
! Quadratic interpolation
!-------------------------
        if(iopt.eq.1) then
!-------------------------

        imin = f2c_ind_unk(1,1,nblk_ind)
        imax = f2c_ind_unk(2,1,nblk_ind)
        jmin = f2c_ind_unk(1,2,nblk_ind)
        jmax = f2c_ind_unk(2,2,nblk_ind)
        kmin = f2c_ind_unk(1,3,nblk_ind)
        kmax = f2c_ind_unk(2,3,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2
      
        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if(iblock.eq.1) then

              unk1(1:nvar,i0,j0,k0,iblock) = 0.

              if (k0 <= kparmin + (kparmax-kparmin)/2) then
                 ks = k0*2-1
                 kw = 1
              else
                 ks = k0*2-(1+k3d)
                 kw = 2
              end if
              if (j0 <= jparmin + (jparmax-jparmin)/2) then
                 js = j0*2-1
                 jw = 1
              else
                 js = j0*2-(1+k2d)
                 jw = 2
              end if
              if (i0 <= iparmin + (iparmax-iparmin)/2) then
                 is = i0*2-1
                 iw = 1
              else
                 is = i0*2-2
                 iw = 2
              end if

              k = ks
              do kkk = -k3d,k3d
              j = js
              do jjj = -k2d,k2d
              i = is
              do iii = -1,1

              if (ndim == 1) then
                 www = weight(iw,iii)
              elseif (ndim == 2) then
                 www = weight(iw,iii)*
     .                 weight(jw,jjj)
              elseif (ndim == 3) then
                 www = weight(iw,iii)*
     .                 weight(jw,jjj)*
     .                 weight(kw,kkk)
              end if

              unk1(1:nvar,i0,j0,k0,iblock) = 
     .           unk1(1:nvar,i0,j0,k0,iblock) +
     .           (www*unk1_fl(1:nvar,i,j,k))

              i = i + 1
              end do
              j = j + 1
              end do
              k = k + 1
              end do

           endif

        enddo
        enddo
        enddo

!-------------------------
        elseif(iopt.gt.1) then
!-------------------------

        imin = f2c_ind_work(1,1,nblk_ind)
        imax = f2c_ind_work(2,1,nblk_ind)
        jmin = f2c_ind_work(1,2,nblk_ind)
        jmax = f2c_ind_work(2,2,nblk_ind)
        kmin = f2c_ind_work(1,3,nblk_ind)
        kmax = f2c_ind_work(2,3,nblk_ind)
           
        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2
           
        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if(iblock.eq.1) then

              work1(i0,j0,k0,iblock) = 0.

              if (k0 <= kparmin + (kparmax-kparmin)/2) then
                 ks = k0*2-1
                 kw = 1
              else
                 ks = k0*2-(1+k3d)
                 kw = 2
              end if
              if (j0 <= jparmin + (jparmax-jparmin)/2) then
                 js = j0*2-1
                 jw = 1
              else
                 js = j0*2-(1+k2d)
                 jw = 2
              end if
              if (i0 <= iparmin + (iparmax-iparmin)/2) then
                 is = i0*2-1
                 iw = 1
              else
                 is = i0*2-2
                 iw = 2
              end if

              k = ks
              do kkk = -k3d,k3d
              j = js
              do jjj = -k2d,k2d
              i = is
              do iii = -1,1

              if (ndim == 1) then
                 www = weight(iw,iii)
              elseif (ndim == 2) then
                 www = weight(iw,iii)*
     .                 weight(jw,jjj)
              elseif (ndim == 3) then
                 www = weight(iw,iii)*
     .                 weight(jw,jjj)*
     .                 weight(kw,kkk)
              end if

              work1(i0,j0,k0,iblock) = 
     .           work1(i0,j0,k0,iblock) +
     .           (www*work1_fl(i,j,k))

              i = i + 1
              end do
              j = j + 1
              end do
              k = k + 1
              end do

           endif

        enddo
        enddo
        enddo


!-------------------------
        endif                         ! end of iopt if test
!-------------------------

        endif                         ! end of lcc if test

!
! cell-face-centered data
        if(lfc) then

        if(iblock.eq.1) then

!!! facevarx !!!

        kmin = f2c_ind_facex(1,3,nblk_ind)
        kmax = f2c_ind_facex(2,3,nblk_ind)
        jmin = f2c_ind_facex(1,2,nblk_ind)
        jmax = f2c_ind_facex(2,2,nblk_ind)
        imin = f2c_ind_facex(1,1,nblk_ind)
        imax = f2c_ind_facex(2,1,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2

        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if (k0 <= kparmin + (kparmax-kparmin)/2) then
              ks = k0*2-1
              kw = 1
           else
              ks = k0*2-(1+k3d)
              kw = 2
           end if
           if (j0 <= jparmin + (jparmax-jparmin)/2) then
              js = j0*2-1
              jw = 1
           else
              js = j0*2-(1+k2d)
              jw = 2
           end if
           i = i0*2-1
           
           if (i >= imin .and. i <= imax) then
              
              facevarx1(1:nbndvar,i0,j0,k0,iblock) = 0.
              
              k = ks
              do kkk = -k3d,k3d
              j = js
              do jjj = -k2d,k2d
                    
                 if (ndim == 1) then
                    www = 1.
                 elseif (ndim == 2) then
                    www = 1.*
     .                   weight(jw,jjj)
                 elseif (ndim == 3) then
                    www = 1.*
     .                   weight(jw,jjj)*
     .                   weight(kw,kkk)
                 end if
                 
                 facevarx1(1:nbndvar,i0,j0,k0,iblock) = 
     .                facevarx1(1:nbndvar,i0,j0,k0,iblock) +
     .                (www*facevarx1_fl(1:nbndvar,i,j,k))
                 
              j = j + 1
              end do
              k = k + 1
              end do
              
           end if

        enddo
        enddo
        enddo

!!! facevary !!!

        if(ndim.ge.2) then

        kmin = f2c_ind_facey(1,3,nblk_ind)
        kmax = f2c_ind_facey(2,3,nblk_ind)
        jmin = f2c_ind_facey(1,2,nblk_ind)
        jmax = f2c_ind_facey(2,2,nblk_ind)
        imin = f2c_ind_facey(1,1,nblk_ind)
        imax = f2c_ind_facey(2,1,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2

        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if (k0 <= kparmin + (kparmax-kparmin)/2) then
              ks = k0*2-1
              kw = 1
           else
              ks = k0*2-(1+k3d)
              kw = 2
           end if
           if (i0 <= iparmin + (iparmax-iparmin)/2) then
              is = i0*2-1
              iw = 1
           else
              is = i0*2-2
              iw = 2
           end if
           j = j0*2-1
           
           if (j >= jmin .and. j <= jmax) then
              
              facevary1(1:nbndvar,i0,j0,k0,iblock) = 0.
              
              k = ks
              do kkk = -k3d,k3d
              i = is
              do iii = -1,1

              if (ndim == 1) then
                 www = 1.
              elseif (ndim == 2) then
                 www = 1.*
     .                 weight(iw,iii)
              elseif (ndim == 3) then
                 www = 1.*
     .                 weight(iw,iii)*
     .                 weight(kw,kkk)
              end if

              facevary1(1:nbndvar,i0,j0,k0,iblock) = 
     .           facevary1(1:nbndvar,i0,j0,k0,iblock) +
     .           (www*facevary1_fl(1:nbndvar,i,j,k))

           i = i + 1
           end do
           k = k + 1
           end do

           end if

        enddo
        enddo
        enddo

        endif ! end if (ndim .ge. 2

!!! facevarz !!!

        if(ndim.eq.3) then

        kmin = f2c_ind_facez(1,3,nblk_ind)
        kmax = f2c_ind_facez(2,3,nblk_ind)
        jmin = f2c_ind_facez(1,2,nblk_ind)
        jmax = f2c_ind_facez(2,2,nblk_ind)
        imin = f2c_ind_facez(1,1,nblk_ind)
        imax = f2c_ind_facez(2,1,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2

        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if (j0 <= jparmin + (jparmax-jparmin)/2) then
              js = j0*2-1
              jw = 1
           else
              js = j0*2-(1+k2d)
              jw = 2
           end if
           if (i0 <= iparmin + (iparmax-iparmin)/2) then
              is = i0*2-1
              iw = 1
           else
              is = i0*2-2
              iw = 2
           end if
           k = k0*2-1
           
           if (k >= kmin .and. k <= kmax) then
              
              facevarz1(1:nbndvar,i0,j0,k0,iblock) = 0.
              
              j = js
              do jjj = -k2d,k2d
              i = is
              do iii = -1,1

              www = 1.*
     .              weight(iw,iii)*
     .              weight(jw,jjj)

              facevarz1(1:nbndvar,i0,j0,k0,iblock) = 
     .           facevarz1(1:nbndvar,i0,j0,k0,iblock) +
     .           (www*facevarz1_fl(1:nbndvar,i,j,k))

              i = i + 1
              end do
              j = j + 1
              end do

           end if

        enddo
        enddo
        enddo

        endif                         ! end if (ndim .eq. 3

        endif                         ! end iblock if test

        endif                         ! end of lfc if test

!
! cell-edge-centered data
        if(lec) then

          if(iblock.eq.1) then

!!! ex !!!

        kmin = f2c_ind_unkex(1,3,nblk_ind)
        kmax = f2c_ind_unkex(2,3,nblk_ind)
        jmin = f2c_ind_unkex(1,2,nblk_ind)
        jmax = f2c_ind_unkex(2,2,nblk_ind)
        imin = f2c_ind_unkex(1,1,nblk_ind)
        imax = f2c_ind_unkex(2,1,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2

        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if (i0 <= iparmin + (iparmax-iparmin)/2) then
              is = i0*2-1
              iw = 1
           else
              is = i0*2-2
              iw = 2
           end if
           k = k0*2-1
           j = j0*2-1
           
           if (k >= kmin .and. k <= kmax) then
           if (j >= jmin .and. j <= jmax) then
              
              unk_e_x1(1:nbndvare,i0,j0,k0,iblock) = 0.
              
              i = is
              do iii = -1,1

              www = 1.*
     .              weight(iw,iii)

              unk_e_x1(1:nbndvare,i0,j0,k0,iblock) = 
     .           unk_e_x1(1:nbndvare,i0,j0,k0,iblock) +
     .           (www*unk_e_x1_fl(1:nbndvare,i,j,k))

              i = i + 1
              end do

           end if
           end if

        enddo
        enddo
        enddo

!!! ey !!!

        kmin = f2c_ind_unkey(1,3,nblk_ind)
        kmax = f2c_ind_unkey(2,3,nblk_ind)
        jmin = f2c_ind_unkey(1,2,nblk_ind)
        jmax = f2c_ind_unkey(2,2,nblk_ind)
        imin = f2c_ind_unkey(1,1,nblk_ind)
        imax = f2c_ind_unkey(2,1,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2

        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if (j0 <= jparmin + (jparmax-jparmin)/2) then
              js = j0*2-1
              jw = 1
           else
              js = j0*2-2
              jw = 2
           end if
           k = k0*2-1
           i = i0*2-1
           
           if (k >= kmin .and. k <= kmax) then
           if (i >= imin .and. i <= imax) then
              
              unk_e_y1(1:nbndvare,i0,j0,k0,iblock) = 0.
              
              j = js
              do jjj = -1,1

              www = 1.*
     .              weight(jw,jjj)

              unk_e_y1(1:nbndvare,i0,j0,k0,iblock) = 
     .           unk_e_y1(1:nbndvare,i0,j0,k0,iblock) +
     .           (www*unk_e_y1_fl(1:nbndvare,i,j,k))

              j = j + 1
              end do

           end if
           end if

        enddo
        enddo
        enddo

!!! ez !!!

        if (ndim == 3) then

        kmin = f2c_ind_unkez(1,3,nblk_ind)
        kmax = f2c_ind_unkez(2,3,nblk_ind)
        jmin = f2c_ind_unkez(1,2,nblk_ind)
        jmax = f2c_ind_unkez(2,2,nblk_ind)
        imin = f2c_ind_unkez(1,1,nblk_ind)
        imax = f2c_ind_unkez(2,1,nblk_ind)

        iparmin = (imin+1)/2
        iparmax = (imax+1)/2
        jparmin = (jmin+1)/2
        jparmax = (jmax+1)/2
        kparmin = (kmin+1)/2
        kparmax = (kmax+1)/2

        do k0 = kparmin,kparmax
        do j0 = jparmin,jparmax
        do i0 = iparmin,iparmax

           if (k0 <= kparmin + (kparmax-kparmin)/2) then
              ks = k0*2-1
              kw = 1
           else
              ks = k0*2-2
              kw = 2
           end if
           j = j0*2-1
           i = i0*2-1
           
           if (j >= jmin .and. j <= jmax) then
           if (i >= imin .and. i <= imax) then
              
              unk_e_z1(1:nbndvare,i0,j0,k0,iblock) = 0.
              
              k = ks
              do kkk = -1,1

              www = 1.*
     .              weight(kw,kkk)

              unk_e_z1(1:nbndvare,i0,j0,k0,iblock) = 
     .           unk_e_z1(1:nbndvare,i0,j0,k0,iblock) +
     .           (www*unk_e_z1_fl(1:nbndvare,i,j,k))

              k = k + 1
              end do

           end if
           end if

        enddo
        enddo
        enddo

        endif                       ! end if ndim.eq. 3

        endif                       ! end if iblock

        endif                       ! end of lec if test

!
! cell-face-centered data
        if(lnc) then

          if(iblock.eq.1) then
            do k = f2c_ind_unkn(1,3,nblk_ind),
     .             f2c_ind_unkn(2,3,nblk_ind)
              k0 = 1 + k3d*k/2
            do j = f2c_ind_unkn(1,2,nblk_ind),
     .             f2c_ind_unkn(2,2,nblk_ind)
              j0 = 1 + k2d*j/2
            do i = f2c_ind_unkn(1,1,nblk_ind),
     .             f2c_ind_unkn(2,1,nblk_ind)
              i0 = i/2 + 1
              if(mod(i,2).eq.1.and.mod(j,2).eq.1
     .                        .and.mod(k,2).eq.1) then
                unk_n1(1:nvarcorn,i0,j0,k0,iblock)
     .              = unk_n1_fl(1:nvarcorn,i,j,k)
              endif
            enddo
            enddo
            enddo
          endif

        endif                         ! end of lnc if test

      end subroutine quadratic_restriction
!-------------------------------------------------------------
!-------------------------------------------------------------



!-------------------------------------------------------------
!-------------------------------------------------------------
      subroutine gr_restriction

! This is the 3 point formula, weighted toward the cell interior
! which gives 3rd order accurate interpolation of cell centered
! data.
!
! It is written for uniformly spaced data. 

! For the guardcell next to the block, the stencil is
!   f =  ( 6 * f1 + 10 * f2 - f3 )/15
! in the direction normal to the boundary, and uses
! linear 2 pt averaging in the other directions.
! It uses diagonal directions as the normal at edges and 
! corners.
! For cells further from the block boundary the stencil
! is
!   f =  ( 3 * f1 + 6 * f2 - f3 )/8
! f is the data midway between points 1 and 2.
! Point 1, is the point furthest from the current block.


      implicit none

! cell centered data
        if(lcc) then

!-------------------------
! Default linear interpolation
!-------------------------
        if(iopt.eq.1) then
!-------------------------

        if(iblock.eq.1) then
        
! First inject unk1 values into first interior layers of unk1_fl
! x faces
        do k = 1+2*nguard*k3d,nzb+(nzb+2*nguard)*k3d
          k0 = (k+1)/2
        do j = 1+2*nguard*k2d,nyb+(nyb+2*nguard)*k2d
          j0 = (j+1)/2
        do i = 1+2*nguard,nxb+(nxb+2*nguard), 2*nxb-1
          i0 = (i+1)/2
          unk1_fl(1:nvar,i,j,k) = unk1(1:nvar,i0,j0,k0,iblock)
        enddo
        enddo
        enddo

! Now y faces
        do k = 1+2*nguard*k3d,nzb+(nzb+2*nguard)*k3d
          k0 = (k+1)/2
        do j = 1+2*nguard*k2d,nyb+(nyb+2*nguard)*k2d,nyb+(nyb-1)*k2d
          j0 = (j+1)/2
        do i = 1+2*nguard,nxb+(nxb+2*nguard)
          i0 = (i+1)/2
          unk1_fl(1:nvar,i,j,k) = unk1(1:nvar,i0,j0,k0,iblock)
        enddo
        enddo
        enddo

! Now z faces
        do k = 1+2*nguard*k3d,nzb+(nzb+2*nguard)*k3d,nzb+(nzb-1)*k3d
          k0 = (k+1)/2
        do j = 1+2*nguard*k2d,nyb+(nyb+2*nguard)*k2d
          j0 = (j+1)/2
        do i = 1+2*nguard,nxb+(nxb+2*nguard)
          i0 = (i+1)/2
          unk1_fl(1:nvar,i,j,k) = unk1(1:nvar,i0,j0,k0,iblock)
        enddo
        enddo
        enddo
! Injection of layer 1 is complete

        c1(1) = 6.
        c2(1) = 10.
        c3(1) = -1.
        denom(1) = 15.

        c1(2) = 3.
        c2(2) = 6.
        c3(2) = -1.
        denom(2) = 8.


! Now fill guardcells
        il = 1+nguard
        iu = nxb+nguard
        jl = 1+nguard*k2d
        ju = nyb + nguard*k2d
        kl = 1+nguard*k3d
        ku = nzb + nguard*k3d
        istr = 1
        jstr = 1
        kstr = 1


        if(nblk_ind.eq.13.or.nblk_ind.eq.15) then
! x faces

          il = 1
          iu = nguard
          if(nblk_ind.eq.15) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
          endif

          do k = kl,ku
            k0 = 2*k-1
          do j = jl,ju
            j0 = 2*j-1
          do i = il,iu,istr
            i0 = 2*i-(1+istr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0,j0,k0) 
     .                 + unk1_fl(:,i0,j0+k2d,k0) 
     .                 + unk1_fl(:,i0,j0,k0+k3d) 
     .                 + unk1_fl(:,i0,j0+k2d,k0+k3d) )*.25
     .        +c2(ic)*(+ unk1_fl(:,i0+istr,j0,k0) 
     .                 + unk1_fl(:,i0+istr,j0+k2d,k0) 
     .                 + unk1_fl(:,i0+istr,j0,k0+k3d) 
     .                 + unk1_fl(:,i0+istr,j0+k2d,k0+k3d) )*.25
     .        +c3(ic)*(+ unk1_fl(:,i0+2*istr,j0,k0) 
     .                 + unk1_fl(:,i0+2*istr,j0+k2d,k0) 
     .                 + unk1_fl(:,i0+2*istr,j0,k0+k3d) 
     .                 + unk1_fl(:,i0+2*istr,j0+k2d,k0+k3d) )*.25
     .                )
          enddo
          enddo
          enddo


        elseif(nblk_ind.eq.11.or.nblk_ind.eq.17) then
! y faces

          jl = 1
          ju = nguard
          if(nblk_ind.eq.17) then
            jl = nyb+2*nguard
            ju = nyb+nguard+1
            jstr = -1
          endif


          do k = kl,ku
            k0 = 2*k-1
          do j = jl,ju,jstr
            j0 = 2*j-(1+jstr)/2
          do i = il,iu
            i0 = 2*i-1
            ic = 2
            if(j.eq.ju) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0  ,j0,k0)
     .                 + unk1_fl(:,i0+1,j0,k0)
     .                 + unk1_fl(:,i0  ,j0,k0+k3d)
     .                 + unk1_fl(:,i0+1,j0,k0+k3d) )*.25
     .        +c2(ic)*(+ unk1_fl(:,i0  ,j0+jstr,k0)
     .                 + unk1_fl(:,i0+1,j0+jstr,k0)
     .                 + unk1_fl(:,i0  ,j0+jstr,k0+k3d)
     .                 + unk1_fl(:,i0+1,j0+jstr,k0+k3d) )*.25
     .        +c3(ic)*(+ unk1_fl(:,i0  ,j0+2*jstr,k0)
     .                 + unk1_fl(:,i0+1,j0+2*jstr,k0)
     .                 + unk1_fl(:,i0  ,j0+2*jstr,k0+k3d)
     .                 + unk1_fl(:,i0+1,j0+2*jstr,k0+k3d) )*.25
     .                )
          enddo
          enddo
          enddo

        elseif(nblk_ind.eq.5.or.nblk_ind.eq.23) then
! z faces

          kl = 1
          ku = nguard
          if(nblk_ind.eq.23) then
            kl = nzb+2*nguard
            ku = nzb+nguard+1
            kstr = -1
          endif


          do k = kl,ku,kstr
            k0 = 2*k-(1+kstr)/2
          do j = jl,ju
            j0 = 2*j-1
          do i = il,iu
            i0 = 2*i-1
            ic = 2
            if(k.eq.ku) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .        c1(ic)*(+ unk1_fl(:,i0  ,j0  ,k0)
     .                 + unk1_fl(:,i0+1,j0  ,k0)
     .                 + unk1_fl(:,i0  ,j0+1,k0)
     .                 + unk1_fl(:,i0+1,j0+1,k0) )*.25
     .        +c2(ic)*(+ unk1_fl(:,i0  ,j0  ,k0+kstr)
     .                 + unk1_fl(:,i0+1,j0  ,k0+kstr)
     .                 + unk1_fl(:,i0  ,j0+1,k0+kstr)
     .                 + unk1_fl(:,i0+1,j0+1,k0+kstr) )*.25
     .        +c3(ic)*(+ unk1_fl(:,i0  ,j0  ,k0+2*kstr)
     .                 + unk1_fl(:,i0+1,j0  ,k0+2*kstr)
     .                 + unk1_fl(:,i0  ,j0+1,k0+2*kstr)
     .                 + unk1_fl(:,i0+1,j0+1,k0+2*kstr) )*.25
     .                )
          enddo
          enddo
          enddo

        elseif(nblk_ind.eq.10.or.nblk_ind.eq.18) then
! bottom left and top right z edges
          il = 1
          iu = nguard
          jside = 0
          if(nblk_ind.eq.18) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            jside = 1
          endif
          jstr = istr

          do k = kl,ku
            k0 = 2*k-1
          do i = il,iu,istr
            i0 = 2*i-(1+istr)/2
            j  = i + (nyb - nxb)*jside
            j0 = 2*j-(1+jstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0,j0,k0)
     .                 + unk1_fl(:,i0,j0,k0+k3d) )*.5
     .        +c2(ic)*(+ unk1_fl(:,i0+istr  ,j0+jstr,k0)
     .                 + unk1_fl(:,i0+istr  ,j0+jstr,k0+k3d) )*.5
     .        +c3(ic)*(+ unk1_fl(:,i0+2*istr,j0+2*jstr,k0)
     .                 + unk1_fl(:,i0+2*istr,j0+2*jstr,k0+k3d) )*.5
     .                )
          enddo
          enddo

        elseif(nblk_ind.eq.16.or.nblk_ind.eq.12) then
! top left and bottom right z edges

          il = 1
          iu = nguard
          nb0 = nyb
          if(nblk_ind.eq.12) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            nb0 = nxb
          endif
          jstr = -istr

          do k = kl,ku
            k0 = 2*k-1
          do i = il,iu,istr
            i0 = 2*i-(1+istr)/2
            j  = nb0 + 2*nguard + 1 - i 
            j0 = 2*j-(1+jstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0,j0,k0)
     .                 + unk1_fl(:,i0,j0,k0+k3d) )*.5
     .        +c2(ic)*(+ unk1_fl(:,i0+istr  ,j0+jstr,k0)
     .                 + unk1_fl(:,i0+istr  ,j0+jstr,k0+k3d) )*.5
     .        +c3(ic)*(+ unk1_fl(:,i0+2*istr,j0+2*jstr,k0)
     .                 + unk1_fl(:,i0+2*istr,j0+2*jstr,k0+k3d) )*.5
     .                )
          enddo
          enddo

        elseif(nblk_ind.eq.4.or.nblk_ind.eq.24) then
! bottom left and top right y edges
          il = 1
          iu = nguard
          kside = 0
          if(nblk_ind.eq.24) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            kside = 1
          endif
          kstr = istr

          do j = jl,ju
            j0 = 2*j-1
          do i = il,iu,istr
            i0 = 2*i-(1+istr)/2
            k  = i + (nzb - nxb)*kside
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0,j0    ,k0)
     .                 + unk1_fl(:,i0,j0+k2d,k0) )*.5
     .        +c2(ic)*(+ unk1_fl(:,i0+istr  ,j0    ,k0+kstr)
     .                 + unk1_fl(:,i0+istr  ,j0+k2d,k0+kstr) )*.5
     .        +c3(ic)*(+ unk1_fl(:,i0+2*istr,j0    ,k0+2*kstr)
     .                 + unk1_fl(:,i0+2*istr,j0+k2d,k0+2*kstr) )*.5
     .                )
          enddo
          enddo

        elseif(nblk_ind.eq.22.or.nblk_ind.eq.6) then
! top left and bottom right y edges
          il = 1
          iu = nguard
          nb0 = nzb
          if(nblk_ind.eq.6) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            nb0 = nxb
          endif
          kstr = -istr

          do j = jl,ju
            j0 = 2*j-1
          do i = il,iu,istr
            i0 = 2*i-(1+istr)/2
            k  = nb0 + 2*nguard + 1 - i 
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0,j0    ,k0)
     .                 + unk1_fl(:,i0,j0+k2d,k0) )*.5
     .        +c2(ic)*(+ unk1_fl(:,i0+istr  ,j0    ,k0+kstr)
     .                 + unk1_fl(:,i0+istr  ,j0+k2d,k0+kstr) )*.5
     .        +c3(ic)*(+ unk1_fl(:,i0+2*istr,j0    ,k0+2*kstr)
     .                 + unk1_fl(:,i0+2*istr,j0+k2d,k0+2*kstr) )*.5
     .                )
          enddo
          enddo

        elseif(nblk_ind.eq.2.or.nblk_ind.eq.26) then
! bottom left and top right y edges
          jl = 1
          ju = nguard
          kside = 0
          if(nblk_ind.eq.26) then
            jl = nyb+2*nguard
            ju = nyb+nguard+1
            jstr = -1
            kside = 1
          endif
          kstr = jstr

          do i = il,iu
            i0 = 2*i-1
          do j = jl,ju,jstr
            j0 = 2*j-(1+jstr)/2
            k = j + (nzb - nyb)*kside
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(j.eq.ju) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0  ,j0,k0)
     .                 + unk1_fl(:,i0+1,j0,k0) )*.5
     .        +c2(ic)*(+ unk1_fl(:,i0  ,j0+jstr  ,k0+kstr)
     .                 + unk1_fl(:,i0+1,j0+jstr  ,k0+kstr) )*.5
     .        +c3(ic)*(+ unk1_fl(:,i0  ,j0+2*jstr,k0+2*kstr)
     .                 + unk1_fl(:,i0+1,j0+2*jstr,k0+2*kstr) )*.5
     .                )
          enddo
          enddo

        elseif(nblk_ind.eq.8.or.nblk_ind.eq.20) then
! bottom left and top right y edges
          jl = 1
          ju = nguard
          nb0 = nzb
          if(nblk_ind.eq.8) then
            jl = nyb+2*nguard
            ju = nyb+nguard+1
            jstr = -1
            nb0 = nyb
          endif
          kstr = -jstr

          do i = il,iu
            i0 = 2*i-1
          do j = jl,ju,jstr
            j0 = 2*j-(1+jstr)/2
            k  = nb0 + 2*nguard + 1 - j 
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(j.eq.ju) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)*(+ unk1_fl(:,i0  ,j0,k0)
     .                 + unk1_fl(:,i0+1,j0,k0) )*.5
     .        +c2(ic)*(+ unk1_fl(:,i0  ,j0+jstr  ,k0+kstr)
     .                 + unk1_fl(:,i0+1,j0+jstr  ,k0+kstr) )*.5
     .        +c3(ic)*(+ unk1_fl(:,i0  ,j0+2*jstr,k0+2*kstr)
     .                 + unk1_fl(:,i0+1,j0+2*jstr,k0+2*kstr) )*.5
     .                )
          enddo
          enddo

        elseif(nblk_ind.eq.1.or.nblk_ind.eq.27) then
! bottom front left corner
          il = 1
          iu = nguard
          jside = 0
          kside = 0
          if(nblk_ind.eq.27) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            jside = 1
            kside = 1
          endif
          jstr = istr
          kstr = istr

          do i = il,iu,istr
            j = i + (nyb - nxb)*jside
            k = i + (nzb - nxb)*kside
            i0 = 2*i-(1+istr)/2
            j0 = 2*j-(1+jstr)/2
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)* unk1_fl(:,i0,j0,k0)
     .        +c2(ic)* unk1_fl(:,i0+istr  ,j0+jstr  ,k0+kstr)
     .        +c3(ic)* unk1_fl(:,i0+2*istr,j0+2*jstr,k0+2*kstr)
     .                )
          enddo

        elseif(nblk_ind.eq.7.or.nblk_ind.eq.21) then
! bottom front left corner
          il = 1
          iu = nguard
          nb0 = nyb
          kside = 0
          if(nblk_ind.eq.21) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            nb0 = nxb
            kside = 1
          endif
          jstr = -istr
          kstr = istr


          do i = il,iu,istr
            j = nb0 + 2*nguard + 1 - i 
            k = i + (nzb - nxb)*kside
            i0 = 2*i-(1+istr)/2
            j0 = 2*j-(1+jstr)/2
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)* unk1_fl(:,i0,j0,k0)
     .        +c2(ic)* unk1_fl(:,i0+istr  ,j0+jstr  ,k0+kstr)
     .        +c3(ic)* unk1_fl(:,i0+2*istr,j0+2*jstr,k0+2*kstr)
     .                )
          enddo

        elseif(nblk_ind.eq.19.or.nblk_ind.eq.9) then
! bottom front left corner
          il = 1
          iu = nguard
          nb0 = nzb
          jside = 0
          if(nblk_ind.eq.9) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            nb0 = nxb
            jside = 1
          endif
          jstr = istr
          kstr = -istr

          do i = il,iu,istr
            j = i + (nyb - nxb)*jside
            k = nb0 + 2*nguard + 1 - i
            i0 = 2*i-(1+istr)/2
            j0 = 2*j-(1+jstr)/2
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .         c1(ic)* unk1_fl(:,i0,j0,k0)
     .        +c2(ic)* unk1_fl(:,i0+istr  ,j0+jstr  ,k0+kstr)
     .        +c3(ic)* unk1_fl(:,i0+2*istr,j0+2*jstr,k0+2*kstr)
     .                )
          enddo

        elseif(nblk_ind.eq.25.or.nblk_ind.eq.3) then
! bottom front left corner
          il = 1
          iu = nguard
          nb1 = nyb
          nb2 = nzb
          if(nblk_ind.eq.3) then
            il = nxb+2*nguard
            iu = nxb+nguard+1
            istr = -1
            nb1 = nxb
            nb2 = nxb
          endif
          jstr = -istr
          kstr = -istr

          do i = il,iu,istr
            j = nb1 + 2*nguard + 1 - i
            k = nb2 + 2*nguard + 1 - i
            i0 = 2*i-(1+istr)/2
            j0 = 2*j-(1+jstr)/2
            k0 = 2*k-(1+kstr)/2
            ic = 2
            if(i.eq.iu) ic = 1
            unk1(:,i,j,k,iblock) = (1./denom(ic))*(
     .          c1(ic)* unk1_fl(:,i0,j0,k0)
     .         +c2(ic)* unk1_fl(:,i0+istr  ,j0+jstr  ,k0+kstr)
     .         +c3(ic)* unk1_fl(:,i0+2*istr,j0+2*jstr,k0+2*kstr)
     .                )
          enddo



        endif                        ! end of nblk_ind iftest

        endif                        ! end of iblock iftest


!-------------------------
        elseif(iopt.gt.1) then
!-------------------------

#ifdef NOTYET

! NEEDS UPDATE FOR GR ALGORITHM
! NEEDS UPDATE FOR GR ALGORITHM
! NEEDS UPDATE FOR GR ALGORITHM

        do k = f2c_ind_work(1,3,nblk_ind),
     .         f2c_ind_work(2,3,nblk_ind),2
          k0 = (k+1)/2
        do j = f2c_ind_work(1,2,nblk_ind),
     .         f2c_ind_work(2,2,nblk_ind),2
          j0 = (j+1)/2
        do i = f2c_ind_work(1,1,nblk_ind),
     .         f2c_ind_work(2,1,nblk_ind),2
          i0 = i/2+1

          if(iblock.eq.1) then
          work1(i0,j0,k0,iblock) =
     .           ( +work1_fl(i  ,j    ,k)
     .             +work1_fl(i+1,j    ,k)
     .             +work1_fl(i  ,j+k2d,k)
     .             +work1_fl(i+1,j+k2d,k)
     .             +work1_fl(i  ,j    ,k+k3d)
     .             +work1_fl(i+1,j    ,k+k3d)
     .             +work1_fl(i  ,j+k2d,k+k3d)
     .             +work1_fl(i+1,j+k2d,k+k3d)
     .           /8.
           endif

        enddo
        enddo
        enddo

#endif /* NOTYET */

!-------------------------
        endif                         ! end of iopt if test
!-------------------------

        endif                         ! end of lcc if test


      end subroutine gr_restriction
!-------------------------------------------------------------
!-------------------------------------------------------------

      end subroutine amr_1blk_guardcell_f_to_c_fil
