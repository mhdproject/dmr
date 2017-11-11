!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_restrict_fulltree(mype,iopt,lcc,lfc,lec,lnc)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

#ifndef LIBRARY

      use paramesh_interfaces, only : amr_1blk_copy_soln,
     &                                amr_1blk_guardcell_reset,
     &                                amr_perm_to_1blk,
     &                                amr_1blk_guardcell,
     &                                amr_block_geometry,
     &                                amr_restrict_unk_fun,
     &                                amr_restrict_unk_fun_recip,
     &                                amr_restrict_fc_fun,
     &                                amr_restrict_ec_fun,
     &                                amr_restrict_work_fun,
     &                                amr_restrict_work_fun_recip,
     &                                comm_int_min_to_all,
     &                                comm_int_max_to_all

!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives interior data
! only, not guard cell data. 
! The parent get data for each child and them applies
! the restriction operator to it. Guardcell data may be needed for the
! child blocks, depending on the particular restriction operator being used.
! Thus amr_1blk_guardcell is called below.
! This routine calls a user provided routine called restrict_fun
! which defines the pattern of restriction which the user wishes to
! apply.
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      implicit none

#include "amr_shmem.fh"


      integer, intent(in)  :: mype,iopt
      logical, intent(in)  :: lcc,lfc,lec,lnc

      integer nguard0,nguard_work0
      parameter(nguard0 = nguard*npgs, nguard_work0 = nguard_work*npgs)

!------------------------------------
! local arrays

      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)

      real temp(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
      real send(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

      real recvn0(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                                      kl_bnd:ku_bnd+k3d)
      real recvn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1+k3d)
      real tempn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1+k3d)
      real sendn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1+k3d)

#ifndef LIBRARY
      integer, parameter :: maxbnd = max(nbndvare,nbndvarc,nbndvar)
      real tempf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+k3d)
      real sendf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+k3d)
#else
      integer ::  maxbnd 
      real,allocatable :: tempf(:,:,:,:)
      real,allocatable :: sendf(:,:,:,:)
#endif


      logical l_srl_only,ldiag

      integer remote_pe,remote_block,icoord
      integer cnodetype,cempty
      save    cnodetype,cempty, recvn0

      integer,save :: llrefine_min,llrefine_max
      integer,save :: llrefine_mint,llrefine_maxt

      integer :: lb,level,ich,jchild,ioff,joff,koff,nlayers
      integer :: idest,i,j,k,ii,jj,kk,ivar,iopt0,jface,ng0
      integer :: ia,ja,ka,ib,jb,kb,isa,isb,jsa,jsb,ksa,ksb

!------------------------------------
#ifdef LIBRARY
      maxbnd = max(nbndvare,nbndvarc,nbndvar)
      allocate(
     . tempf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+k3d))
      allocate(
     . sendf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+k3d))
#endif
#ifndef DIAGONALS
      write(*,*) 'amr_1blk_restrict:  diagonals off'
#endif

! For cell-corner data, during a restriction operation, the
! data on a block boundary shared with a neighbor at the same
! refinement level, needs to be acquired during the operation
! of amr_1blk_guardcell_srl for the parent during the call to
! amr_1blk_guardcell_c_to_f. This flag tells amr_1blk_guardcell_srl
! to get this data. Ordinarily amr_1blk_guardcell_srl does not 
! get this data.
      lrestrict_in_progress = .true.

      call amr_1blk_guardcell_reset()

!
! Make sure the gt_unk, gt_facevarx, etc copy of the solution exists
! if NO_PERMANENT_GUARDCELL is defined. This may be needed to fill
! guardcells if the restriction operator needs guardcell data.
      level = -1
      call amr_1blk_copy_soln(level)


! Cycle through parents in decreasing order of refinement
      llrefine_max = maxval(lrefine)
      llrefine_min = llrefine_max
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      llrefine_min = min(lrefine(lb),llrefine_min)
      enddo
      endif
      llrefine_maxt = llrefine_max
      llrefine_mint = llrefine_min
      call shmem_barrier_all()
      call comm_int_max_to_all (llrefine_max,llrefine_maxt)
      call comm_int_min_to_all (llrefine_min,llrefine_mint)
      call shmem_barrier_all()


      if(llrefine_max.gt.llrefine_min) then
      do level = llrefine_max-1,llrefine_min,-1
      call shmem_barrier_all()

! Now parents of leaf nodes get data
! from their children and then perform restriction on it.


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks


! Is this a parent block of at least one leaf node?
      if(nodetype(lb).ge.2.and.lrefine(lb).eq.level) then


! If yes then cycle through its children.
      do ich=1,nchild

        jchild = ich


! Is this child a leaf block? If it is then fetch its data.
        remote_pe     = child(2,ich,lb)
        remote_block  = child(1,ich,lb)
        call shmem_integer_get(cnodetype,nodetype(remote_block),
     .                         1,remote_pe)
        cempty=1
        call shmem_integer_get(cempty,empty(remote_block),
     .                         1,remote_pe)



! compute the offset in the parent block appropriate for this child
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2



! Get the child blocks data and fill its guardcells, putting the result
! into the current working block
           if(iopt.eq.1) then
                nlayers = nguard
           elseif(iopt.ge.2) then
                nlayers = nguard_work
           endif


           if((.not.lnc).and.(.not.lec)) then
            if( lcc.or.lfc.or.(iopt.ge.2)) then
!
! Put child blocks data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding. Note, for even grid sizes the guardcells
! do not need to be filled.
!#ifndef NO_PERMANENT_GUARDCELLS
             idest = 1
             call amr_perm_to_1blk(lcc,lfc,lec,lnc,
     .                             remote_block,remote_pe,
     .                             iopt,idest)

!#endif /* NO_PERMANENT_GUARDCELLS */
            endif
           elseif(lnc.or.lec)  then
!#else /* NO_PERMANENT_GUARDCELLS */
!
! Put child blocks data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding, and fill its guardcells with valid
! data.

             ldiag = .false.
#ifdef DIAGONALS
             ldiag = .true.
#endif
             icoord = 0                          ! fill in all coord directions
!             l_srl_only = .true.
             l_srl_only = .false.
             call amr_1blk_guardcell(mype,iopt,nlayers,
     .                               remote_block,remote_pe,
     .                               lcc,lfc,lec,lnc,
     .                               l_srl_only,icoord,ldiag)
!#endif /* NO_PERMANENT_GUARDCELLS */
           endif


!-----------------------
#ifdef CURVILINEAR
! compute geometry variables for the child block (remote_block,remote_pe)
         call amr_block_geometry(remote_block,remote_pe)

         if(iopt.eq.1) then

! Compute volume weighted cell center data for conservative restriction
           do ivar = 1,nvar
             unk1(ivar,:,:,:,1) = unk1(ivar,:,:,:,1)*cell_vol(:,:,:)
           enddo
! Compute area weighted cell face-center data for conservative restriction
           do ivar = 1,nfacevar
             facevarx1(ivar,:,:,:,1) = facevarx1(ivar,:,:,:,1)
     .                            *cell_area1(:,:,:)
             facevary1(ivar,:,:,:,1) = facevary1(ivar,:,:,:,1)
     .                            *cell_area2(:,:,:)
             facevarz1(ivar,:,:,:,1) = facevarz1(ivar,:,:,:,1)
     .                            *cell_area3(:,:,:)
           enddo
! Compute distance weighted cell edge-center data for conservative restriction
           do ivar = 1,nvaredge
             unk_e_x1(ivar,:,:,:,1) = unk_e_x1(ivar,:,:,:,1)
     .                            *cell_leng1(:,:,:)
             unk_e_y1(ivar,:,:,:,1) = unk_e_y1(ivar,:,:,:,1)
     .                            *cell_leng2(:,:,:)
             unk_e_z1(ivar,:,:,:,1) = unk_e_z1(ivar,:,:,:,1)
     .                            *cell_leng3(:,:,:)
           enddo

         else

! Compute volume weighted cell center data for conservative restriction
!  of work1.
           ndel = nguard_work - nguard
           do k=kl_bnd1+nguard*k3d,ku_bnd1+nguard*k3d
           do j=jl_bnd1+nguard*k2d,ju_bnd1+nguard*k2d
           do i=il_bnd1+nguard    ,iu_bnd1+nguard
             work1(i+ndel,j+ndel*k2d,k+ndel*k3d,1) = 
     .            work1(i+ndel,j+ndel*k2d,k+ndel*k3d,1) 
     .                                        *cell_vol(i,j,k)
           enddo
           enddo
           enddo

         endif

! now reset geometry factors to appropriate values for the current block lb
         call amr_block_geometry(lb,mype)

#endif /* CURVILINEAR */
!-----------------------



         if(iopt.eq.1) then

! Compute restricted cell-centered data from the data in the buffer
           if(lcc) then
             call amr_restrict_unk_fun(unk1(:,:,:,:,1),temp,lb)

           do k=1+nguard*k3d,nzb+nguard*k3d,2
             kk = (k-nguard*k3d)/2+1+nguard*k3d
             do j=1+nguard*k2d,nyb+nguard*k2d,2
               jj = (j-nguard*k2d)/2+1+nguard*k2d
               do i=1+nguard,nxb+nguard,2
                 ii = (i-nguard)/2+1+nguard
                 do ivar=1,nvar
                   send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
                 enddo
               enddo
             enddo
           enddo

! update the parent block
           do k=1,nzb+(-nzb/2)*k3d
             do j=1,nyb+(-nyb/2)*k2d
               do i=1,nxb-nxb/2
                 do ivar=1,nvar
                   unk(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                      k+nguard0*k3d+koff,lb) =
     .                send(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_vol(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
                 enddo
               enddo
             enddo
           enddo

           endif                ! end of lcc iftest


! Compute restricted cell corner data from the data in the buffer
           if(lnc) then

             call amr_restrict_unk_fun_recip( unk_n1(:,:,:,:,1),
     .                                        tempn )

           do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
             kk = (k-nguard*k3d)/2+1+nguard*k3d
             do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
               jj = (j-nguard*k2d)/2+1+nguard*k2d
               do i=1+nguard,nxb+nguard+1,2
                 ii = (i-nguard)/2+1+nguard
                 do ivar=1,nvarcorn
                   sendn(ivar,ii,jj,kk) = tempn(ivar,i,j,k)
                 enddo
               enddo
             enddo
           enddo


! update the parent block
           do k=1,nzb+(-nzb/2+1)*k3d
             do j=1,nyb+(-nyb/2+1)*k2d
               do i=1,nxb-nxb/2+1
                 do ivar=1,nvarcorn
                   unk_n(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                           k+nguard0*k3d+koff,lb)=
     .              sendn(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                 enddo
               enddo
             enddo
           enddo


           endif                 ! end of lnc if test
           endif                ! end of iopt iftest


! Compute restricted cell-face-centered data from the data in the buffer
           if(lfc) then
! Compute restricted data from the data in the buffer
       sendf(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1)
     .   =  facevarx1(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1,1)

       call amr_restrict_fc_fun(sendf,tempf,1)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nfacevar
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2)*k3d
         do j=1,nyb+(-nyb/2)*k2d
           do i=1,nxb-nxb/2+1
             do ivar=1,nfacevar
               facevarx(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                       k+nguard0*k3d+koff,lb)=
     .          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_area1(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
             enddo
           enddo
         enddo
       enddo


! y face next
       if(ndim.ge.2) then

! Compute restricted data from the data in the buffer
       sendf(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                        kl_bnd1:ku_bnd1)
     .   = facevary1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                 kl_bnd1:ku_bnd1,1)
       call amr_restrict_fc_fun(sendf,tempf,2)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nfacevar
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2)*k3d
         do j=1,nyb+(-nyb/2+1)*k2d
           do i=1,nxb-nxb/2
             do ivar=1,nfacevar
               facevary(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                       k+nguard0*k3d+koff,lb)=
     .          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_area2(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
             enddo
           enddo
         enddo
       enddo
       endif                      ! end of ndim>=2 test

! z face last
       if(ndim.eq.3) then

! Compute restricted data from the data in the buffer
       sendf(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                 kl_bnd1:ku_bnd1+k3d)
     .    = facevarz1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1+k3d,1)
       call amr_restrict_fc_fun(sendf,tempf,3)

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nfacevar
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2+1)*k3d
         do j=1,nyb+(-nyb/2)*k2d
           do i=1,nxb-nxb/2
             do ivar=1,nfacevar
               facevarz(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                       k+nguard0*k3d+koff,lb)=
     .          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_area3(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
             enddo
           enddo
         enddo
       enddo

       endif                      ! end of ndim=3 test

           endif                 ! end of lfc if test

           if (ndim > 1) then
! Compute restricted cell-edge-centered data from the data in the buffer
           if(lec) then

! Compute restricted data from the data in the buffer
       sendf(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                         kl_bnd1:ku_bnd1+k3d)
     .   =  unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                  kl_bnd1:ku_bnd1+k3d,1)
       call amr_restrict_ec_fun(sendf,tempf,1)


       sendf = 0.
       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvaredge
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo


! update the parent block
       do k=1,nzb+(-nzb/2+1)*k3d
         do j=1,nyb+(-nyb/2+1)*k2d
           do i=1,nxb-nxb/2
             do ivar=1,nvaredge
               unk_e_x(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                       k+nguard0*k3d+koff,lb)=
     .          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_leng1(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
             enddo
           enddo
         enddo
       enddo


! y edge next
! Compute restricted data from the data in the buffer
       sendf(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                 kl_bnd1:ku_bnd1+k3d)
     .   = unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                 kl_bnd1:ku_bnd1+k3d,1)
       call amr_restrict_ec_fun(sendf,tempf,2)

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvaredge
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2+1)*k3d
         do j=1,nyb+(-nyb/2)*k2d
           do i=1,nxb-nxb/2+1
             do ivar=1,nvaredge
               unk_e_y(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                       k+nguard0*k3d+koff,lb)=
     .          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_leng2(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
             enddo
           enddo
         enddo
       enddo

       if (ndim == 3) Then
! z edge last
! Compute restricted data from the data in the buffer
       sendf(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1)
     .    = unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1,
     .       jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1,1)
       call amr_restrict_ec_fun(sendf,tempf,3)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvaredge
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2)*k3d
         do j=1,nyb+(-nyb/2+1)*k2d
           do i=1,nxb-nxb/2+1
             do ivar=1,nvaredge
               unk_e_z(ivar,i+nguard0+ioff,j+nguard0*k2d+joff,
     .                       k+nguard0*k3d+koff,lb)=
     .          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
#ifdef CURVILINEAR
     .           / cell_leng3(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
             enddo
           enddo
         enddo
       enddo
       end if

       endif                 ! end of lec if test
       end if


         if(iopt.ge.2) then
           iopt0 = iopt-1


! Compute restricted cell-centered workspace data from the data in the buffer


           if(mod(nxb,2).eq.0) then
             call amr_restrict_work_fun(work1(:,:,:,1),tempw1)
           else
             call amr_restrict_work_fun_recip(work1(:,:,:,1),tempw1)
           endif


           do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
             kk = (k-nguard_work*k3d)/2+1+nguard_work*k3d
             do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
               jj = (j-nguard_work*k2d)/2+1+nguard_work*k2d
               do i=1+nguard_work,nxb+nguard_work,2
                 ii = (i-nguard_work)/2+1+nguard_work
                 recvw1(ii,jj,kk,1) = tempw1(i,j,k)
               enddo
             enddo
           enddo


! update the parent block
           do k=1,nzb+(-nzb/2)*k3d
             do j=1,nyb+(-nyb/2)*k2d
               do i=1,nxb-nxb/2
                 work(i+nguard_work0+ioff,j+nguard_work0*k2d+joff,
     .                k+nguard_work0*k3d+koff,lb,iopt0) =  
     .             recvw1(i+nguard_work,j+nguard_work*k2d,
     .                    k+nguard_work*k3d,1)
#ifdef CURVILINEAR
     .           / cell_vol(i+nguard+ioff,j+nguard*k2d+joff,
     .                      k+nguard*k3d+koff)
#endif /* CURVILINEAR */
               enddo
             enddo
           enddo

       endif                      ! end of iopt if test


      enddo                       ! end of loop over children


! If using odd sized grid blocks then parent copies any face bounding
! a leaf block
      if(iopt.eq.1) then
      if(lnc) then


! cycle through parents neighbors
        do jface = 1,nfaces

! get this neighbors nodetype
        remote_pe     = neigh(2,jface,lb)
        remote_block  = neigh(1,jface,lb)
        if(remote_block.gt.0) then

        call shmem_integer_get(cnodetype,nodetype(remote_block),
     .                         1,remote_pe)

        if(cnodetype.eq.1) then

        if(iopt.eq.1) then
             ng0 = nguard0
        elseif(iopt.ge.2) then
             ng0 = nguard_work0
        endif
        ia = 1+ng0
        ib = nxb+ng0+1
        ja = 1+ng0*k2d
        jb = nyb+(ng0+1)*k2d
        ka = 1+ng0*k3d
        kb = nzb+(ng0+1)*k3d
        isa = 1+ng0
        isb = nxb+ng0+1
        jsa = 1+ng0*k2d
        jsb = nyb+(ng0+1)*k2d
        ksa = 1+ng0*k3d
        ksb = nzb+(ng0+1)*k3d

         if(jface.eq.1) then
           ib  = ia
           isa = isb
         elseif(jface.eq.2) then
           ia  = ib
           isb = isa
         elseif(jface.eq.3) then
           jb  = ja
           jsa = jsb
         elseif(jface.eq.4) then
           ja  = jb
           jsb = jsa
         elseif(jface.eq.5) then
           kb  = ka
           ksa = ksb
         elseif(jface.eq.6) then
           ka  = kb
           ksb = ksa
         endif


! Copy neighbor face into this block
         if(iopt.eq.1) then
           call shmem_real_get( recvn0(1,1,1,1),
     .                          unk_n(1,1,1,1,remote_block),
     .                          len_blockn*nvarcorn,remote_pe)
           unk_n(:,ia:ib,ja:jb,ka:kb,lb) = 
     .                      recvn0(:,isa:isb,jsa:jsb,ksa:ksb)

         endif


        endif                     ! end of cnodetype if test

      endif                       ! end of remote_block if test
      enddo                       ! end of jface loop
      endif                       ! end of lnc if test
       endif                      ! end of iopt if test


      endif
      enddo                       ! end of loop over blocks

      endif



! Make sure that the global copy of the newly restricted data is
! up to date.
      call amr_1blk_copy_soln(level)

      enddo                      ! end of loop over refinement levels
      endif

      lrestrict_in_progress = .false.

#ifdef LIBRARY
      deallocate(tempf)
      deallocate(sendf)
#endif

#endif

      call shmem_barrier_all()

      return
      end subroutine amr_restrict_fulltree
