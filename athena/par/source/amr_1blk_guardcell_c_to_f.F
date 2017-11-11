!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_guardcell_c_to_f(mype,lb,pe,iopt,nlayers,
     .                       surrblks,lcc,lfc,lec,lnc,icoord,ldiag,
     .                       nlayers0x,nlayers0y,nlayers0z)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : amr_1blk_cc_prol_gen_work_fun,
     &                                amr_1blk_cc_prol_gen_unk_fun,
     &                                amr_1blk_fc_prol_gen_fun,
     &                                amr_1blk_ec_prol_gen_fun,
     &                                amr_1blk_nc_prol_gen_fun

!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks required to fill guard cells on block lb, assuming that 
! exchange is only required from surrounding blocks at a coarser 
! refinement level.
!
! If you are using an odd number of grid cells, then the interface
! condition implicit in this routine is that the finer block uses
! data prolonged from its parent on the faces which border coarser
! grid blocks.
!
!
! Written :     Peter MacNeice          July 1998
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor number
!      lb             the selected block
!      pe             processor storing the selected block
!      iopt           a switch to control which data source is to be used
!                      iopt=1 will use 'unk'
!                      iopt>=2 will use 'work'
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

      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayers0x,nlayers0y,nlayers0z


!------------------------------------
! local arrays


        common/blk_pattern/psurrblks(3,3,3,3),pcsurrblks(3,3,3,3)
     .                    ,cparent(2)
        integer psurrblks,pcsurrblks,cparent

        integer :: csurr_blks(3,3,1+2*k2d,1+2*k3d)
        integer, parameter :: maxbnd = max(nbndvar,nbndvare,nbndvarc)
        save csurr_blks

        integer :: cwhich_child
        save cwhich_child

        real recv(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
        real recvfx(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .       kl_bnd1:ku_bnd1)
        real recvfy(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1)
        real recvfz(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .       kl_bnd1:ku_bnd1+k3d)
        real recvex(nbndvare,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+k3d)
        real recvey(nbndvare,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .       kl_bnd1:ku_bnd1+k3d)
        real recvez(nbndvare,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1)
        real recvn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+1)
        real recvff(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .       kl_bnd1:ku_bnd1+k3d)
        save recv,recvfx,recvfy,recvfz,recvff
        save recvex,recvey,recvez,recvn

        logical lcoarse_neigh

        integer parent_blk,parent_pe
        integer remote_blk,remote_pe
        save    remote_blk,remote_pe

        integer :: igc_off_x,igc_off_y,igc_off_z,i,j,k
        integer :: ioff,joff,koff,nsurr_blks,nguard0,nlayers0
        integer :: jfl,jfu,ia,ib,ja,jb,ka,kb
        integer :: ionea,ioneb,jonea,joneb,konea,koneb,ii,jj,kk
        integer :: imod,jmod,kmod,ia1,ib1,ja1,jb1,ka1,kb1

        integer,save :: idataout(1)

!------------------------------------

       igc_off_x = gc_off_x
       igc_off_y = gc_off_y
       igc_off_z = gc_off_z

!
! Does this block have any neighbors at lower refinement level?
      lcoarse_neigh = .false.
      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
         if(surrblks(1,i,j,k).lt.1.and.surrblks(1,i,j,k).gt.-20)
     .        lcoarse_neigh = .true.
      enddo
      enddo
      enddo
      if(.not.lcoarse_neigh) return

!------------------------------------

!      pcache_pe = -1
!      pcache_blk = -1

!
! Does current block have a parent?
      if(lnew_parent) then
        call shmem_integer_get(cparent,parent(1,lb),2,pe)
      endif
      
      if(cparent(1).gt.-1) then

!
! Get parent's address
      parent_blk = cparent(1)
      parent_pe  = cparent(2)

      call shmem_integer_get(idataout(1),which_child(lb),1,pe)
      cwhich_child = idataout(1)
!
! compute the offset in the parent block appropriate for this child
       ioff = mod(cwhich_child-1,2)*nxb/2
       joff = mod((cwhich_child-1)/2,2)*nyb/2
       koff = mod((cwhich_child-1)/4,2)*nzb/2

!
! Get list of blocks surrounding parent
      if(.not.lnew_parent) then
          psurrblks(:,:,:,:) = pcsurrblks(:,:,:,:)
      else
          nsurr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
          call shmem_integer_get(csurr_blks,
     .                           surr_blks(1,1,1,1,parent_blk),
     .                           nsurr_blks,parent_pe)
          psurrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
     .                csurr_blks(:,:,1:1+2*k2d,1:1+2*k3d)

          lnew_parent = .false.
          pcsurrblks(:,:,:,:) = psurrblks(:,:,:,:)
      endif

      endif                               ! end of parent if test

!------------------------------------

      nguard0 = nguard
      nlayers0 = nguard

      if(iopt.ge.2) then 
           nguard0 = nguard_work
!           nlayers0 = nlayers
           nlayers0 = max(nlayers0x,nlayers0y,nlayers0z)
      endif


!------------------------------------


!
! Get data from parent's working block
          if(iopt.eq.1.and.lcc) then
              call shmem_real_get(recv(1,1,1,1),unk1(1,1,1,1,2),
     .                            len_block1,mype)
          elseif(iopt.ge.2) then
              call shmem_real_get(recvw1(1,1,1,1),work1(1,1,1,2),
     .                            len_wblock1,mype)
          endif

          if(lfc.and.iopt.eq.1) then

            call shmem_real_get( recvfx(1,1,1,1),facevarx1(1,1,1,1,2),
     .                           len_blockfx1*nbndvar,mype )

            if(ndim.ge.2)
     .      call shmem_real_get( recvfy(1,1,1,1),facevary1(1,1,1,1,2),
     .                           len_blockfy1*nbndvar,mype )

            if(ndim.eq.3) 
     .      call shmem_real_get( recvfz(1,1,1,1),facevarz1(1,1,1,1,2),
     .                           len_blockfz1*nbndvar,mype )

          endif
          if(lec.and.iopt.eq.1) then
!!!#if N_DIM > 1
           if (ndim > 1) then
           call shmem_real_get( recvex(1,1,1,1),unk_e_x1(1,1,1,1,2),
     .                           len_blockex1*nbndvare,mype )
           call shmem_real_get( recvey(1,1,1,1),unk_e_y1(1,1,1,1,2),
     .                           len_blockey1*nbndvare,mype )
!!!#endif /* N_DIM > 1 */
           end if
!!!#if N_DIM == 3
           if (ndim == 3) then
           call shmem_real_get( recvez(1,1,1,1),unk_e_z1(1,1,1,1,2),
     .                           len_blockez1*nbndvare,mype )
!!!#endif /* N_DIM == 3 */
           end if
          endif
          if(lnc.and.iopt.eq.1) then
            call shmem_real_get( recvn(1,1,1,1),unk_n1(1,1,1,1,2),
     .                           len_blockn1*nvarcorn,mype )
          endif


!------------------------------------

!
! First deal with block's regular faces

      jfl = 1
      jfu = nfaces
      if(icoord.gt.0) then
        jfl = 1 + 2*(icoord-1)
        jfu = jfl + 1
      endif

! cycle through block faces
      do jface = jfl,jfu


        ia = 1 + nguard0
        ib = nxb + nguard0
        ja = 1 + nguard0*k2d
        jb = nyb + nguard0*k2d
        ka = 1 + nguard0*k3d
        kb = nzb + nguard0*k3d
        ionea = 0
        ioneb = 1
        jonea = 0
        joneb = k2d
        konea = 0
        koneb = k3d

        if(jface.eq.1) then
          remote_blk = surrblks(1,1,2,2)
          remote_pe  = surrblks(2,1,2,2)
          ia   = 1 + nguard0 - nlayers0
          ib   = nguard0 + igc_off_x
        elseif(jface.eq.2) then
          remote_blk = surrblks(1,3,2,2)
          remote_pe  = surrblks(2,3,2,2)
          ia   = 1 + nxb + nguard0 - igc_off_x
          ib   = nxb + 2*nguard0
        elseif(jface.eq.3) then
          remote_blk = surrblks(1,2,1,2)
          remote_pe  = surrblks(2,2,1,2)
          ja   = 1 + nguard0 - nlayers0
          jb   = nguard0 + igc_off_y
        elseif(jface.eq.4) then
          remote_blk = surrblks(1,2,3,2)
          remote_pe  = surrblks(2,2,3,2)
          ja   = 1 + nyb + nguard0 - igc_off_y
          jb   = nyb + 2*nguard0
        elseif(jface.eq.5) then
          remote_blk = surrblks(1,2,2,1)
          remote_pe  = surrblks(2,2,2,1)
          ka   = 1 + nguard0 - nlayers0
          kb   = nguard0 + igc_off_z
        elseif(jface.eq.6) then
          remote_blk = surrblks(1,2,2,3)
          remote_pe  = surrblks(2,2,2,3)
          ka   = 1 + nzb + nguard0 - igc_off_z
          kb   = nzb + 2*nguard0
        endif

!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.-20.and.remote_blk.lt.0) then


! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb,
     .                                  1,ioff,joff,koff,mype,
     .                                  lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(recvw1(:,:,:,1),
     .                                  ia,ib,ja,jb,ka,kb,
     .                                  1,ioff,joff,koff,mype,
     .                                  lb,parent_pe,parent_blk,
     .                                  interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

            recvff(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                       kl_bnd1:ku_bnd1)
     .      = recvfx(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1)

            if(jface.eq.2) ionea = 1
            if(jface.eq.1) ioneb = 0
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia+ionea,ib+ioneb,ja,jb,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            recvff(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
     .      = recvfy(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
            if(jface.eq.4) jonea = 1
            if(jface.eq.3) joneb = 0
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja+jonea,jb+joneb,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            recvff(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                           kl_bnd1:ku_bnd1+k3d)
     .      = recvfz(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                           kl_bnd1:ku_bnd1+k3d)

            if(jface.eq.6) konea = 1
            if(jface.eq.5) koneb = 0
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb,ka+konea,kb+koneb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,3)
            endif

          endif

!!!#if N_DIM > 1
          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            recvff(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvex(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)


            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb+joneb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,1)


            recvff(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                 kl_bnd1:ku_bnd1+k3d)
     .      = recvey(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,2)


!!!#if N_DIM == 3
            if (ndim == 3) then
            recvff(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
     .     =recvez(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb+joneb,ka,kb,1,
     .                    ioff,joff,koff,mype,3)
!!!#endif /* N_DIM == 3 */
            end if

          endif

!!!#endif /* N_DIM > 1 */
          end if

          if(lnc.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            recvff(:nvarcorn,il_bnd1:iu_bnd1+1,
     .               jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d)
     .      = recvn(:nvarcorn,il_bnd1:iu_bnd1+1,
     .               jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d)

            if(jface.eq.2) ionea = 1
            if(jface.eq.1) ioneb = 0
            if(jface.eq.4) jonea = 1
            if(jface.eq.3) joneb = 0
            if(jface.eq.6) konea = 1
            if(jface.eq.5) koneb = 0
            call amr_1blk_nc_prol_gen_fun( recvff,
     .                    ia+ionea,ib+ioneb,ja+jonea,jb+joneb,
     .                    ka+konea,kb+koneb,1,
     .                    ioff,joff,koff,mype)

          endif

        endif

      enddo                           ! end of loop over block faces

!------------------------------------

      if(ldiag) then

!------------------------------------
      if(ndim.ge.2) then
      if(icoord.ne.3) then
!
! Now fill from edges along the z axis.

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      do jj = 1,3,2
      do ii = 1,3,2
    
          remote_blk = surrblks(1,ii,jj,2)
          remote_pe  = surrblks(2,ii,jj,2)

          imod = ii/2
          jmod = jj/2*k2d

          ia = 1 + (nguard0 - nlayers0) + 
     .             (nxb + nlayers0 - igc_off_x )*imod
          ib = nguard0 + igc_off_x + 
     .             (nxb + nlayers0 - igc_off_x )*imod
          ja = 1 + (nguard0 - nlayers0) + 
     .             (nyb + nlayers0 - igc_off_y )*jmod
          jb = nguard0 + igc_off_y +
     .             (nyb + nlayers0 - igc_off_y )*jmod

          ia1 = imod
          ib1 = imod
          ja1 = jmod
          jb1 = jmod


        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb,
     .                                    1,ioff,joff,koff,mype,
     .                                    lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(recvw1(:,:,:,1),
     .                                    ia,ib,ja,jb,ka,kb,
     .                                    1,ioff,joff,koff,mype,
     .                                    lb,parent_pe,parent_blk,
     .                                    interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

            recvff(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                        kl_bnd1:ku_bnd1)
     .      = recvfx(1:nfacevar,il_bnd1:iu_bnd1+1,
     .                  jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia+ia1,ib+ib1,ja,jb,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            recvff(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .               kl_bnd1:ku_bnd1)
     .      = recvfy(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja+ja1,jb+jb1,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            recvff(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
     .      = recvfz(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb,ka,kb+1,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,3)
            endif

          endif

!!!#if N_DIM > 1
          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            recvff(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvex(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb+joneb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,1)

            recvff(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                 kl_bnd1:ku_bnd1+k3d)
     .      = recvey(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,2)

!!!#if N_DIM == 3
            if (ndim == 3) then
            recvff(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
     .     =recvez(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb+joneb,ka,kb,1,
     .                    ioff,joff,koff,mype,3)
!!!#endif /* N_DIM == 3 */
            end if

          endif

!!!#endif /* N_DIM > 1 */
          end if

          if(lnc.and.iopt.eq.1) then

            recvff(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .     = recvn(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_nc_prol_gen_fun( recvff,
     .                    ia+imod,ib+imod,ja+jmod,jb+jmod,
     .                    ka,kb+k3d,1,
     .                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo


      endif                                   ! end of if test icoord=3
      endif                                   ! end of if test for 2D corners


      if(ndim.eq.3) then
      if(icoord.ne.2) then

!
! Now fill from edges along the y axis.

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      do kk = 1,3,2
      do ii = 1,3,2
    
          remote_blk = surrblks(1,ii,2,kk)
          remote_pe  = surrblks(2,ii,2,kk)

          imod = ii/2
          kmod = kk/2*k3d

          ia = 1 + (nguard0 - nlayers0) + 
     .             (nxb + nlayers0 - igc_off_x )*imod
          ib = nguard0 + igc_off_x + 
     .             (nxb + nlayers0 - igc_off_x )*imod
          ka = 1 + (nguard0 - nlayers0) + 
     .             (nzb + nlayers0 - igc_off_z )*kmod
          kb = nguard0 + igc_off_z +
     .             (nzb + nlayers0 - igc_off_z )*kmod

          ia1 = imod
          ib1 = imod
          ka1 = kmod
          kb1 = kmod


        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb,
     .                                   1,ioff,joff,koff,mype,
     .                                   lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(recvw1(:,:,:,1),
     .                                   ia,ib,ja,jb,ka,kb,
     .                                   1,ioff,joff,koff,mype,
     .                                   lb,parent_pe,parent_blk,
     .                                   interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

            recvff(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                          kl_bnd1:ku_bnd1)
     .      = recvfx(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia+ia1,ib+ib1,ja,jb,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            recvff(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .               kl_bnd1:ku_bnd1)
     .      = recvfy(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb+1,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            recvff(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
     .      = recvfz(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb,ka+ka1,kb+kb1,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,3)
            endif

          endif



!!!#if N_DIM > 1
          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d

            recvff(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvex(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb+joneb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,1)

            recvff(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                 kl_bnd1:ku_bnd1+k3d)
     .      = recvey(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,2)

!!!#if N_DIM == 3
            if (ndim == 3) then
            recvff(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
     .      = recvez(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb+joneb,ka,kb,1,
     .                    ioff,joff,koff,mype,3)
!!!#endif /* N_DIM == 3 */
            end if

          endif

!!!#endif /* N_DIM > 1 */
          end if

          if(lnc.and.iopt.eq.1) then

            recvff(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvn(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_nc_prol_gen_fun( recvff,
     .                    ia+imod,ib+imod,ja,jb+k2d,
     .                    ka+kmod,kb+kmod,1,
     .                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo

      endif                                    ! end of icoord=2 if test
!
! Now fill from edges along the x axis.

      if(icoord.ne.1) then

      ia = 1 + nguard0
      ib = nxb + nguard0
      ja = 1 + nguard0*k2d
      jb = nyb + nguard0*k2d
      ka = 1 + nguard0*k3d
      kb = nzb + nguard0*k3d

! Loop over the 4 corners
      do kk = 1,3,2
      do jj = 1,3,2
    
          remote_blk = surrblks(1,2,jj,kk)
          remote_pe  = surrblks(2,2,jj,kk)

          jmod = jj/2*k2d
          kmod = kk/2*k3d

          ja = 1 + (nguard0 - nlayers0) + 
     .             (nyb + nlayers0 - igc_off_y )*jmod
          jb = nguard0 + igc_off_y +
     .             (nyb + nlayers0 - igc_off_y )*jmod
          ka = 1 + (nguard0 - nlayers0) + 
     .             (nzb + nlayers0 - igc_off_z )*kmod
          kb = nguard0 + igc_off_z +
     .             (nzb + nlayers0 - igc_off_z )*kmod

          ja1 = jmod
          jb1 = jmod
          ka1 = kmod
          kb1 = kmod

        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb,
     .                                    1,ioff,joff,koff,mype,
     .                                    lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(recvw1(:,:,:,1),
     .                                     ia,ib,ja,jb,ka,kb,
     .                                     1,ioff,joff,koff,mype,
     .                                     lb,parent_pe,parent_blk,
     .                                     interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

            recvff(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1)
     .      = recvfx(1:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib+1,ja,jb,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            recvff(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .               kl_bnd1:ku_bnd1)
     .      = recvfy(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja+ja1,jb+jb1,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            recvff(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
     .      = recvfz(1:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb,ka+ka1,kb+kb1,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,3)
            endif

          endif

!!!#if N_DIM > 1
          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d
            recvff(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvex(1:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb+joneb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,1)

            recvff(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                 kl_bnd1:ku_bnd1+k3d)
     .      = recvey(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,2)

!!!#if N_DIM == 3
            if (ndim == 3) then
            recvff(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
     .     =recvez(1:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb+joneb,ka,kb,1,
     .                    ioff,joff,koff,mype,3)
!!!#endif /* N_DIM == 3 */
            end if

          endif

!!!#endif /* N_DIM > 1 */
          end if

          if(lnc.and.iopt.eq.1) then

            recvff(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .     = recvn(1:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_nc_prol_gen_fun( recvff,
     .                    ia,ib+1,ja+jmod,jb+jmod,
     .                    ka+kmod,kb+kmod,1,
     .                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo

      endif                                 ! end of if test for icoord=1

      endif                                 ! end of if test for x and y edges

!------------------------------------

! Finally fill corners in 3D.

!------------------------------------
      if(ndim.eq.3) then


! Loop over the 4 corners
      do kk = 1,3,2
      do jj = 1,3,2
      do ii = 1,3,2
    
          remote_blk = surrblks(1,ii,jj,kk)
          remote_pe  = surrblks(2,ii,jj,kk)


          imod = ii/2
          jmod = jj/2
          kmod = kk/2

          ia = 1 + (nguard0 - nlayers0) + 
     .             (nxb + nlayers0 - igc_off_x )*imod
          ib = nguard0 + igc_off_x + 
     .             (nxb + nlayers0 - igc_off_x )*imod
          ja = 1 + (nguard0 - nlayers0) + 
     .             (nyb + nlayers0 - igc_off_y )*jmod
          jb = nguard0 + igc_off_y +
     .             (nyb + nlayers0 - igc_off_y )*jmod
          ka = 1 + (nguard0 - nlayers0) + 
     .             (nzb + nlayers0 - igc_off_z )*kmod
          kb = nguard0 + igc_off_z +
     .             (nzb + nlayers0 - igc_off_z )*kmod

          ia1 = imod
          ib1 = imod
          ja1 = jmod
          jb1 = jmod
          ka1 = kmod
          kb1 = kmod

        if(remote_blk.gt.-20.and.remote_blk.lt.0) then

! interpolate(prolongate) data from the parent to the child
          if(iopt.eq.1.and.lcc) then
             call amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb,
     .                                    1,ioff,joff,koff,mype,
     .                                    lb,parent_pe,parent_blk)
          elseif(iopt.ge.2) then
             call amr_1blk_cc_prol_gen_work_fun(recvw1(:,:,:,1),
     .                                    ia,ib,ja,jb,ka,kb,
     .                                    1,ioff,joff,koff,mype,
     .                                    lb,parent_pe,parent_blk,
     .                                    interp_mask_work(iopt-1))
          endif

          if(lfc.and.iopt.eq.1) then

            recvff(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                         kl_bnd1:ku_bnd1)
     .      = recvfx(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia+ia1,ib+ib1,ja,jb,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,1)

            if(ndim.ge.2) then
            recvff(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .               kl_bnd1:ku_bnd1)
     .      = recvfy(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja+ja1,jb+jb1,ka,kb,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,2)
            endif

            if(ndim.eq.3) then
            recvff(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
     .      = recvfz(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_fc_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb,ka+ka1,kb+kb1,1,
     .                    ioff,joff,koff,
     .                    mype,lb,parent_pe,parent_blk,3)
            endif

          endif


!!!#if N_DIM > 1
          if (ndim > 1) then
          if(lec.and.iopt.eq.1) then

            ionea = 0
            ioneb = 1
            jonea = 0
            joneb = k2d
            konea = 0
            koneb = k3d
            recvff(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvex(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib,ja,jb+joneb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,1)

            recvff(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                 kl_bnd1:ku_bnd1+k3d)
     .      = recvey(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                   kl_bnd1:ku_bnd1+k3d)
            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb,ka,kb+koneb,1,
     .                    ioff,joff,koff,mype,2)

!!!#if N_DIM == 3
            if (ndim == 3) then
            recvff(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)
     .      = recvez(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                           kl_bnd1:ku_bnd1)

            call amr_1blk_ec_prol_gen_fun( recvff,
     .                    ia,ib+ioneb,ja,jb+joneb,ka,kb,1,
     .                    ioff,joff,koff,mype,3)
!!!#endif /* N_DIM == 3 */
            end if

          endif

!!!#endif /* N_DIM > 1 */
          end if

          if(lnc.and.iopt.eq.1) then

            recvff(:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                               kl_bnd1:ku_bnd1+k3d)
     .      = recvn(:nvarcorn,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                   kl_bnd1:ku_bnd1+k3d)

            call amr_1blk_nc_prol_gen_fun( recvff,
     .                    ia+imod,ib+imod,ja+jmod,jb+jmod,
     .                    ka+kmod,kb+kmod,1,
     .                    ioff,joff,koff,mype)

          endif

        endif


      enddo
      enddo
      enddo

      endif                                 ! end of loop over 3D corners

!------------------------------------

      endif                                 ! end of ldiag if test

!------------------------------------

      return
      end subroutine amr_1blk_guardcell_c_to_f
