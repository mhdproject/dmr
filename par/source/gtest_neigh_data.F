!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#define DRL
      subroutine test_neigh_data(mype,istep)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : comm_int_sum_to_all


!------------------------------------------------------------------------
!
! This routine compares face centered data on neighboring blocks at
! a selected point in space. The chosen point (x0,y0,z0) must be within
! one cell width of a block boundary. If there is no refinement jump
! across the block boundary the field value on the chosen cell face
! for both blocks are printed. If there is a refinement jump then the
! 4 face values on the finer block are summed, and the sum is printed
! along with the corresponding face value on the coarse neighbor.
!
! Arguments:
!      mype             local processor
!
!
! Written :     Peter MacNeice          June 1998
!------------------------------------------------------------------------

#include "amr_shmem.fh"


#ifdef TIMINGS
#include "timer.fh"
#endif

!-------------------------

      integer remote_pe,remote_block

! local arrays
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd+k3d)
      save recvx,recvy,recvz


      integer	itface(maxblocks)
      integer   cnodetype
      save cnodetype

      integer shmem_n_pes
      integer par_neigh(2,mfaces)
      integer par_child(2,mchild)
      integer parent_blk,parent_pe
      save par_neigh,par_child

      character*6 pattern,pattern1
      integer     pattern_buf(4)
      integer     lb_buf(4),iprint,jprint,kprint,level_buf(4)
      real        data_buf(4,4)
      save        pattern_buf,lb_buf,data_buf,level_buf
      save        iprint,jprint,kprint

      real :: divbmax,divbnlmax,dx,dy,dz,deli,eps,divbl,divbnl
      real :: divbm,bxsum,bxsum1,bysum,bysum1,bzsum,bzsum1

      integer :: nprocs,lb,iface,neigh_blk,neigh_pe
      integer :: i0,i1,j0,j1,k0,k1,i2,j2,k2,ii,jj,kk
      integer :: ii1,jj1,kk1,ii0,jj0,kk0,ioff,joff,koff
      integer :: jchild,ich,iproc,ip

!-------------------------

! Test point coords.

! srl y
      x0 = .374
      y0 = .25
      z0 = .124
      icoord0 = 2

! srl z
      x0 = .124
      y0 = .376
      z0 = .25
      icoord0 = 3

! srl x 
      x0 = .25
      y0 = .151
      z0 = .067+.25
      icoord0 = 1

! srl y 
      x0 = .124
      y0 = .124
      z0 = .124
      icoord0 = 3

#ifdef DRL
! drl x+
      x0 = .25
      y0 = .151
      z0 = .067
      icoord0 = 1

! drl x-
      x0 = .25
      y0 = .25+.06
      z0 = .25+.06
      icoord0 = 1

! drl y-
      y0 = .25
      x0 = .063+.25
      z0 = .063+.25
      icoord0 = 2

! drl y+
      y0 = .25
      x0 = .06
      z0 = .063
      icoord0 = 2

! drl z-
      z0 = .25
      x0 = .06+.25
      y0 = .063+.25
      icoord0 = 3

! drl z+
      z0 = .25
      x0 = .06
      y0 = .063
      icoord0 = 3

! drl x+
      x0 = .25
      y0 = .06
      z0 = .25
      icoord0 = 1
#endif

      nprocs = shmem_n_pes()

      iprint = 0
      jprint = 0
      kprint = 0
      lb_buf(:) = 0
      level_buf(:) = 0
      pattern_buf(:) = 0
      data_buf(:,:) = 0.

      itface(:) = 0

! Loop over blocks.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then


! Is this the block we wish to test?
        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        dz = bsize(3,lb)/real(nzb)
        deli = 1./dx
        eps = .01*dx
        eps = .5*dx
        xlo = bnd_box(1,1,lb)-eps
        xhi = bnd_box(2,1,lb)+eps
        xtest = (xhi-x0)*(xlo-x0)
        ylo = bnd_box(1,2,lb)-eps
        yhi = bnd_box(2,2,lb)+eps
        ytest = (yhi-y0)*(ylo-y0)
        zlo = bnd_box(1,3,lb)-eps
        zhi = bnd_box(2,3,lb)+eps
        ztest = (zhi-z0)*(zlo-z0)
! Is selected point in a cell next to a block boundary?
        if(xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) then
          do if = 1,nfaces
           xtest = 1.e10
           ytest = 1.e10
           ztest = 1.e10
           if(if.eq.1.and.icoord0.eq.1) then
             xtest = abs(bnd_box(1,1,lb)-x0)
             if(xtest.lt.dx) itface(lb) = if
           elseif(if.eq.2.and.icoord0.eq.1) then
             xtest = abs(bnd_box(2,1,lb)-x0)
             if(xtest.lt.dx) itface(lb) = if
           elseif(if.eq.3.and.icoord0.eq.2) then
             ytest = abs(bnd_box(1,2,lb)-y0)
             if(ytest.lt.dy) itface(lb) = if
           elseif(if.eq.4.and.icoord0.eq.2) then
             ytest = abs(bnd_box(2,2,lb)-y0)
             if(ytest.lt.dy) itface(lb) = if
           elseif(if.eq.5.and.icoord0.eq.3) then
             ztest = abs(bnd_box(1,3,lb)-z0)
             if(ztest.lt.dz) itface(lb) = if
           elseif(if.eq.6.and.icoord0.eq.3) then
             ztest = abs(bnd_box(2,3,lb)-z0)
             if(ztest.lt.dz) itface(lb) = if
           endif
          enddo
        endif


! Locate the required block face.
      if(itface(lb).gt.0) then

! Compute data index in local block.

         i0 = int((x0+eps-bnd_box(1,1,lb))/dx)+1+nguard
         j0 = int((y0+eps-bnd_box(1,2,lb))/dy)+1+nguard
         k0 = int((z0+eps-bnd_box(1,3,lb))/dz)+1+nguard
          

         neigh_blk = neigh(1,itface(lb),lb)
         neigh_pe  = neigh(2,itface(lb),lb)


! Is the neighbor at the same refinement level?
      if(neigh_blk.gt.0) then

          call shmem_integer_get(cnodetype,nodetype(neigh_blk),
     .       1,neigh_pe)


! Only consider neighbor if it is a leaf node.
          if(cnodetype.eq.1) then


! Compute data index for corresponding cell on neighbor.
         i1 = i0
         j1 = j0
         k1 = k0
         if(itface(lb).eq.1) then
           i1 = nxb+nguard+1
         elseif(itface(lb).eq.2) then
           i1 = nguard+1
         elseif(itface(lb).eq.3) then
           j1 = nyb+nguard+1
         elseif(itface(lb).eq.4) then
           j1 = nguard+1
         elseif(itface(lb).eq.5) then
           k1 = nzb+nguard+1
         elseif(itface(lb).eq.6) then
           k1 = nguard+1
         endif


! Copy complete remote block into a buffer block called recv1.
! Only consider faces 2, 4 and 6, to avoid redundant check.
         if(itface(lb).eq.2) then

           call shmem_real_get(recvx(1,1,1,1),
     .                          facevarx(1,1,1,1,neigh_blk),
     .                          len_blockfx*nbndvar,neigh_pe) 
           divbl = deli*(
     .       facevarx(1,i0,j0,k0,lb)    -facevarx(1,i0-1,j0,k0,lb)
     .     + facevary(1,i0-1,j0+1,k0,lb)-facevary(1,i0-1,j0,k0,lb)
     .     + facevarz(1,i0-1,j0,k0+1,lb)-facevarz(1,i0-1,j0,k0,lb))

           divbnl = deli*(
     .       recvx(1,i1,j1,k1)          -facevarx(1,i0-1,j0,k0,lb)
     .     + facevary(1,i0-1,j0+1,k0,lb)-facevary(1,i0-1,j0,k0,lb)
     .     + facevarz(1,i0-1,j0,k0+1,lb)-facevarz(1,i0-1,j0,k0,lb))

           iprint = iprint+1
           pattern_buf(iprint) = 1
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = facevarx(1,i0,j0,k0,lb)
           data_buf(2,iprint) = recvx(1,i1,j1,k1)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.4) then

           call shmem_real_get(recvy(1,1,1,1),
     .                          facevary(1,1,1,1,neigh_blk),
     .                          len_blockfy*nbndvar,neigh_pe) 
           divbl = deli*(
     .       facevarx(1,i0+1,j0-1,k0,lb)-facevarx(1,i0,j0-1,k0,lb)
     .     + facevary(1,i0,j0,k0,lb)    -facevary(1,i0,j0-1,k0,lb)
     .     + facevarz(1,i0,j0-1,k0+1,lb)-facevarz(1,i0,j0-1,k0,lb))

           divbnl = deli*(
     .       facevarx(1,i0+1,j0-1,k0,lb)-facevarx(1,i0,j0-1,k0,lb)
     .     + recvy(1,i1,j1,k1)          -facevary(1,i0,j0-1,k0,lb)
     .     + facevarz(1,i0,j0-1,k0+1,lb)-facevarz(1,i0,j0-1,k0,lb))

           iprint = iprint+1
           pattern_buf(iprint) = 2
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = facevary(1,i0,j0,k0,lb)
           data_buf(2,iprint) = recvy(1,i1,j1,k1)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.6) then

           call shmem_real_get(recvz(1,1,1,1),
     .                          facevarz(1,1,1,1,neigh_blk),
     .                          len_blockfz*nbndvar,neigh_pe) 
           divbl = deli*(
     .       facevarx(1,i0+1,j0,k0-1,lb)-facevarx(1,i0,j0,k0-1,lb)
     .     + facevary(1,i0,j0+1,k0-1,lb)-facevary(1,i0,j0,k0-1,lb)
     .     + facevarz(1,i0,j0,k0,lb)    -facevarz(1,i0,j0,k0-1,lb))

           divbnl = deli*(
     .       facevarx(1,i0+1,j0,k0-1,lb)-facevarx(1,i0,j0,k0-1,lb)
     .     + facevary(1,i0,j0+1,k0-1,lb)-facevary(1,i0,j0,k0-1,lb)
     .     + recvz(1,i1,j1,k1)          -facevarz(1,i0,j0,k0-1,lb))

           iprint = iprint+1
           pattern_buf(iprint) = 3
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = facevarz(1,i0,j0,k0,lb)
           data_buf(2,iprint) = recvz(1,i1,j1,k1)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl


         endif


         endif                          ! end of neigh leaf node test

      else

! Now deal with neighbor at coarser level

! Get parents neighbor and child info
        parent_blk = parent(1,lb)
        parent_pe  = parent(2,lb)
        call shmem_real_get(par_neigh,neigh(1,1,parent_blk),
     .                                 2*mfaces,parent_pe)
        call shmem_real_get(par_child,child(1,1,parent_blk),
     .                                 2*mchild,parent_pe)

! Select the correct neighbor
        neigh_blk = par_neigh(1,itface(lb))
        neigh_pe  = par_neigh(2,itface(lb))


! Compute data (index-nguard) range in local block for 4 cell faces.
         i0 = int((x0+eps-bnd_box(1,1,lb))/dx)+1
         j0 = int((y0+eps-bnd_box(1,2,lb))/dy)+1
         k0 = int((z0+eps-bnd_box(1,3,lb))/dz)+1
         i0 = min(i0,nxb)
         j0 = min(j0,nyb)
         k0 = min(k0,nzb)
         i1 = i0-1
         if(mod(i0,2).eq.1) i1 = i0+1
         j1 = j0-1
         if(mod(j0,2).eq.1) j1 = j0+1
         k1 = k0-1
         if(mod(k0,2).eq.1) k1 = k0+1

         i2 = (i0-1)/2+1 + nguard
         j2 = (j0-1)/2+1 + nguard
         k2 = (k0-1)/2+1 + nguard

! Add nguard
         i0 = i0+nguard
         i1 = i1+nguard
         j0 = j0+nguard
         j1 = j1+nguard
         k0 = k0+nguard
         k1 = k1+nguard

! Set block boundary index
         if(itface(lb).eq.1) then
           i0 = nguard+1
         elseif(itface(lb).eq.2) then
           i0 = nxb+nguard+1
         elseif(itface(lb).eq.3) then
           j0 = nguard+1
         elseif(itface(lb).eq.4) then
           j0 = nyb+nguard+1
         elseif(itface(lb).eq.5) then
           k0 = nguard+1
         elseif(itface(lb).eq.6) then
           k0 = nzb+nguard+1
         endif


! identify which child the local block represents
         jchild = 0
         do ich=1,nchild
           if( (par_child(1,ich).eq.lb)
     .       .and.(par_child(2,ich).eq.mype) )
     .       jchild=ich
         enddo

! compute the offset in the parent block appropriate for this child
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2

! flip the offset in the appropriate axis
         if(itface(lb).le.2) then
           ioff = mod(jchild,2)*nxb/2
         elseif(itface(lb).eq.3.or.itface(lb).eq.4) then
           joff = mod((jchild)/2,2)*nyb/2
         elseif(itface(lb).gt.4) then
           koff = mod((jchild)/4,2)*nzb/2
         endif

         i2 = i2+ioff
         j2 = j2+joff
         k2 = k2+koff


! Copy complete remote block into a buffer block called recv1.
         if(itface(lb).eq.1) then

            call shmem_real_get(recvx(1,1,1,1),
     .                          facevarx(1,1,1,1,neigh_blk),
     .                          len_blockfx*nbndvar,neigh_pe)
            i2 = nxb+nguard+1
            bxsum = facevarx(1,i0,j0,k0,lb)+facevarx(1,i0,j1,k0,lb)
     .             +facevarx(1,i0,j0,k1,lb)+facevarx(1,i0,j1,k1,lb)

            bxsum1=facevarx(1,i0+2,j0,k0,lb)+facevarx(1,i0+2,j1,k0,lb)
     .            +facevarx(1,i0+2,j0,k1,lb)+facevarx(1,i0+2,j1,k1,lb)
            ii = i0
            ii1 = ii+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            kk0 = min(k0,k1)
            kk1 = kk0+1
            bysum =facevary(1,ii,jj0,kk0,lb)+facevary(1,ii1,jj0,kk0,lb)
     .            +facevary(1,ii,jj0,kk1,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1= 
     .         facevary(1,ii,jj0+2,kk0,lb)+facevary(1,ii1,jj0+2,kk0,lb)
     .        +facevary(1,ii,jj0+2,kk1,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            bzsum =facevarz(1,ii,jj0,kk0,lb)+facevarz(1,ii1,jj0,kk0,lb)
     .            +facevarz(1,ii,jj1,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1=
     .         facevarz(1,ii,jj0,kk0+2,lb)+facevarz(1,ii1,jj0,kk0+2,lb)
     .        +facevarz(1,ii,jj1,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)

            divbl  = (bxsum1 - bxsum + bysum1 - bysum + bzsum1 - bzsum)
     .               *deli*.125
            divbnl = (bxsum1 - 4.*recvx(1,i2,j2,k2) + bysum1 - bysum + 
     .               bzsum1 - bzsum)*deli*.125

           iprint = iprint+1
           pattern_buf(iprint) = 4
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bxsum
           data_buf(2,iprint) = recvx(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.2) then

            call shmem_real_get(recvx(1,1,1,1),
     .                          facevarx(1,1,1,1,neigh_blk),
     .                          len_blockfx*nbndvar,neigh_pe)
            i2 = nguard+1
            bxsum = facevarx(1,i0,j0,k0,lb)+facevarx(1,i0,j1,k0,lb)
     .             +facevarx(1,i0,j0,k1,lb)+facevarx(1,i0,j1,k1,lb)

            bxsum1=facevarx(1,i0-2,j0,k0,lb)+facevarx(1,i0-2,j1,k0,lb)
     .            +facevarx(1,i0-2,j0,k1,lb)+facevarx(1,i0-2,j1,k1,lb)
            ii = i0-2
            ii1 = ii+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            kk0 = min(k0,k1)
            kk1 = kk0+1
            bysum =facevary(1,ii,jj0,kk0,lb)+facevary(1,ii1,jj0,kk0,lb)
     .            +facevary(1,ii,jj0,kk1,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1=
     .         facevary(1,ii,jj0+2,kk0,lb)+facevary(1,ii1,jj0+2,kk0,lb)
     .        +facevary(1,ii,jj0+2,kk1,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            bzsum =facevarz(1,ii,jj0,kk0,lb)+facevarz(1,ii1,jj0,kk0,lb)
     .            +facevarz(1,ii,jj1,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1=
     .         facevarz(1,ii,jj0,kk0+2,lb)+facevarz(1,ii1,jj0,kk0+2,lb)
     .        +facevarz(1,ii,jj1,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)

            divbl  = (bxsum - bxsum1 + bysum1 - bysum + bzsum1 - bzsum)
     .               *deli*.125
            divbnl = (4.*recvx(1,i2,j2,k2) - bxsum1 + bysum1 - bysum + 
     .               bzsum1 - bzsum)*deli*.125

           iprint = iprint+1
           pattern_buf(iprint) = 4
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bxsum
           data_buf(2,iprint) = recvx(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.3) then
            call shmem_real_get(recvy(1,1,1,1),
     .                          facevary(1,1,1,1,neigh_blk),
     .                          len_blockfy*nbndvar,neigh_pe)
            j2 = nyb+nguard+1
            bysum = facevary(1,i0,j0,k0,lb)+facevary(1,i1,j0,k0,lb)
     .        +facevary(1,i0,j0,k1,lb)+facevary(1,i1,j0,k1,lb)

            bysum1= facevary(1,i0,j0+2,k0,lb)+facevary(1,i1,j0+2,k0,lb)
     .        +facevary(1,i0,j0+2,k1,lb)+facevary(1,i1,j0+2,k1,lb)
            jj = j0
            jj1 = jj+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            kk0 = min(k0,k1)
            kk1 = kk0+1
            bxsum =facevarx(1,ii0,jj,kk0,lb)+facevarx(1,ii0,jj1,kk0,lb)
     .            +facevarx(1,ii0,jj,kk1,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1=
     .         facevarx(1,ii0+2,jj,kk0,lb)+facevarx(1,ii0+2,jj1,kk0,lb)
     .        +facevarx(1,ii0+2,jj,kk1,lb)+facevarx(1,ii0+2,jj1,kk1,lb)
            bzsum =facevarz(1,ii0,jj,kk0,lb)+facevarz(1,ii0,jj1,kk0,lb)
     .            +facevarz(1,ii1,jj,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1=
     .         facevarz(1,ii0,jj,kk0+2,lb)+facevarz(1,ii0,jj1,kk0+2,lb)
     .        +facevarz(1,ii1,jj,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)


            divbl  = (bxsum1 - bxsum + bysum1 - bysum + bzsum1 - bzsum)
     .               *deli*.125
            divbnl = (bxsum1 - bxsum + bysum1 - 4.*recvy(1,i2,j2,k2) + 
     .               bzsum1 - bzsum)*deli*.125


           iprint = iprint+1
           pattern_buf(iprint) = 5
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bysum
           data_buf(2,iprint) = recvy(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.4) then
            call shmem_real_get(recvy(1,1,1,1),
     .                          facevary(1,1,1,1,neigh_blk),
     .                          len_blockfy*nbndvar,neigh_pe)
            j2 = nguard+1
            bysum = facevary(1,i0,j0,k0,lb)+facevary(1,i1,j0,k0,lb)
     .             +facevary(1,i0,j0,k1,lb)+facevary(1,i1,j0,k1,lb)

            bysum1= facevary(1,i0,j0-2,k0,lb)+facevary(1,i1,j0-2,k0,lb)
     .             +facevary(1,i0,j0-2,k1,lb)+facevary(1,i1,j0-2,k1,lb)
            jj = j0-2
            jj1 = jj+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            kk0 = min(k0,k1)
            kk1 = kk0+1
            bxsum =facevarx(1,ii0,jj,kk0,lb)+facevarx(1,ii0,jj1,kk0,lb)
     .            +facevarx(1,ii0,jj,kk1,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1=
     .         facevarx(1,ii0+2,jj,kk0,lb)+facevarx(1,ii0+2,jj1,kk0,lb)
     .        +facevarx(1,ii0+2,jj,kk1,lb)+facevarx(1,ii0+2,jj1,kk1,lb)
            bzsum =facevarz(1,ii0,jj,kk0,lb)+facevarz(1,ii0,jj1,kk0,lb)
     .            +facevarz(1,ii1,jj,kk0,lb)+facevarz(1,ii1,jj1,kk0,lb)
            bzsum1=
     .         facevarz(1,ii0,jj,kk0+2,lb)+facevarz(1,ii0,jj1,kk0+2,lb)
     .        +facevarz(1,ii1,jj,kk0+2,lb)+facevarz(1,ii1,jj1,kk0+2,lb)


            divbl  = (bxsum1 - bxsum + bysum - bysum1 + bzsum1 - bzsum)
     .               *deli*.125
            divbnl = (bxsum1 - bxsum + 4.*recvy(1,i2,j2,k2) -bysum1 + 
     .               bzsum1 - bzsum)*deli*.125


           iprint = iprint+1
           pattern_buf(iprint) = 5
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bysum
           data_buf(2,iprint) = recvy(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.5) then
            call shmem_real_get(recvz(1,1,1,1),
     .                          facevarz(1,1,1,1,neigh_blk),
     .                          len_blockfz*nbndvar,neigh_pe)
            k2 = nzb+nguard+1
            bzsum = facevarz(1,i0,j0,k0,lb)+facevarz(1,i1,j0,k0,lb)
     .             +facevarz(1,i0,j1,k0,lb)+facevarz(1,i1,j1,k0,lb)

            bzsum1= facevarz(1,i0,j0,k0+2,lb)+facevarz(1,i1,j0,k0+2,lb)
     .             +facevarz(1,i0,j1,k0+2,lb)+facevarz(1,i1,j1,k0+2,lb)
            kk = k0
            kk1 = kk+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            bysum =facevary(1,ii0,jj0,kk,lb)+facevary(1,ii0,jj0,kk1,lb)
     .            +facevary(1,ii1,jj0,kk,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1=
     .         facevary(1,ii0,jj0+2,kk,lb)+facevary(1,ii0,jj0+2,kk1,lb)
     .        +facevary(1,ii1,jj0+2,kk,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            bxsum =facevarx(1,ii0,jj0,kk,lb)+facevarx(1,ii0,jj0,kk1,lb)
     .            +facevarx(1,ii0,jj1,kk,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1=
     .         facevarx(1,ii0+2,jj0,kk,lb)+facevarx(1,ii0+2,jj0,kk1,lb)
     .        +facevarx(1,ii0+2,jj1,kk,lb)+facevarx(1,ii0+2,jj1,kk1,lb)

            divbl  = (bxsum1 - bxsum + bysum1 - bysum + bzsum1 - bzsum)
     .               *deli*.125
            divbnl = (bxsum1 - bxsum + bysum1 - bysum + bzsum1 - 
     .                  4.*recvz(1,i2,j2,k2))*deli*.125

           iprint = iprint+1
           pattern_buf(iprint) = 6
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bzsum
           data_buf(2,iprint) = recvz(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         elseif(itface(lb).eq.6) then
            call shmem_real_get(recvz(1,1,1,1),
     .                          facevarz(1,1,1,1,neigh_blk),
     .                          len_blockfz*nbndvar,neigh_pe)
            k2 = nguard+1
            bzsum = facevarz(1,i0,j0,k0,lb)+facevarz(1,i1,j0,k0,lb)
     .             +facevarz(1,i0,j1,k0,lb)+facevarz(1,i1,j1,k0,lb)

            bzsum1= facevarz(1,i0,j0,k0-2,lb)+facevarz(1,i1,j0,k0-2,lb)
     .             +facevarz(1,i0,j1,k0-2,lb)+facevarz(1,i1,j1,k0-2,lb)

            kk = k0-2
            kk1 = kk+1
            jj0 = min(j0,j1)
            jj1 = jj0+1
            ii0 = min(i0,i1)
            ii1 = ii0+1
            bysum =facevary(1,ii0,jj0,kk,lb)+facevary(1,ii0,jj0,kk1,lb)
     .            +facevary(1,ii1,jj0,kk,lb)+facevary(1,ii1,jj0,kk1,lb)
            bysum1=
     .         facevary(1,ii0,jj0+2,kk,lb)+facevary(1,ii0,jj0+2,kk1,lb)
     .        +facevary(1,ii1,jj0+2,kk,lb)+facevary(1,ii1,jj0+2,kk1,lb)
            bxsum =facevarx(1,ii0,jj0,kk,lb)+facevarx(1,ii0,jj0,kk1,lb)
     .            +facevarx(1,ii0,jj1,kk,lb)+facevarx(1,ii0,jj1,kk1,lb)
            bxsum1=
     .         facevarx(1,ii0+2,jj0,kk,lb)+facevarx(1,ii0+2,jj0,kk1,lb)
     .        +facevarx(1,ii0+2,jj1,kk,lb)+facevarx(1,ii0+2,jj1,kk1,lb)

            divbl  = (bxsum1 - bxsum + bysum1 - bysum + bzsum - bzsum1)
     .               *deli*.125
            divbnl = (bxsum1 - bxsum + bysum1 - bysum +  
     .                  4.*recvz(1,i2,j2,k2) - bzsum1)*deli*.125

           iprint = iprint+1
           pattern_buf(iprint) = 6
           lb_buf(iprint) = lb
           level_buf(iprint) = lrefine(lb)
           data_buf(1,iprint) = .25*bzsum
           data_buf(2,iprint) = recvz(1,i2,j2,k2)
           data_buf(3,iprint) = divbl
           data_buf(4,iprint) = divbnl

         endif


      endif                             ! end of neigh type iftest


      endif                             ! loop over block faces

      endif
      enddo                             ! end of loop over blocks
      endif

      call shmem_barrier_all()

!
! Is there any output to collect?
      jprint = 0
      call comm_int_sum_to_all(jprint,iprint)

      if(mype.eq.0.and.jprint.gt.0) then

        open(unit=22,status='unknown',position='append',
     .                                  file='facevar.dbg')
        do iproc=0,nprocs-1
          kprint = 0
          call shmem_integer_get(kprint,iprint,1,iproc)
          if(kprint.gt.0) then
            call shmem_integer_get(lb_buf,lb_buf,
     .                                        kprint,iproc)
            call shmem_integer_get(level_buf,level_buf,
     .                                        kprint,iproc)
            call shmem_integer_get(pattern_buf,pattern_buf,
     .                                        kprint,iproc)
            call shmem_real_get(data_buf(1,1),data_buf(1,1),
     .                          4*kprint,iproc)
            do ip = 1,kprint
              if(pattern_buf(ip).eq.1) then
                pattern = 'srl x '
              elseif(pattern_buf(ip).eq.2) then
                pattern = 'srl y '
              elseif(pattern_buf(ip).eq.3) then
                pattern = 'srl z '
              elseif(pattern_buf(ip).eq.4) then
                pattern = 'drl x '
              elseif(pattern_buf(ip).eq.5) then
                pattern = 'drl y '
              elseif(pattern_buf(ip).eq.6) then
                pattern = 'drl z '
              endif
              pattern1 = 'div B '
              write(22,700) pattern,istep,iproc,lb_buf(ip),
     .             x0,y0,z0,data_buf(1,ip),data_buf(2,ip),level_buf(ip)
              write(22,700) pattern1,istep,iproc,lb_buf(ip),
     .             x0,y0,z0,data_buf(3,ip),data_buf(4,ip),level_buf(ip)
            enddo
          endif
700      format(a6,i5,1x,i4,1x,i4,1x,5(1x,1pe12.5),1x,i4)

        enddo

        close(unit=22)

      endif


      call shmem_barrier_all()

      return
      end subroutine test_neigh_data
