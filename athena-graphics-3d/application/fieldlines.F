#include "paramesh_preprocessor.fh"

      subroutine fieldlines(x0,y0,z0,xf,yf,zf,
     .                      sign,istep_type,npts,mype,ldebug)

!
! This routine traces a fieldline begining at (x0,y0,z0).
!
!
! Arguments :
!  x0          real           starting x point 
!  y0          real           starting y point 
!  z0          real           starting z point 
!  xf          real array     x coords of fieldline points
!  yf          real array     y coords of fieldline points
!  zf          real array     z coords of fieldline points
!  sign        real           field direction to follow
!  istep_type  integer        algorithm to use in tracing fieldline
!                              1 = fixed stepsize
!                              2 = variable stepsize based on ??
!                              3 = Richardson extrapolation
!  npts        integer        no of points on this fieldline
!  mype        integer        local processor number
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */




      integer :: istep_type,npts
      logical :: ldebug

!---------------------------------------------------


#include "amr_shmem.fh"

#include "b_corner.common"


      integer, parameter :: nrich = 5
      real, dimension(nrich) :: x_rich,y_rich,z_rich,steps


#ifdef MPI_USED
      include "mpif.h"
!     integer :: errcode
#endif /* MPI_USED */

!---------------------------------------------------

! local amr variables
      integer nprocs
      integer shmem_my_pe,shmem_n_pes

      integer mype,errcode

      integer npt_max
c      parameter(npt_max = 1500)
      parameter(npt_max = 3500)

      real  xf(npt_max),yf(npt_max),zf(npt_max)
      real  bbox(2,3)
      real  bx_f(iu_bnd+1,ju_bnd+k2d,ku_bnd+k3d)
      real  by_f(iu_bnd+1,ju_bnd+k2d,ku_bnd+k3d)
      real  bz_f(iu_bnd+1,ju_bnd+k2d,ku_bnd+k3d)
      integer  blk_in_cache(2)

      logical lfound,lcutdx

      integer nblocks,tot_blocks,cnodetype,itot,pe,ib,lb
      save nblocks,tot_blocks,itot,bbox,cnodetype

      real, allocatable :: bounds(:,:,:)
      integer, allocatable :: block_add(:,:)

!---------------------------------------------------


      nprocs = shmem_n_pes()
      blk_in_cache(:) = -1

      stepsize0 = 100./100.
      stepsize0 = 200./100.
      stepsize = stepsize0



!---------------------------------------------------
!
! Step 1. Construct a list of the bounding boxes for all blocks.
!

      itot = 0
      do lb = 1,lnblocks
        if(nodetype(lb).eq.1) itot = itot+1
      enddo

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

! Find total number of leaf blocks
      call comm_int_sum_to_all(tot_blocks,itot)
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


      if(mype.eq.0) then

        allocate(bounds(2,3,tot_blocks))
        allocate(block_add(2,tot_blocks))
     
!
! construct block bounds for all blocks
      icount=0
      do iproc = 0,nprocs-1
        call shmem_integer_get(nblocks,lnblocks,1,iproc)
        if(nblocks.gt.0) then
          do lb=1,nblocks
            call shmem_integer_get(cnodetype,nodetype(lb),1,iproc)
            if(cnodetype.eq.1) then
              call shmem_real_get(bbox(1,1),bnd_box(1,1,lb),6,
     .                                                  iproc)
              icount = icount+1
              bounds(:,:,icount) = bbox(:,:)
              block_add(1,icount) = lb
              block_add(2,icount) = iproc
            endif
          enddo 
        endif
      enddo

!---------------------------------------------------
!
! Step 2. Convert starting point to (x,y,z).


! initial point
      xx = x0
      yy = y0
      zz = z0

      npt = 1
      bmax = 1.

      xf(npt) = xx
      yf(npt) = yy
      zf(npt) = zz

!      write(*,*) 'initial point ',x0,y0,z0

!---------------------------------------------------
      if(istep_type.eq.1) then
!---------------------------------------------------
!
! Use a constant stepsize when tracing B.
!
!
! Step 3. Trace fieldline.


      do while (bmax.gt.0.and.npt.lt.npt_max)


!-------
! Step 3.1  find bx,by,bz for current location
!-------


! first locate grid block containing (xx,yy,zz).
        lfound=.false.
        pe = 0
        ib = 0
        lb = 0
        do while (.not.lfound.and.lb.lt.tot_blocks)
          lb = lb+1

          xtest = (bounds(1,1,lb)-xx)*(bounds(2,1,lb)-xx)
          ytest = (bounds(1,2,lb)-yy)*(bounds(2,2,lb)-yy)
          if(ndim.gt.2) then
            ztest = (bounds(1,3,lb)-zz)*(bounds(2,3,lb)-zz)
          else
            ztest = -1.
          endif
          if(xtest.le.0..and.ytest.le.0..and.ztest.le.0.) then
            ib = block_add(1,lb)
            pe = block_add(2,lb)
            lfound = .true.
          endif
        enddo

        if(.not.lfound) then
!            write(*,*) 'fieldline search failed or line has exited'
!            write(*,*) 'line has ',npt,' points'
            npts = npt
            return
        endif


! fetch field for block lb on processor pe
! note : these fields are assumed to be known at the corners of cells
        if(ib.ne.blk_in_cache(1).or.pe.ne.blk_in_cache(2)) then

            call shmem_real_get(bx_f(1,1,1),bbx(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            call shmem_real_get(by_f(1,1,1),bby(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            call shmem_real_get(bz_f(1,1,1),bbz(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            blk_in_cache(1) = ib
            blk_in_cache(2) = pe

        endif


! find local field values
        dx = (bounds(2,1,lb)-bounds(1,1,lb))/real(nxb)
        dy = (bounds(2,2,lb)-bounds(1,2,lb))/real(nyb)
        dz = (bounds(2,3,lb)-bounds(1,3,lb))/real(nzb)
        if(ndim.gt.2) then
          dstep = min(dx,dy,dz)
        else
          dstep = min(dx,dy)
        endif
        stepsize = min(stepsize,.25*dstep)

        ii = 1
        do i=1,nxb
          xlo = bounds(1,1,lb)+real(i-1)*dx
          xhi = xlo+dx
          if(i.eq.1) xlo = bounds(1,1,lb)
          if(i.eq.nxb) xhi = bounds(2,1,lb)
          xtest = (xhi-xx)*(xlo-xx)
          if(xtest.le.0.) then
            ii = i
            x1 = xlo
            x2 = xhi
          endif
        enddo
        jj = 1
        do j=1,nyb
          ylo = bounds(1,2,lb)+real(j-1)*dy
          yhi = ylo+dy
          if(j.eq.1) ylo = bounds(1,2,lb)
          if(j.eq.nyb) yhi = bounds(2,2,lb)
          ytest = (yhi-yy)*(ylo-yy)
          if(ytest.le.0.) then
            jj = j
            y1 = ylo
            y2 = yhi
          endif
        enddo
        kk = 1
        do k=1,nzb
          zlo = bounds(1,3,lb)+real(k-1)*dz
          zhi = zlo+dz
          if(k.eq.1) zlo = bounds(1,3,lb)
          if(k.eq.nzb) zhi = bounds(2,3,lb)
          ztest = (zhi-zz)*(zlo-zz)
          if(ztest.le.0.) then
            kk = k
            z1 = zlo
            z2 = zhi
          endif
        enddo
        if(ndim.eq.2) then
          kk=1
          z1 = zz
          z2 = z1
        endif

        dx1 = abs((xx-x1)/(x2-x1))
        dx2 = 1.-dx1
        dy1 = abs((yy-y1)/(y2-y1))
        dy2 = 1.-dy1
        if(ndim.eq.3) then
          dz1 = abs((zz-z1)/(z2-z1))
        else
          dz1 = 1.
        endif
        dz2 = 1.-dz1


!-------
! Step 3.2  find bx,by,bz for current location
!-------


        bx = dz1*(
     . +dy1*(dx1*bx_f(ii+1,jj+k2d,kk+k3d)+dx2*bx_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*bx_f(ii+1,jj,kk+k3d)+dx2*bx_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*bx_f(ii+1,jj+k2d,kk)+dx2*bx_f(ii,jj+k2d,kk))
     . +dy2*(dx1*bx_f(ii+1,jj,kk)+dx2*bx_f(ii,jj,kk))
     .           )
        by = dz1*(
     . +dy1*(dx1*by_f(ii+1,jj+k2d,kk+k3d)+dx2*by_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*by_f(ii+1,jj,kk+k3d)+dx2*by_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*by_f(ii+1,jj+k2d,kk)+dx2*by_f(ii,jj+k2d,kk))
     . +dy2*(dx1*by_f(ii+1,jj,kk)+dx2*by_f(ii,jj,kk))
     .           )
        if(ndim.gt.2) then
        bz = dz1*(
     . +dy1*(dx1*bz_f(ii+1,jj+k2d,kk+k3d)+dx2*bz_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*bz_f(ii+1,jj,kk+k3d)+dx2*bz_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*bz_f(ii+1,jj+k2d,kk)+dx2*bz_f(ii,jj+k2d,kk))
     . +dy2*(dx1*bz_f(ii+1,jj,kk)+dx2*bz_f(ii,jj,kk))
     .           )
        else
          bz = 0.
        endif
          

        bmax = max(abs(bx),abs(by),abs(bz))

!-------
! Step 3.3  compute next point (x1,y1,z1) along fieldline
!-------
! Is field non-zero ?
! If  not compute next fieldline element.
        if(bmax.gt.0) then
          dx = bx/bmax
          dy = by/bmax
          dz = bz/bmax
          dnorm = sqrt(dx*dx+dy*dy+dz*dz)
          dx = sign*stepsize*dx/dnorm
          dy = sign*stepsize*dy/dnorm
          dz = sign*stepsize*dz/dnorm
          x1 = xf(npt) + dx
          y1 = yf(npt) + dy
          z1 = zf(npt) + dz


! compute next point along fieldline
          xf(npt+1) = x1
          yf(npt+1) = y1
          zf(npt+1) = z1
          npt = npt+1

!-------
! Step 3.4  update current point.
!-------
          xx = x1
          yy = y1
          zz = z1

!-------

        endif

      enddo


!---------------------------------------------------
      elseif(istep_type.eq.2) then
!---------------------------------------------------
!
!  Control stepsize using change in B.
!

! Step 3. Trace fieldline.

      iloop = 0
      do while (bmax.gt.0.and.npt.lt.npt_max)


!-------
! Step 3.1  find bx,by,bz for current location
!-------


! first locate grid block containing (xx,yy,zz).
        lfound=.false.
        pe = 0
        ib = 0
        lb = 0
        do while (.not.lfound.and.lb.lt.tot_blocks)
          lb = lb+1

          xtest = (bounds(1,1,lb)-xx)*(bounds(2,1,lb)-xx)
          ytest = (bounds(1,2,lb)-yy)*(bounds(2,2,lb)-yy)
          if(ndim.gt.2) then
            ztest = (bounds(1,3,lb)-zz)*(bounds(2,3,lb)-zz)
          else
            ztest = -1.
          endif
          if(xtest.le.0..and.ytest.le.0..and.ztest.le.0.) then
            ib = block_add(1,lb)
            pe = block_add(2,lb)
            lfound = .true.
          endif
        enddo

        if(.not.lfound) then
!            write(*,*) 'fieldline search failed or line has exited'
!            write(*,*) 'at location ',xx,yy,zz
            write(*,*) 'line has ',npt,' points'
            npts = npt
            return
        endif

!        write(*,*) 'pt in block ',ib,pe

! fetch field for block lb on processor pe
! note : these fields are assumed to be known at the corners of cells
        if(ib.ne.blk_in_cache(1).or.pe.ne.blk_in_cache(2)) then

            call shmem_real_get(bx_f(1,1,1),bbx(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            call shmem_real_get(by_f(1,1,1),bby(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            call shmem_real_get(bz_f(1,1,1),bbz(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            blk_in_cache(1) = ib
            blk_in_cache(2) = pe

        endif

! find local field values
        dx = (bounds(2,1,lb)-bounds(1,1,lb))/real(nxb)
        dy = (bounds(2,2,lb)-bounds(1,2,lb))/real(nyb)
        dz = (bounds(2,3,lb)-bounds(1,3,lb))/real(nzb)
        dstep = min(dx,dy,dz)
        stepsize = min(stepsize,.25*dstep)

        ii = 1
        do i=1,nxb
          xlo = bounds(1,1,lb)+real(i-1)*dx
          xhi = xlo+dx
          if(i.eq.1) xlo = bounds(1,1,lb)
          if(i.eq.nxb) xhi = bounds(2,1,lb)
          xtest = (xhi-xx)*(xlo-xx)
          if(xtest.le.0.) then
            ii = i
            x1 = xlo
            x2 = xhi
          endif
        enddo
        jj = 1
        do j=1,nyb
          ylo = bounds(1,2,lb)+real(j-1)*dy
          yhi = ylo+dy
          if(j.eq.1) ylo = bounds(1,2,lb)
          if(j.eq.nyb) yhi = bounds(2,2,lb)
          ytest = (yhi-yy)*(ylo-yy)
          if(ytest.le.0.) then
            jj = j
            y1 = ylo
            y2 = yhi
          endif
        enddo
        kk = 1
        do k=1,nzb
          zlo = bounds(1,3,lb)+real(k-1)*dz
          zhi = zlo+dz
          if(k.eq.1) zlo = bounds(1,3,lb)
          if(k.eq.nzb) zhi = bounds(2,3,lb)
          ztest = (zhi-zz)*(zlo-zz)
          if(ztest.le.0.) then
            kk = k
            z1 = zlo
            z2 = zhi
          endif
        enddo
        if(ndim.eq.2) then
          kk=1
          z1 = zz
          z2 = z1
        endif

        dx1 = abs((xx-x1)/(x2-x1))
        dx2 = 1.-dx1
        dy1 = abs((yy-y1)/(y2-y1))
        dy2 = 1.-dy1
        if(ndim.eq.3) then
          dz1 = abs((zz-z1)/(z2-z1))
        else
          dz1 = 1.
        endif
        dz2 = 1.-dz1

!        write(*,*) 'pt in cell ',ii,jj,kk

!-------
! Step 3.2  find bx,by,bz for current location
!-------

        bx = dz1*(
     . +dy1*(dx1*bx_f(ii+1,jj+k2d,kk+k3d)+dx2*bx_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*bx_f(ii+1,jj,kk+k3d)+dx2*bx_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*bx_f(ii+1,jj+k2d,kk)+dx2*bx_f(ii,jj+k2d,kk))
     . +dy2*(dx1*bx_f(ii+1,jj,kk)+dx2*bx_f(ii,jj,kk))
     .           )
        by = dz1*(
     . +dy1*(dx1*by_f(ii+1,jj+k2d,kk+k3d)+dx2*by_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*by_f(ii+1,jj,kk+k3d)+dx2*by_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*by_f(ii+1,jj+k2d,kk)+dx2*by_f(ii,jj+k2d,kk))
     . +dy2*(dx1*by_f(ii+1,jj,kk)+dx2*by_f(ii,jj,kk))
     .           )
        if(ndim.gt.2) then
          bz = dz1*(
     . +dy1*(dx1*bz_f(ii+1,jj+k2d,kk+k3d)+dx2*bz_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*bz_f(ii+1,jj,kk+k3d)+dx2*bz_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*bz_f(ii+1,jj+k2d,kk)+dx2*bz_f(ii,jj+k2d,kk))
     . +dy2*(dx1*bz_f(ii+1,jj,kk)+dx2*bz_f(ii,jj,kk))
     .           )
        else
          bz = 0.
        endif


          
! Compute field components in x,y,z coords.
        bx_new = bx
        by_new = by
        bz_new = bz

!        write(*,*) 'est new field ',bx_new,by_new,bz_new


        if(iloop.eq.0) then
          bx_old = bx_new
          by_old = by_new
          bz_old = bz_new
        endif

!-------
! Step 3.3  compute next point (xnew,ynew,znew) along fieldline
!-------

! Is field non-zero ?
! If  not compute next fieldline element.
        bmax = max(abs(bx_new),abs(by_new),abs(bz_new))
        if(bmax.gt.0) then

        if(iloop.eq.0) then
          dx = bx_new/bmax
          dy = by_new/bmax
          dz = bz_new/bmax
          dnorm = sqrt(dx*dx+dy*dy+dz*dz)
          dx = sign*stepsize*dx/dnorm
          dy = sign*stepsize*dy/dnorm
          dz = sign*stepsize*dz/dnorm
          xnew = xf(npt) + dx
          ynew = yf(npt) + dy
          znew = zf(npt) + dz

        elseif(iloop.gt.0) then

! Is dB too large?
          lcutdx = .false.
          bfactor = abs(
     .          (bx_new*bx_old+by_new*by_old+bz_new*bz_old)
     .         /(bx_old*bx_old+by_old*by_old+bz_old*bz_old) - 1.)
!          write(*,*) 'bfactor = ',bfactor
          if(bfactor.gt..1) lcutdx = .true.

          if(lcutdx) then
!            write(*,*) 'bfactor = ',bfactor
            stepsize = .9*stepsize
            write(*,*) 'new stepsize = ',stepsize
            dx = bx_new/bmax
            dy = by_new/bmax
            dz = bz_new/bmax
            dnorm = sqrt(dx*dx+dy*dy+dz*dz)
            dx = sign*stepsize*dx/dnorm
            dy = sign*stepsize*dy/dnorm
            dz = sign*stepsize*dz/dnorm
            xnew = xf(npt) + dx
            ynew = yf(npt) + dy
            znew = zf(npt) + dz
          else
! save this point
            xf(npt+1) = xnew
            yf(npt+1) = ynew
            zf(npt+1) = znew
            if(ldebug)
     .       write(*,*) 'save pt x y z npt+1 ',xnew,ynew,znew,npt+1
! save this field
            bx_old = bx_new
            by_old = by_new
            bz_old = bz_new
! take a new step
            dx = bx_new/bmax
            dy = by_new/bmax
            dz = bz_new/bmax
            dnorm = sqrt(dx*dx+dy*dy+dz*dz)
            dx = sign*stepsize*dx/dnorm
            dy = sign*stepsize*dy/dnorm
            dz = sign*stepsize*dz/dnorm
!            write(*,*) 'new dx dy dz ',dx,dy,dz
! compute the new location
            xnew = xf(npt+1) + dx
            ynew = yf(npt+1) + dy
            znew = zf(npt+1) + dz
!            write(*,*) 'new x y z ',xnew,ynew,znew

            npt = npt+1
          endif

        endif                    ! end of iloop if test


        endif                    ! end of bmax if test

!-------
! Step 3.4  update current point.
!-------
        xx = xnew
        yy = ynew
        zz = znew
!-------


      iloop = iloop + 1

      enddo


!---------------------------------------------------
      elseif(istep_type.eq.3) then
!---------------------------------------------------
!
! Use Richardson extrapolation to trace field.
!
!

! Step 3. Trace fieldline.

      do while (bmax.gt.0.and.npt.lt.npt_max)


!-------
! Step 3.1  find bx,by,bz for current location
!-------

! save coords before begining Richardson extrapolation for this step
        x_old = xx
        y_old = yy
        z_old = zz
        xnew = xf(npt)
        ynew = yf(npt)
        znew = zf(npt)


! evaluate new point on fieldline nrich times using different step sizes
        do nr = 1,nrich

        xnew = xf(npt)
        ynew = yf(npt)
        znew = zf(npt)
!        write(*,*) 'new x y z ',xnew,ynew,znew

        xx = x_old
        yy = y_old
        zz = z_old
        stepsize = stepsize0/2**(nrich-nr) 

        do iter = 1,2**(nrich-nr)

!        write(*,*) 'nr iter ',nr,iter



! first locate grid block containing (xx,yy,zz).
        lfound=.false.
        pe = 0
        ib = 0
        lb = 0
        do while (.not.lfound.and.lb.lt.tot_blocks)
          lb = lb+1

          xtest = (bounds(1,1,lb)-xx)*(bounds(2,1,lb)-xx)
          ytest = (bounds(1,2,lb)-yy)*(bounds(2,2,lb)-yy)
          if(ndim.gt.2) then
            ztest = (bounds(1,3,lb)-zz)*(bounds(2,3,lb)-zz)
          else
            ztest = -1.
          endif
          if(xtest.le.0..and.ytest.le.0..and.ztest.le.0.) then
            ib = block_add(1,lb)
            pe = block_add(2,lb)
            lfound = .true.
          endif
        enddo

        if(.not.lfound) then
!            write(*,*) 'fieldline search failed or line has exited'
!            write(*,*) 'at location ',xx,yy,zz
!            write(*,*) 'line has ',npt,' points'
            npts = npt
            return
        endif

!        write(*,*) 'pt in block ',ib,pe

! fetch field for block lb on processor pe
! note : these fields are assumed to be known at the corners of cells
        if(ib.ne.blk_in_cache(1).or.pe.ne.blk_in_cache(2)) then

            call shmem_real_get(bx_f(1,1,1),bbx(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            call shmem_real_get(by_f(1,1,1),bby(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            call shmem_real_get(bz_f(1,1,1),bbz(1,1,1,ib),
     .            (iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d),pe)
            blk_in_cache(1) = ib
            blk_in_cache(2) = pe

        endif

! find local field values
        dx = (bounds(2,1,lb)-bounds(1,1,lb))/real(nxb)
        dy = (bounds(2,2,lb)-bounds(1,2,lb))/real(nyb)
        dz = (bounds(2,3,lb)-bounds(1,3,lb))/real(nzb)

        ii = 1
        do i=1,nxb
          xlo = bounds(1,1,lb)+real(i-1)*dx
          xhi = xlo+dx
          if(i.eq.1) xlo = bounds(1,1,lb)
          if(i.eq.nxb) xhi = bounds(2,1,lb)
          xtest = (xhi-xx)*(xlo-xx)
          if(xtest.le.0.) then
            ii = i
            x1 = xlo
            x2 = xhi
          endif
        enddo
        jj = 1
        do j=1,nyb
          ylo = bounds(1,2,lb)+real(j-1)*dy
          yhi = ylo+dy
          if(j.eq.1) ylo = bounds(1,2,lb)
          if(j.eq.nyb) yhi = bounds(2,2,lb)
          ytest = (yhi-yy)*(ylo-yy)
          if(ytest.le.0.) then
            jj = j
            y1 = ylo
            y2 = yhi
          endif
        enddo
        kk = 1
        do k=1,nzb
          zlo = bounds(1,3,lb)+real(k-1)*dz
          zhi = zlo+dz
          if(k.eq.1) zlo = bounds(1,3,lb)
          if(k.eq.nzb) zhi = bounds(2,3,lb)
          ztest = (zhi-zz)*(zlo-zz)
          if(ztest.le.0.) then
            kk = k
            z1 = zlo
            z2 = zhi
          endif
        enddo
        if(ndim.eq.2) then
          kk=1
          z1 = zz
          z2 = z1
        endif

        dx1 = abs((xx-x1)/(x2-x1))
        dx2 = 1.-dx1
        dy1 = abs((yy-y1)/(y2-y1))
        dy2 = 1.-dy1
        if(ndim.eq.3) then
          dz1 = abs((zz-z1)/(z2-z1))
        else
          dz1 = 1.
        endif
        dz2 = 1.-dz1

!        write(*,*) 'pt in cell ',ii,jj,kk

!-------
! Step 3.2  find bx,by,bz for current location
!-------

        bx = dz1*(
     . +dy1*(dx1*bx_f(ii+1,jj+k2d,kk+k3d)+dx2*bx_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*bx_f(ii+1,jj,kk+k3d)+dx2*bx_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*bx_f(ii+1,jj+k2d,kk)+dx2*bx_f(ii,jj+k2d,kk))
     . +dy2*(dx1*bx_f(ii+1,jj,kk)+dx2*bx_f(ii,jj,kk))
     .           )
        by = dz1*(
     . +dy1*(dx1*by_f(ii+1,jj+k2d,kk+k3d)+dx2*by_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*by_f(ii+1,jj,kk+k3d)+dx2*by_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*by_f(ii+1,jj+k2d,kk)+dx2*by_f(ii,jj+k2d,kk))
     . +dy2*(dx1*by_f(ii+1,jj,kk)+dx2*by_f(ii,jj,kk))
     .           )
        if(ndim.gt.2) then
        bz = dz1*(
     . +dy1*(dx1*bz_f(ii+1,jj+k2d,kk+k3d)+dx2*bz_f(ii,jj+k2d,kk+k3d))
     . +dy2*(dx1*bz_f(ii+1,jj,kk+k3d)+dx2*bz_f(ii,jj,kk+k3d))
     .           ) + dz2*(
     . +dy1*(dx1*bz_f(ii+1,jj+k2d,kk)+dx2*bz_f(ii,jj+k2d,kk))
     . +dy2*(dx1*bz_f(ii+1,jj,kk)+dx2*bz_f(ii,jj,kk))
     .           )
        else
          bz = 0.
        endif

          
!        write(*,*) 'est new field ',bx,by,bz


!-------
! Step 3.3  compute next point (xnew,ynew,znew) along fieldline
!-------

! Is field non-zero ?
! If  not compute next fieldline element.
        bmax = max(abs(bx),abs(by),abs(bz))
        if(bmax.gt.0) then

          dx = bx/bmax
          dy = by/bmax
          dz = bz/bmax
          dnorm = sqrt(dx*dx+dy*dy+dz*dz)
          dx = sign*stepsize*dx/dnorm
          dy = sign*stepsize*dy/dnorm
          dz = sign*stepsize*dz/dnorm
          xnew = xnew + dx
          ynew = ynew + dy
          znew = znew + dz

        endif                    ! end of bmax if test

!-------
! Step 3.4  update current point.
!-------
        xx = xnew
        yy = ynew
        zz = znew

!-------
!      write(*,*) 'new x y z ',xnew,ynew,znew

      enddo                      ! end of iter do loop

! store results for this step size
      x_rich(nrich-nr+1) = xnew
      y_rich(nrich-nr+1) = ynew
      z_rich(nrich-nr+1) = znew
      steps(nrich-nr+1)  = stepsize
!      write(*,*) 'step x y z ',stepsize,xnew,ynew,znew

      enddo                      ! end of nr do loop

! extrapolate results to zero stepsize
      call polint(steps,x_rich,nrich,0.,x_int,errorx)
      call polint(steps,y_rich,nrich,0.,y_int,errory)
      call polint(steps,z_rich,nrich,0.,z_int,errorz)

      xf(npt+1) = x_int
      yf(npt+1) = y_int
      zf(npt+1) = z_int
      npt = npt+1
!      write(*,*) 'intr x y z ',x_int,y_int,z_int

      enddo




!---------------------------------------------------
      endif
!---------------------------------------------------


      deallocate( bounds )
      deallocate( block_add )

      endif                      ! end of mype=0 if test

      npts = npt

!     write(*,*) 'line has ',npts,' points'

      return
      end
