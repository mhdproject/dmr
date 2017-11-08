!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#define DEBUG_MORTON



      integer mdim,mchild
      parameter(mdim = 2)
      parameter(mchild = 2**mdim)

      integer integer :: 
     .  mort_child(mchild),mort_neigh(3,3,3)
      integer :: mort_p
      real    coord(3),bbsize(3)

      nchild = mchild
      ndim = mdim
      lrefine_max = 3

      mort_p = -1

      lrefine = 3
      mort = 4
      mort_p = morton_parent(mort,lrefine,lrefine_max,ndim)
      write(*,*) 'parent of ',mort,' is ',mort_p

      lrefine = 3
      mort =5 
      mort_p = morton_parent(mort,lrefine,lrefine_max,ndim)
      write(*,*) 'parent of ',mort,' is ',mort_p

      lrefine = 3
      mort = 6
      mort_p = morton_parent(mort,lrefine,lrefine_max,ndim)
      write(*,*) 'parent of ',mort,' is ',mort_p

      lrefine = 2
      mort = 12
      mort_p = morton_parent(mort,lrefine,lrefine_max,ndim)
      write(*,*) 'parent of ',mort,' is ',mort_p


      mort = 0
      lrefine = 1
      call morton_children(mort,lrefine,lrefine_max,ndim,
     .                           mort_child,nchild)
      write(*,*) 'children ',mort_child

      mort = 4
      lrefine = 2
      call morton_children(mort,lrefine,lrefine_max,ndim,
     .                           mort_child,nchild)
      write(*,*) 'children ',mort_child


      xmin = 0.
      ymin = 0.
      zmin = 0.
      coord(1) = .625
      coord(2) = .375
      coord(3) = .0
      bbsize(1) = .25
      bbsize(2) = .25
      bbsize(3) = 1.0
      lrefine = 3
      lrefine_max = 3
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                            .true.,.true.,.true.,
     .                            coord,size,ndim,
     .                            lrefine,lrefine_max,mort_neigh)
      do k = 1,3
      do j = 1,3
      write(*,*) mort_neigh(:,j,k)
      enddo
      enddo

      stop
      end




      integer 
     .    function morton_parent(mort,lrefine,lrefine_max,ndim)

!---------------------------------------------------------------------
!
! This function returns the morton number for the parent of the
! tree node MORT.
!
! Written :     April 1999     Peter MacNeice
!---------------------------------------------------------------------

      integer lrefine,lrefine_max,ndim,nbits_zero
      integer :: mort
      integer :: mort_p
      integer :: lbit

      nbits_zero = (lrefine_max+1-lrefine)*ndim - 1

      mort_p = mort
      do lbit = 0,nbits_zero
       mort_p = ibclr(mort_p,lbit)
      enddo

      morton_parent = mort_p

      return
      end




      subroutine morton_children(mort,lrefine,lrefine_max,ndim,
     .                           mort_child,nchild)
!---------------------------------------------------------------------
! Computes the morton numbers for children of a node with given
! morton number, refinement level, and where the maximum refinement
! level in use is known.
!
! Arguments:
! mort          integer        morton number of current node
! lrefine       integer        refinement level of current node
! lrefine_max   integer        highest refinement level currently in use
! ndim          integer        physical dimension
! mort_child    integer array  array of morton numbers of this nodes 
!                                children
! nchild        integer        number of children of any node
! 
! Written :     April 1999     Peter MacNeice
!
!---------------------------------------------------------------------

      integer mort,lrefine,lrefine_max,ndim
      integer lbit,nchild,ichild

      integer :: mort_child(nchild)
      integer :: mort
      integer :: i

      if(lrefine.lt.lrefine_max) then
        lbit = (lrefine_max-lrefine-1)*ndim
        do ichild = 0,nchild-1
          i = ishft(ichild,lbit)
          mort_child(ichild+1) = ior(mort,i)
        enddo
      endif

      return
      end




      subroutine morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                            lperiodicx,lperiodicy,lperiodz,
     .                            coord,size,ndim,
     .                            lrefine,lrefine_max,mort_neigh)

!---------------------------------------------------------------------
!
! This routine computes the morton numbers which all possible
! neighbors of the current block, at the same refinement level
! would have. If a neighbor would lie below the minimun coordinate
! value for any axis, the morton number returned defaults to -1.
!
! Arguments:
! xmin          real           lower x range of computational domain
! ymin          real           lower y range of computational domain
! zmin          real           lower z range of computational domain
! xmax          real           upper x range of computational domain
! ymax          real           upper y range of computational domain
! zmax          real           upper z range of computational domain
! lperiodicx    logical        .true. if periodic in x-direction, otherwise
!                                set to .false.
! lperiodicy    logical        .true. if periodic in y-direction, otherwise
!                                set to .false.
! lperiodicz    logical        .true. if periodic in z-direction, otherwise
!                                set to .false.
! coord         real array     array storing coordinates of current block
! size          real array     array storing size of current block
! ndim          integer        physical dimension
! lrefine       integer        refinement level of current block
! lrefine_max   integer        highest refinement level currently in use
!                                in the calculation
! mort_neigh    integer array  array storing morton numbers of all
!                                this blocks possible neighbors at its
!                                own refinement level.
!
! Written :     April 1999     Peter MacNeice
!
!---------------------------------------------------------------------

      real    xmin,ymin,zmin,xmax,ymax,zmax,coord(3),bbsize(3)
      logical lperiodicx,lperiodicy,lperiodicz
      integer ndim,mort_neigh(3,3,3)
      integer lrefine,lrefine_max
      integer i,j,k,nbits,ipos,ipos2
      real    x0,y0,z0,xtest,ytest,ztest,test
      integer :: ix,iy,iz
      integer :: mort_neigh(3,3,3)
      integer :: mort

#ifdef DEBUG_MORTON
      integer inxt_bits(3),inxt_bit,nbitshft
#endif

      nbits = bit_size(mort)
      mort_neigh(:,:,:) = -1
      k2d = 1
      k3d = 1

      if(ndim.eq.2) then
          k3d = 0
          bbsize(3) = 1.
      elseif(ndim.eq.1) then
          k2d = 0
          bbsize(2) = 1.
      endif

! cycle over possible neighbors
      do k = -k3d,k3d
        z0 = coord(3) + real(k)*bbsize(3) - zmin 
        if(lperiodicz) z0 = mod(z0,zmax-zmin)
        do j = -k2d,k2d
          y0 = coord(2) + real(j)*bbsize(2) - ymin
          if(lperiodicy) y0 = mod(y0,ymax-ymin)
          do i = -1,1
            x0 = coord(1) + real(i)*bbsize(1) - xmin
            if(lperiodicx) x0 = mod(x0,xmax-xmin)

! compute morton number
            ix = int(x0/bbsize(1))
            if (ndim.ge.2) then
               iy = int(y0/bbsize(2))
            else
               iy = 0
            end if
            if (ndim.eq.3) then
               iz = int(z0/bbsize(3))
            else
               iz = 0
            end if

            ipos  = 0
            ipos2 = 0
            mort = 0
            do while (ipos.lt. nbits-1-3)
              call mvbits (ix,ipos2,1,mort,ipos)
              ipos = ipos + 1
              if(ndim.ge.2) then
                 call mvbits (iy,ipos2,1,mort,ipos)
                 ipos = ipos + 1
              endif
              if(ndim.eq.3) then
                 call mvbits (iz,ipos2,1,mort,ipos)
                 ipos = ipos + 1
              endif
              ipos2 = ipos2 + 1
            end do
#ifdef DEBUG_MORTON
            inxt_bits(1) = ibits(ix,ipos2,1)
            inxt_bits(2) = ibits(iy,ipos2,1)
            inxt_bits(3) = ibits(iz,ipos2,1)
            inxt_bit = maxval(inxt_bits)
            if(inxt_bit.gt.0) then
                write(*,*) 'Error 1 in amr_morton_neighbor'
                call amr_abort()
            endif
#endif

            mort = ishft(mort,ndim*(lrefine_max-lrefine))
            
#ifdef DEBUG_MORTON
            nbitshft = ndim*(lrefine_max-lrefine)
            if(nbitshft.gt.0) then
            inxt_bit = ibits(mort,nbits-1-(nbitshft-1),
     .                       nbitshft-1)
            if(inxt_bit.gt.0) then
                write(*,*) 'Error 2 in morton_neighbor'
                call amr_abort()
            endif
            endif
#endif

            mort_neigh(i+2,j+2,k+2) = mort

! if coord is below minimum coords set mort to -1 
            xtest = x0*(x0-xmax)
            ytest = y0*(y0-ymax)
            ztest = z0*(z0-zmax)
            test = max(xtest,ytest)
            test = max( test,ztest)
            
            if(test.gt.0.) mort_neigh(i+2,j+2,k+2) = -1

          enddo
        enddo
      enddo

      return
      end
