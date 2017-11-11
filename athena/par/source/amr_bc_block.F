!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_bc_block(jface,ibc,iopt,l,mype)





!------------------------------------------------------------------------

! This routine sets the guard cell elements of the solution arrays 
! or workspace array on face jface of block l, using a boundary
! condition algorithm which the user must insert. Places where user
! editing is required are clearly indicated.
!
! Note, if your algorithm uses diagonal elements (see the users guide for
! a definition of diagonal elements) then you are responsible for setting
! the data values in these diagonal guard cells on the corners and edges
! of external non-periodic boundaries.
!
!
! Written :     Peter MacNeice          September 1997
!------------------------------------------------------------------------
!
! Arguments:
!       jface           number designating selected block face
!       ibc             the boundary condition flag. This should be
!                       set equal to neigh(1,jface,l) by the calling
!                       routine. Its value can be used to vary the
!                       boundary conditions applied on different boundaries.
!       iopt            a switch to control which data is updated
!                       iopt=1 will use 'unk' and/or 'facevarx(y)(z)'
!                       and/or 'unk_e_x(y)(z)' and or 'unk_n'.
!                       iopt=2 will use 'work'
!       l               number designating selected block 
!       mype            local processor number
!
!------------------------------------



! include file to define physical qualities of the model and mesh

      use paramesh_dimensions
      use physicaldata
      use workspace

! include file defining the tree
      use tree

      use paramesh_interfaces, only :  amr_1blk_bcset


! include file required for shmem library.
#include "amr_shmem.fh"
      include 'shmem_reduce.fh'

      integer, intent(in) :: jface,ibc,iopt,l,mype
 
      integer :: surrblks(3,3,3,3)

!------------------------------------

#ifndef NO_PERMANENT_GUARDCELLS 


      ip1 = 0
      jp1 = 0
      kp1 = 0

      surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) = 
     .         surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,l)

      if(iopt.eq.1) then

      if(jface.eq.2) ip1 = 1
      if(jface.eq.4) jp1 = 1
      if(jface.eq.6) kp1 = 1

! The loop bounds set here include the diagonal elements. If these are
! not required then these bounds can be changed to ilbnd=1+nguard,
! iubnd=nxb+nguard, and similarily for y and z.
                ilbnd = 1
                iubnd = nxb+2*nguard
                jlbnd = 1
                jubnd = nyb+2*nguard*k2d
                klbnd = 1
                kubnd = nzb+2*nguard*k3d

                if(jface.eq.1) iubnd = nguard
                if(jface.eq.2) ilbnd = 1+nxb+nguard
                if(jface.eq.3) jubnd = nguard
                if(jface.eq.4) jlbnd = 1+nyb+nguard
                if(jface.eq.5) kubnd = nguard
                if(jface.eq.6) klbnd = 1+nzb+nguard


                ibnd = 0
                jbnd = 0
                kbnd = 0
                if(jface.eq.1) ibnd = -1
                if(jface.eq.2) ibnd = 1
                if(jface.eq.3) jbnd = -1
                if(jface.eq.4) jbnd = 1
                if(jface.eq.5) kbnd = -1
                if(jface.eq.6) kbnd = 1

! Limit the application of this routine to leaf blocks
! or the parents of leaf blocks.
      if(nodetype(l).eq.1.or.nodetype(l).eq.2) then


       idest = 1
       if(nvar.gt.0) unk1(:,:,:,:,idest) = unk(:,:,:,:,l)
       if(nfacevar.gt.0) then
         facevarx1(:,:,:,:,idest) = facevarx(:,:,:,:,l)
         if(ndim.gt.1) 
     .   facevary1(:,:,:,:,idest) = facevary(:,:,:,:,l)
         if(ndim.eq.3) 
     .   facevarz1(:,:,:,:,idest) = facevarz(:,:,:,:,l)
       endif
       if(nvaredge.gt.0) then
         unk_e_x1(:,:,:,:,idest) = unk_e_x(:,:,:,:,l)
         unk_e_y1(:,:,:,:,idest) = unk_e_y(:,:,:,:,l)
         unk_e_z1(:,:,:,:,idest) = unk_e_z(:,:,:,:,l)
       endif
       if(nvarcorn.gt.0) unk_n1(:,:,:,:,idest) = 
     .                                unk_n(:,:,:,:,l)



      ilays = iubnd - ilbnd + 1
      jlays = jubnd - jlbnd + 1
      klays = kubnd - klbnd + 1
      call amr_1blk_bcset(mype,ibc,l,mype,
     .                    idest,iopt,ibnd,jbnd,kbnd,surrblks)

       if(nvar.gt.0) unk(:,:,:,:,l) = unk1(:,:,:,:,idest)
       if(nfacevar.gt.0) then
         facevarx(:,:,:,:,l) = facevarx1(:,:,:,:,idest)
         if(ndim.gt.1) 
     .   facevary(:,:,:,:,l) = facevary1(:,:,:,:,idest)
         if(ndim.eq.3) 
     .   facevarz(:,:,:,:,l) = facevarz1(:,:,:,:,idest)
       endif
       if(nvaredge.gt.0) then
         unk_e_x(:,:,:,:,l) = unk_e_x1(:,:,:,:,idest)
         unk_e_y(:,:,:,:,l) = unk_e_y1(:,:,:,:,idest)
         unk_e_z(:,:,:,:,l) = unk_e_z1(:,:,:,:,idest)
       endif
       if(nvarcorn.gt.0) unk_n(:,:,:,:,l) = 
     .                                unk_n1(:,:,:,:,idest)


      endif

      elseif(iopt.eq.2) then



! The loop bounds set here include the diagonal elements. If these are
! not required then these bounds can be changed to ilbnd=1+nguard_work,
! iubnd=nxb+nguard_work, and similarily for y and z.
                ilbnd = 1
                iubnd = nxb+2*nguard_work
                jlbnd = 1
                jubnd = nyb+2*nguard_work*k2d
                klbnd = 1
                kubnd = nzb+2*nguard_work*k3d

                if(jface.eq.1) iubnd = nguard_work
                if(jface.eq.2) ilbnd = 1+nxb+nguard_work
                if(jface.eq.3) jubnd = nguard_work
                if(jface.eq.4) jlbnd = 1+nyb+nguard_work
                if(jface.eq.5) kubnd = nguard_work
                if(jface.eq.6) klbnd = 1+nzb+nguard_work

! Limit the application of this routine to leaf blocks
! or the parents of leaf blocks.
       if(nodetype(l).eq.1.or.nodetype(l).eq.2) then

       idest = 1
       work1(:,:,:,idest) = work(:,:,:,l,iopt-1)

      ilays = iubnd - ilbnd + 1
      jlays = jubnd - jlbnd + 1
      klays = kubnd - klbnd + 1
      call amr_1blk_bcset(mype,ibc,l,mype,
     .                    idest,iopt,ibnd,jbnd,kbnd,surrblks)

       work(:,:,:,l,iopt-1) = work1(:,:,:,idest)

       endif


      endif

#endif /* NO_PERMANENT_GUARDCELLS */

      return
      end subroutine amr_bc_block
