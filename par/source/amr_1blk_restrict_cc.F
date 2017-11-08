!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_restrict_cc(mype,iopt,lb)

???? Not yet written to accomodate nxb odd ???
The difficulty is the overlap between restricted data computed
by siblings. How do we provide storage for this.




!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives interior data
! only, not guard cell data, and it receives it in a pattern which
! needs reassembly after all children have reported.
!
! This routine calls a user provided routine called restrict_fun
! which defines the pattern of restriction which the user wishes to
! apply.
!
! Written :     Peter MacNeice          January 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

#include "amr_shmem.fh"


      integer nxbs,nybs,nzbs
      parameter(nxbs = (nxb+1)/2)
      parameter(nybs = (nyb+1)/2)
      parameter(nzbs = (nzb+1)/2)
      integer len_send_cc
      parameter(len_send_cc=nxbs*nybs*nzbs)

      

!------------------------------------
! local arrays
      real send(nvar,1:nxbs,1:nybs,1:nzbs)
      real sendws(nvar,1:nxbs,1:nybs,1:nzbs)
      real recv(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
      real temp(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

      integer parent_pe,parent_blk,lb
      integer jchild
      save    send,sendws

!------------------------------------

       

! identify which child this leaf block represents
       jchild = which_child(lb)

! compute the offset in the parent block appropriate for this child
       ioff = mod(jchild-1,2)*nxb/2
       joff = mod((jchild-1)/2,2)*nyb/2
       koff = mod((jchild-1)/4,2)*nzb/2

       parent_blk = parent(1,lb)
       parent_pe  = parent(2,lb)


! compute its restricted data needed by its parent.
       if(iopt.eq.1) then


       if(mod(nxb,2).eq.0) then
         call amr_restrict_unk_fun(unk1(:,:,:,:,1),temp)
       else
         call amr_restrict_unk_fun_recip(unk(:,:,:,:,1),temp)
       endif


       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1
             do ivar=1,nvar
               send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

       ioffset = (jchild-1)*nvar*len_send_cc
       call shmem_real_put(unk(1+ioffset,1,1,1,parent_blk)
     .                    ,send(1,1,1,1),nvar*len_send_cc, parent_pe)


       elseif(iopt.eq.2) then

       if(mod(nxb,2).eq.0) then
         call amr_restrict_work_fun(work1(:,:,:,1),tempw1)
       else
         call amr_restrict_work_fun_recip(work1(:,:,:,1),tempw1)
       endif

       do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
         kk = (k-nguard_work*k3d)/2+1
         do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
           jj = (j-nguard_work*k2d)/2+1
           do i=1+nguard_work,nxb+nguard_work,2
             ii = (i-nguard_work)/2+1
             sendws(ii,jj,kk) = tempw1(i,j,k)
           enddo
         enddo
       enddo

       ioffset = (jchild-1)*len_send_cc
       call shmem_real_put(work(1+ioffset,1,1,1,parent_blk)
     .                    ,sendws(1,1,1,1),len_send_cc, parent_pe)


      endif                        ! end of iopt if test

      return
      end
