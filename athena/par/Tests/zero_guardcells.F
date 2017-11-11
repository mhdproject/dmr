!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine zero_guardcells(ioptw)



! This routine sets all guard cell data values to zero.
!
! Written:      Peter MacNeice August 1997
!
!


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace


      integer ioptw

      integer ka,kb,kc,kd
      integer kbw,kcw,kdw
#ifndef LIBRARY
      parameter(ka=1)
      parameter(kb=1+(nguard-1)*k3d)
      parameter(kc=1+(nzb+nguard)*k3d)
      parameter(kd=nzb+2*nguard*k3d)
      parameter(kbw=1+(nguard_work-1)*k3d)
      parameter(kcw=1+(nzb+nguard_work)*k3d)
      parameter(kdw=nzb+2*nguard_work*k3d)
#else
      ka=1
      kb=1+(nguard-1)*k3d
      kc=1+(nzb+nguard)*k3d
      kd=nzb+2*nguard*k3d
      kbw=1+(nguard_work-1)*k3d
      kcw=1+(nzb+nguard_work)*k3d
      kdw=nzb+2*nguard_work*k3d
#endif

#ifndef NO_PERMANENT_GUARDCELLS

      unk(:,1:nguard,:,:,:) = 0.
      unk(:,1+nxb+nguard:nxb+2*nguard,:,:,:) = 0.
      unk(:,:,1:nguard,:,:) = 0.
      unk(:,:,1+nyb+nguard:nyb+2*nguard,:,:) = 0.
      if(ndim.eq.3) then
      unk(:,:,:,ka:kb,:) = 0.
      unk(:,:,:,kc:kd,:) = 0.
      endif

!!!      write(*,*) 'unk zeroed'

      unk_n(:,1:nguard,:,:,:) = 0.
      unk_n(:,2+nxb+nguard:nxb+2*nguard+1,:,:,:) = 0.
      unk_n(:,:,1:nguard,:,:) = 0.
      unk_n(:,:,2+nyb+nguard:nyb+2*nguard+1,:,:) = 0.
      if(ndim.eq.3) then
      unk_n(:,:,:,ka:kb,:) = 0.
      unk_n(:,:,:,kc+k3d:kd+k3d,:) = 0.
      endif

!!!      write(*,*) 'unk_n zeroed'

      unk_e_x(:,1:nguard,:,:,:) = 0.
      unk_e_x(:,1+nxb+nguard:nxb+2*nguard,:,:,:) = 0.
      unk_e_x(:,:,1:nguard,:,:) = 0.
      unk_e_x(:,:,2+nyb+nguard:nyb+2*nguard+1,:,:) = 0.
      if(ndim.eq.3) then
      unk_e_x(:,:,:,ka:kb,:) = 0.
      unk_e_x(:,:,:,kc+k3d:kd+k3d,:) = 0.
      endif

      unk_e_y(:,1:nguard,:,:,:) = 0.
      unk_e_y(:,2+nxb+nguard:nxb+2*nguard+1,:,:,:) = 0.
      unk_e_y(:,:,1:nguard,:,:) = 0.
      unk_e_y(:,:,1+nyb+nguard:nyb+2*nguard,:,:) = 0.
      if(ndim.eq.3) then
      unk_e_y(:,:,:,ka:kb,:) = 0.
      unk_e_y(:,:,:,kc+k3d:kd+k3d,:) = 0.
      endif

      unk_e_z(:,1:nguard,:,:,:) = 0.
      unk_e_z(:,2+nxb+nguard:nxb+2*nguard+1,:,:,:) = 0.
      unk_e_z(:,:,1:nguard,:,:) = 0.
      unk_e_z(:,:,2+nyb+nguard:nyb+2*nguard+1,:,:) = 0.
      if(ndim.eq.3) then
      unk_e_z(:,:,:,ka:kb,:) = 0.
      unk_e_z(:,:,:,kc:kd,:) = 0.
      endif

!!!      write(*,*) 'unk_e zeroed'

      work(1:nguard_work,:,:,:,ioptw-1) = 0.
      work(1+nxb+nguard_work:nxb+2*nguard_work,
     .                                     :,:,:,ioptw-1) = 0.
      work(:,1:nguard_work,:,:,ioptw-1) = 0.
      work(:,1+nyb+nguard_work:nyb+2*nguard_work,
     .                                       :,:,ioptw-1) = 0.
      if(ndim.eq.3) then
      work(:,:,1:kbw,:,ioptw-1) = 0.
      work(:,:,kcw:kdw,:,ioptw-1) = 0.
      endif
!!!      write(*,*) 'work zeroed'

      facevarx(:,1:nguard,:,:,:) = 0.
      facevarx(:,2+nxb+nguard:1+nxb+2*nguard,:,:,:) = 0.
      facevarx(:,:,1:nguard,:,:) = 0.
      facevarx(:,:,1+nyb+nguard:nyb+2*nguard,:,:) = 0.
      if(ndim.eq.3) then
      facevarx(:,:,:,ka:kb,:) = 0.
      facevarx(:,:,:,kc:kd,:) = 0.
      endif

      facevary(:,1:nguard,:,:,:) = 0.
      facevary(:,1+nxb+nguard:nxb+2*nguard,:,:,:) = 0.
      facevary(:,:,1:nguard,:,:) = 0.
      facevary(:,:,2+nyb+nguard:1+nyb+2*nguard,:,:) = 0.
      if(ndim.eq.3) then
      facevary(:,:,:,ka:kb,:) = 0.
      facevary(:,:,:,kc:kd,:) = 0.
      endif

      facevarz(:,1:nguard,:,:,:) = 0.
      facevarz(:,1+nxb+nguard:nxb+2*nguard,:,:,:) = 0.
      facevarz(:,:,1:nguard,:,:) = 0.
      facevarz(:,:,1+nyb+nguard:nyb+2*nguard,:,:) = 0.
      if(ndim.eq.3) then
      facevarz(:,:,:,ka:kb,:) = 0.
      facevarz(:,:,:,1+(1+nzb+nguard)*k3d:1+(nzb+2*nguard)*k3d,:)=0.
      endif
!!!      write(*,*) 'facevar zeroed'

#endif

      return
      end
