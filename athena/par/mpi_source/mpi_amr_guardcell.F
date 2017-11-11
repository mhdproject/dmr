!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_guardcell(mype,iopt,nlayers,
     .                         nlayersx,nlayersy,nlayersz)



!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information.
! It cycles through the blocks, applying the routine amr_1blk_guardcell
! to each block in turn.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------
!
! Arguments:
!      mype         integer          local processor number
!      iopt         integer          indicates which data arrays are to
!                                    have their guardcells filled 
!                                      iopt = 1    unk, facevar[x][y][z]
!                                      iopt >= 2   work
!      nlayers      integer          no. of guardcell layers at each block
!                                    boundary
!      nprocs       integer          no. of processors
!
!------------------------------------

      use paramesh_dimensions
      use physicaldata
      use workspace
      use tree

      use paramesh_interfaces, only : amr_1blk_guardcell_reset,
     &                                amr_restrict,
     &                                amr_1blk_guardcell

      use paramesh_mpi_interfaces, only : mpi_amr_comm_setup

      implicit none

      integer, intent(in) :: mype,iopt,nlayers
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

      logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag,l_force_consist
      integer :: lb,icoord
      integer :: id,jd,kd
      integer :: ilays,jlays,klays
      integer :: nlayers0x, nlayers0y, nlayers0z, nguard0
      integer :: i,j,k,ivar
      integer :: ip1,ip2,jp1,jp2,kp1,kp2
      integer :: ilp,iup,jlp,jup,klp,kup
      integer :: nprocs, tag_offset, iempty, iu, ju, ku, iopt0
      integer,external :: shmem_n_pes


#ifdef TIMINGS
#include "timer.fh"
#endif

!------------------------------------
#ifndef DIAGONALS
      write(*,*) 'amr_guardcell:  diagonals off'
#endif

       if (iopt == 1) then

! set users selections of guardcell variables
         int_gcell_on_cc = gcell_on_cc

         if (.not.present(nlayersx)) then
            nlayers0x = nguard
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard
         else
            nlayers0z = nlayersz
         end if
      else
         if (.not.present(nlayersx)) then
            nlayers0x = nguard_work
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard_work
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard_work
         else
            nlayers0z = nlayersz
         end if
      end if

      if (iopt == 1) then
        nguard0 = nguard
      else
        nguard0 = nguard_work
      end if

      nprocs = shmem_n_pes()

#ifdef NO_PERMANENT_GUARDCELLS
      if(mype.eq.0) then
        write(*,*) 'amr_guardcell call ignored!'
        write(*,*) 'NO_PERMANENT_GUARDCELLS is defined'
      endif
      return
#else /* NO_PERMANENT_GUARDCELLS */

! make sure that nlayers and iopt are set consistently.
      if(iopt.eq.1.and.nlayers.ne.nguard) then
        if(mype.eq.0) then
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'Error in guardcell - iopt and nlayers'
          write(*,*) 'are not consistent. For iopt=1 you must'
          write(*,*) 'set nlayers=nguard.'
        endif
        call amr_abort
      elseif(iopt.ge.2.and.nlayers.gt.nguard_work) then
        if(mype.eq.0) then
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'Error in guardcell - iopt and nlayers'
          write(*,*) 'are not consistent. For iopt>=2 you must'
          write(*,*) 'set nlayers le nguard_work.'
        endif
        call amr_abort
      endif

! reinitialize addresses of cached parent blocks
      call amr_1blk_guardcell_reset

#ifdef TIMINGS
      call shmem_barrier_all()
      itimer1 = irtc()
#endif



      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(iopt.eq.1) then
        if(nvar.gt.0) lcc = .true.
        if(nfacevar.gt.0) lfc = .true.
        if(nvaredge.gt.0) lec = .true.
        if(nvarcorn.gt.0) lnc = .true.
      elseif(iopt.ge.2) then
        lcc = .true.
      endif

!
! Restrict solution to parent blocks
#ifndef ADVANCE_ALL_LEVELS

      iempty = 0
      call amr_restrict(mype,iopt,iempty,.true.)


      call amr_1blk_guardcell_reset

#endif

      l_force_consist = .false.
#ifdef FORCE_CONSISTENCY_AT_SRL_INTERFACES
      l_force_consist = .true.
      if(lfc) then
        do lb = 1,lnblocks
          gt_facevarx(:,1,:,:,lb) = facevarx(:,1+nguard,:,:,lb)
          gt_facevarx(:,2,:,:,lb) = facevarx(:,nxb+1+nguard,:,:,lb)
          if(ndim.ge.2) then
          gt_facevary(:,:,1,:,lb) = facevary(:,:,1+nguard*k2d,:,lb)
          gt_facevary(:,:,1+k2d,:,lb) =  
     .                        facevary(:,:,nyb+(1+nguard)*k2d,:,lb)
          endif
          if(ndim.eq.3) then
          gt_facevarz(:,:,:,1,lb) = facevarz(:,:,:,1+nguard*k3d,lb)
          gt_facevarz(:,:,:,1+k3d,lb) =  
     .                        facevarz(:,:,:,nzb+(1+nguard)*k3d,lb)
          endif
        enddo
      endif
#endif  /* FORCE_CONSISTENCY_AT_SRL_INTERFACES */


      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,ledge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset,
     .                        nlayersx,nlayersy,nlayersz)


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb) == 1 .or. nodetype(lb) == 2) then
#endif

! Copy this blocks data into the working block, and fill its guardcells
        ldiag = .false.
#ifdef DIAGONALS
        ldiag = .true.
#endif


        l_srl_only = .false.                     ! fill srl and coarse
        icoord = 0                               ! fill in all coord directions
        call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,
     .                          lcc,lfc,lec,lnc,
     .                          l_srl_only,icoord,ldiag,
     .                          nlayersx,nlayersy,nlayersz)

        do k = 1,1+2*k3d
         klp = 0
         kup = 0
         if (k == 1) then
           klays = nlayers0z*k3d
           kd = nguard0*k3d+1 - klays
           kp1 = 0
           kp2 = 0
           if(l_force_consist) kup = k3d
         else if (k == 2) then
           klays = nzb*k3d
           kd = nguard0*k3d+1
           kp1 = 0
           kp2 = k3d
         else if (k == 3) then
           klays = nlayers0z*k3d
           kd = (nguard0+nzb)*k3d + 1
           kp1 = k3d
           kp2 = k3d
           if(l_force_consist) klp = -k3d
         end if
         ku = kd + klays - k3d
        do j = 1,1+2*k2d
         jlp = 0
         jup = 0
         if (j == 1) then
           jlays = nlayers0y*k2d
           jd = nguard0*k2d+1 - jlays
           jp1 = 0
           jp2 = 0
           if(l_force_consist) jup = k2d
         else if (j == 2) then
           jlays = nyb*k2d
           jd = nguard0*k2d+1
           jp1 = 0
           jp2 = k2d
         else if (j == 3) then
           jlays = nlayers0y*k2d
           jd = (nguard0+nyb)*k2d + 1
           jp1 = k2d
           jp2 = k2d
           if(l_force_consist) jlp = -k2d
         end if
         ju = jd + jlays - k2d
       do i = 1,3
         ilp = 0
         iup = 0
         if (i == 1) then
           ilays = nlayers0x
           id = nguard0+1 - ilays
           ip1 = 0
           ip2 = 0
           if(l_force_consist) iup = 1
           
         else if (i == 2) then
           ilays = nxb
           id = nguard0+1
           ip1 = 0
           ip2 = 1
         else if (i == 3) then
           ilays = nlayers0x
           id = nguard0+nxb + 1
           ip1 = 1
           ip2 = 1
           if(l_force_consist) ilp = -1
         else if (i == 3) then
         end if
         iu = id + ilays - 1

         if (i .eq. 2 .and. j .eq. 1+k2d .and. k .eq. 1+k3d) then

         else

        if (lcc) then
         if (iopt == 1) then
           do ivar=1,nvar
             if(int_gcell_on_cc(ivar)) then
               unk(ivar,id:iu,jd:ju,kd:ku,lb) = 
     .           unk1(ivar,id:iu,jd:ju,kd:ku,1)
             endif
           enddo
         else
          iopt0 = iopt-1
          work(id:iu,jd:ju,kd:ku,lb,iopt0) = 
     .           work1(id:iu,jd:ju,kd:ku,1)
         end if
        end if

        if (lfc) then
          facevarx( 1:nfacevar,id+ip1+ilp:iu+ip2+iup,
     .                         jd:ju,kd:ku,lb) =
     .    facevarx1(1:nfacevar,id+ip1+ilp:iu+ip2+iup,
     .                         jd:ju,kd:ku,1)

          if (ndim > 1) then
          facevary( 1:nfacevar,id:iu,jd+jp1+jlp:ju+jp2+jup,
     .                               kd:ku,lb) =
     .    facevary1(1:nfacevar,id:iu,jd+jp1+jlp:ju+jp2+jup,
     .                               kd:ku,1)
          end if

          if (ndim == 3) then
          facevarz( 1:nfacevar,id:iu,jd:ju,kd+kp1+klp:ku+kp2+kup,lb) =
     .    facevarz1(1:nfacevar,id:iu,jd:ju,kd+kp1+klp:ku+kp2+kup,1)
          end if

        end if

        if (lec) then
           if (ndim > 1) then

           unk_e_x( 1:nvaredge,id:iu,jd+jp1:ju+jp2,kd+kp1:ku+kp2,lb) =
     .     unk_e_x1(1:nvaredge,id:iu,jd+jp1:ju+jp2,kd+kp1:ku+kp2,1)

           unk_e_y( 1:nvaredge,id+ip1:iu+ip2,jd:ju,kd+kp1:ku+kp2,lb) =
     .     unk_e_y1(1:nvaredge,id+ip1:iu+ip2,jd:ju,kd+kp1:ku+kp2,1)

           if (ndim == 3) then
           unk_e_z( 1:nvaredge,id+ip1:iu+ip2,jd+jp1:ju+jp2,kd:ku,lb) =
     .     unk_e_z1(1:nvaredge,id+ip1:iu+ip2,jd+jp1:ju+jp2,kd:ku,1)
           end if

           end if
        end if

        if (lnc) then
           unk_n( 1:nvarcorn,
     .            id+ip1:iu+ip2,
     .            jd+jp1:ju+jp2,
     .            kd+kp1:ku+kp2,lb) = 
     .     unk_n1(1:nvarcorn,
     .            id+ip1:iu+ip2,
     .            jd+jp1:ju+jp2,
     .            kd+kp1:ku+kp2,1)
        end if

        end if ! if (i .ne. 2   

      enddo
      enddo
      enddo

 20   continue

#ifndef ADVANCE_ALL_LEVELS
      end if ! if (nodetype == 1
#endif

      enddo                                 ! end of loop over blocks
      endif



! reinitialize addresses of cached parent blocks
      call amr_1blk_guardcell_reset


! reset selections of guardcell variables to default
      int_gcell_on_cc(:) = .true.

#ifdef TIMINGS
      itimer2 = irtc()
      call shmem_barrier_all()
      irtc_guard = itimer2-itimer1+irtc_guard
#endif

#endif /* NO_PERMANENT_GUARDCELLS */


      return
      end subroutine amr_guardcell


