!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_guardcell_srl(mype,pe,lb,iblock,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag,nlayers0x,
     .                                  nlayers0y,nlayers0z)

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace

      use paramesh_interfaces, only : amr_1blk_cc_cp_remote,
     &                                amr_1blk_fc_cp_remote,
     &                                amr_1blk_ec_cp_remote,
     &                                amr_1blk_nc_cp_remote,
     &                                amr_1blk_bcset

!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks required to fill guard cells on block (pe,lb), assuming that 
! exchange is only required between blocks at the same refinement level.
! The actual exchanges are performed with calls to the routines 
! amr_1blk_cc_cp_remote and amr_1blk_fc_cp_remote.
!
!
! Written :     Peter MacNeice          July 1998
! Modified:     Rick DeVore             February 2001
! Modified:     Peter MacNeice          February 2001
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

      implicit none

#include "amr_shmem.fh"

#ifdef TIMING_MPI
      include 'mpif.h'
#endif /* TIMING_MPI */

      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayers0x,nlayers0y,nlayers0z

!------------------------------------
! local arrays

        integer,save :: remote_blk,remote_pe,remote_type
        integer :: nguard0,nlayers0,ng0,jfl,jfu,jface
        integer :: ilays,jlays,klays,id,jd,kd,is,js,ks
        integer :: ip1,jp1,kp1,ip2,jp2,kp2,ibc,ii,jj,kk
        integer :: local_blk_type
        integer :: nblk_ind
        integer :: ip3, jp3, kp3, ip4, jp4, kp4
        integer :: ip5, jp5, kp5, ip6, jp6, kp6
        integer :: ibnd, jbnd, kbnd
#ifdef TIMING_MPI
        double precision :: time1
        double precision :: time2
#endif /* TIMING_MPI */

!------------------------------------


#ifdef TIMING_MPI
        time1 = mpi_wtime()
#endif /* TIMING_MPI */

! nblk_ind is not needed for any of the cc_XX type calls in this routine
      nblk_ind = -1

      local_blk_type = surrblks(3,2,2,2)

      if(iopt.eq.1) then
           nguard0 = nguard
           nlayers0 = nguard
      elseif(iopt.ge.2) then 
           nguard0 = nguard_work
!           nlayers0 = nlayers
           nlayers0 = max(nlayers0x,nlayers0y,nlayers0z)   ! looks like nlayers is obsolete
      endif

      ng0 = nguard0*npgs


!
! error trapping!
      if(nlayers0.gt.nguard0) then
         write(*,*) ' nguard = ',nguard
         write(*,*) ' nguard_work = ',nguard_work
         write(*,*) ' nlayers0 ',nlayers0,' nguard0 ',nguard0,
     .              ' iopt ',iopt
         write(*,*) 'amr_1blk_guardcell_srl : Too many guardcell ',
     .              'layers requested to be filled'
         call amr_abort()
      endif


      jfl = 1
      jfu = nfaces
      if(icoord.gt.0) then
        jfl = 1 + 2*(icoord-1)
        jfu = jfl + 1
      endif



!
! cycle through block faces
      do jface = jfl,jfu


! Default array index limits

! Range - source indeces are initially computed as though there
! are no permanent guardcells.
        ilays = nxb
        jlays = nyb*k2d
        klays = nzb*k3d
! Starting indeces on destination working block
        id = 1 + nguard0 
        jd = 1 + nguard0*k2d
        kd = 1 + nguard0*k3d
! Starting indeces on source block
        is = nxb
        js = nyb
        ks = nzb

        ip1 = 0
        jp1 = 0
        kp1 = 0
        ip2 = 0
        jp2 = 0
        kp2 = 0
        ip3 = 0
        jp3 = 0
        kp3 = 0
        ip4 = 0
        jp4 = 0
        kp4 = 0
        ip5 = 0
        jp5 = 0
        kp5 = 0
        ip6 = 0
        jp6 = 0
        kp6 = 0

        ibnd = 0
        jbnd = 0
        kbnd = 0
        if(jface.eq.1) ibnd = -1
        if(jface.eq.2) ibnd =  1
        if(jface.eq.3) jbnd = -1
        if(jface.eq.4) jbnd =  1
        if(jface.eq.5) kbnd = -1
        if(jface.eq.6) kbnd =  1

        if(jface.eq.1) then
          remote_blk = surrblks(1,1,2,2)
          remote_pe  = surrblks(2,1,2,2)
          remote_type = surrblks(3,1,2,2)
          id   = 1 + nguard0 - nlayers0x
          is   = 1 + nxb - nlayers0x - gc_off_x
          js   = 1
          ks   = 1
          ip3 = 1
          jp2 = 1
          kp2 = 1
          ilays = nlayers0x
          if(lrestrict_in_progress) ip5 = 1
        elseif(jface.eq.2) then
          remote_blk = surrblks(1,3,2,2)
          remote_pe  = surrblks(2,3,2,2)
          remote_type  = surrblks(3,3,2,2)
          id   = 1 + nxb + nguard0 
          is   = 1 + gc_off_x
          js   = 1
          ks   = 1
          ip3 = 1
          ip4 = 1
          ip1 = 1
          jp2 = 1
          kp2 = 1
          ilays = nlayers0x
          if(lrestrict_in_progress) ip5 = 1
          if(lrestrict_in_progress) ip6 = 1
        elseif(jface.eq.3) then
          remote_blk = surrblks(1,2,1,2)
          remote_pe  = surrblks(2,2,1,2)
          remote_type  = surrblks(3,2,1,2)
          jd   = 1 + nguard0 - nlayers0y
          js   = 1 + nyb - nlayers0y - gc_off_y
          is   = 1
          ks   = 1
          jp3 = 1
          ip2 = 1
          kp2 = 1
          jlays = nlayers0y
          if(lrestrict_in_progress) jp5 = 1
        elseif(jface.eq.4) then
          remote_blk = surrblks(1,2,3,2)
          remote_pe  = surrblks(2,2,3,2)
          remote_type  = surrblks(3,2,3,2)
          jd   = 1 + nyb + nguard0 
          js   = 1 + gc_off_y
          is   = 1
          ks   = 1
          jp3 = 1
          jp4 = 1
          jp1 = 1
          ip2 = 1
          kp2 = 1
          jlays = nlayers0y
          if(lrestrict_in_progress) jp5 = 1
          if(lrestrict_in_progress) jp6 = 1
        elseif(jface.eq.5) then
          remote_blk = surrblks(1,2,2,1)
          remote_pe  = surrblks(2,2,2,1)
          remote_type  = surrblks(3,2,2,1)
          kd   = 1 + nguard0 - nlayers0z
          ks   = 1 + nzb - nlayers0z - gc_off_z
          is   = 1
          js   = 1
          kp3 = 1
          ip2 = 1
          jp2 = 1
          klays = nlayers0z
          if(lrestrict_in_progress) kp5 = 1
        elseif(jface.eq.6) then
          remote_blk = surrblks(1,2,2,3)
          remote_pe  = surrblks(2,2,2,3)
          remote_type  = surrblks(3,2,2,3)
          kd   = 1 + nzb + nguard0
          ks   = 1 + gc_off_z
          is   = 1
          js   = 1
          kp3 = 1
          kp4 = 1
          kp1 = 1
          ip2 = 1
          jp2 = 1
          klays = nlayers0z
          if(lrestrict_in_progress) kp5 = 1
          if(lrestrict_in_progress) kp6 = 1

        endif


! Offset source indeces by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.0) then

! do not execute this section if local block is a leaf block and
! neighbor block is a parent. This will be taken care of by the
! routine amr_1blk_guardcell_f_to_c.
#ifndef ADVANCE_ALL_LEVELS
!         if(.not.(local_blk_type.eq.1.and.remote_type.eq.2)) then
#endif

#ifdef TIMING_MPIX
        time2 = mpi_wtime()
#endif /* TIMING_MPI */



          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,remote_pe,remote_blk,iblock,iopt,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,
     .                 nblk_ind)
#ifdef TIMING_MPIX
        timer_amr_1blk_cc_cp_remote(0) = 
     .     timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
#endif /* TIMING_MPI */

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 ip2,jp2,kp2,jface,
     .                 nblk_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 ip3,jp3,kp3,ip3,jp3,kp3,jface,
     .                 nblk_ind)
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 nblk_ind)

#ifndef ADVANCE_ALL_LEVELS
!         endif                  ! end of local_blk_type if test
#endif

        elseif(remote_blk.le.-20) then
          ibc = remote_blk
             

          call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,
     .                         ibnd,jbnd,kbnd,surrblks)



        endif

      enddo


!------------------------------------

      if(ldiag) then

!------------------------------------
      if(ndim.ge.2) then
      if(icoord.ne.3) then
!
! Now fill from edges along the z axis.

! Loop over the 4 corners
      do jj = 1,3,2
      do ii = 1,3,2

!
!
! Reset default index ranges
          klays = nzb*k3d
          kd = 1 + nguard0*k3d
          ks = 1

          ip1 = 0
          jp1 = 0
          kp1 = 0

          ip2 = 1
          jp2 = 1
          kp2 = 1

          ip3 = 1
          jp3 = 1
          kp3 = 0
          kp4 = 0

          ip5 = 1
          jp5 = 1
          kp5 = 0
          ip6 = 0
          jp6 = 0
          kp6 = 0
          if(lrestrict_in_progress) then
            if(ii.eq.3 )ip6 = 1
            if(jj.eq.3 )jp6 = 1
          endif
    
          remote_blk = surrblks(1,ii,jj,2)
          remote_pe  = surrblks(2,ii,jj,2)
          remote_type  = surrblks(3,ii,jj,2)


          ilays = nlayers0x
          jlays = nlayers0y

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x) + (ii-2)*gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard0
          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)*gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard0

          ip4 = mod(ii/2,2)
          jp4 = mod(jj/2,2)

! Offset source indeces by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        if(ii.eq.1) ibnd = -1
        if(ii.eq.3) ibnd =  1
        if(jj.eq.1) jbnd = -1
        if(jj.eq.3) jbnd =  1
!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.0) then
#ifndef ADVANCE_ALL_LEVELS
!        if(.not.(local_blk_type.eq.1.and.remote_type.eq.2)) then
#endif
#ifdef TIMING_MPIX
        time2 = mpi_wtime()
#endif /* TIMING_MPI */
          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,remote_pe,remote_blk,iblock,iopt,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,
     .                 nblk_ind)
#ifdef TIMING_MPIX
        timer_amr_1blk_cc_cp_remote(0) = 
     .     timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
#endif /* TIMING_MPI */

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblk_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblk_ind)
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,0,0,0,
     .                 nblk_ind)

#ifndef ADVANCE_ALL_LEVELS
!         endif                  ! end of local_blk_type if test
#endif
        elseif(remote_blk.le.-20) then
          ibc = remote_blk
          call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,
     .                         ibnd,jbnd,kbnd,surrblks)
        endif


      enddo
      enddo

      endif                                 ! end of if test icoord=3
      endif                                 ! end of if test for 2D corners


      if(ndim.eq.3) then
      if(icoord.ne.2) then
!
! Now fill from edges along the y axis.

! Loop over the 4 corners
      do kk = 1,3,2
      do ii = 1,3,2
    
!
! Reset default index ranges
          jlays = nyb*k2d
          jd = 1 + nguard0*k2d
          js = 1

          ip1 = 0
          jp1 = 0
          kp1 = 0

          ip2 = 1
          jp2 = 1
          kp2 = 1

          ip3 = 1
          jp3 = 0
          kp3 = 1
          jp4 = 0
    
          remote_blk = surrblks(1,ii,2,kk)
          remote_pe  = surrblks(2,ii,2,kk)
          remote_type  = surrblks(3,ii,2,kk)
          ilays = nlayers0x
          klays = nlayers0z

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x)  + (ii-2)*gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard0
          ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z)  + (kk-2)*gc_off_z
          kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard0

          ip4 = mod(ii/2,2)
          kp4 = mod(kk/2,2)

! Offset source indeces by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        if(ii.eq.1) ibnd = -1
        if(ii.eq.3) ibnd =  1
        if(kk.eq.1) kbnd = -1
        if(kk.eq.3) kbnd =  1

!
!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.0) then
#ifndef ADVANCE_ALL_LEVELS
!         if(.not.(local_blk_type.eq.1.and.remote_type.eq.2)) then
#endif
#ifdef TIMING_MPIX
        time2 = mpi_wtime()
#endif /* TIMING_MPI */
          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,remote_pe,remote_blk,iblock,iopt,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,
     .                 nblk_ind)
#ifdef TIMING_MPIX
        timer_amr_1blk_cc_cp_remote(0) = 
     .     timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
#endif /* TIMING_MPI */

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblk_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblk_ind)

          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,0,0,0,
     .                 nblk_ind)

#ifndef ADVANCE_ALL_LEVELS
!       endif                  ! end of local_blk_type if test
#endif
        elseif(remote_blk.le.-20) then
          ibc = remote_blk
          call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,
     .                         ibnd,jbnd,kbnd,surrblks)
        endif


      enddo
      enddo

      endif                                 ! end of if test icoord=2

!
! Now fill from edges along the x axis.
      if(icoord.ne.1) then

! Loop over the 4 corners
      do kk = 1,3,2
      do jj = 1,3,2
    
!
! Reset default index ranges
          ilays = nxb
          id = 1 + nguard0
          is = 1

          ip1 = 0
          jp1 = 0
          kp1 = 0

          jp2 = 1
          kp2 = 1
          ip2 = 1
    
          jp3 = 1
          kp3 = 1
          ip3 = 0
          ip4 = 0
    

          remote_blk = surrblks(1,2,jj,kk)
          remote_pe  = surrblks(2,2,jj,kk)
          remote_type  = surrblks(3,2,jj,kk)
          jlays = nlayers0y
          klays = nlayers0z

          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)*gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard0
          ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z) + (kk-2)*gc_off_z
          kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard0

          jp4 = mod(jj/2,2)
          kp4 = mod(kk/2,2)

! Offset source indeces by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d


        ibnd = 0
        jbnd = 0
        kbnd = 0
        if(jj.eq.1) jbnd = -1
        if(jj.eq.3) jbnd =  1
        if(kk.eq.1) kbnd = -1
        if(kk.eq.3) kbnd =  1

!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.0) then
#ifndef ADVANCE_ALL_LEVELS
!         if(.not.(local_blk_type.eq.1.and.remote_type.eq.2)) then
#endif
#ifdef TIMING_MPIX
        time2 = mpi_wtime()
#endif /* TIMING_MPI */
          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,remote_pe,remote_blk,iblock,iopt,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,
     .                 nblk_ind)
#ifdef TIMING_MPIX
        timer_amr_1blk_cc_cp_remote(0) = 
     .     timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
#endif /* TIMING_MPI */

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblk_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblk_ind)

          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,0,0,0,
     .                 nblk_ind)

#ifndef ADVANCE_ALL_LEVELS
!         endif                  ! end of local_blk_type if test
#endif
        elseif(remote_blk.le.-20) then
          ibc = remote_blk
          call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,
     .                         ibnd,jbnd,kbnd,surrblks)
        endif


      enddo
      enddo

      endif                                 ! end of if test icoord=1

      endif                                 ! end of if test for x and y edges

!------------------------------------

! Finally fill corners in 3D.

!------------------------------------
      if(ndim.eq.3) then
!
! Loop over the 8 corners
      do kk = 1,3,2
      do jj = 1,3,2
      do ii = 1,3,2
    

          remote_blk = surrblks(1,ii,jj,kk)
          remote_pe  = surrblks(2,ii,jj,kk)
          remote_type  = surrblks(3,ii,jj,kk)

          ilays = nlayers0x
          jlays = nlayers0y
          klays = nlayers0z

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x) + (ii-2)* gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard0
          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)* gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard0
          ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z) + (kk-2)* gc_off_z
          kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard0

          ip1 = 0
          jp1 = 0
          kp1 = 0

          ip4 = mod(ii/2,2)
          jp4 = mod(jj/2,2)
          kp4 = mod(kk/2,2)

          ip2 = 1
          jp2 = 1
          kp2 = 1

          ip3 = 1
          jp3 = 1
          kp3 = 1
    
! Offset source indeces by the no. of permanent guardcells
        is = is + ng0
        js = js + ng0*k2d
        ks = ks + ng0*k3d

        ibnd = 0
        jbnd = 0
        kbnd = 0
        if(ii.eq.1) ibnd = -1
        if(ii.eq.3) ibnd =  1
        if(jj.eq.1) jbnd = -1
        if(jj.eq.3) jbnd =  1
        if(kk.eq.1) kbnd = -1
        if(kk.eq.3) kbnd =  1

!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        if(remote_blk.gt.0) then
#ifndef ADVANCE_ALL_LEVELS
!         if(.not.(local_blk_type.eq.1.and.remote_type.eq.2)) then
#endif
#ifdef TIMING_MPIX
        time2 = mpi_wtime()
#endif /* TIMING_MPI */
          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,remote_pe,remote_blk,iblock,iopt,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,
     .                 nblk_ind)
#ifdef TIMING_MPIX
        timer_amr_1blk_cc_cp_remote(0) = 
     .     timer_amr_1blk_cc_cp_remote(0)+ mpi_wtime() - time2
#endif /* TIMING_MPI */

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblk_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblk_ind)
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,remote_pe,remote_blk,iblock,
     .                 id,jd,kd,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,0,0,0,
     .                 nblk_ind)
#ifndef ADVANCE_ALL_LEVELS
!         endif                  ! end of local_blk_type if test
#endif
        elseif(remote_blk.le.-20) then

          ibc = remote_blk
          call amr_1blk_bcset( mype,ibc,lb,pe,iblock,iopt,
     .                         ibnd,jbnd,kbnd,surrblks)
        endif


      enddo
      enddo
      enddo

      endif                                 ! end of loop over 3D corners

!------------------------------------

      endif                                 ! end of ldiag if test

!------------------------------------
#ifdef TIMING_MPI
       timer_amr_1blk_guardcell_srl =  timer_amr_1blk_guardcell_srl
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */


      return
      end subroutine amr_1blk_guardcell_srl
