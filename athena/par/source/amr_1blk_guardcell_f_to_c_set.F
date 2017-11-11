!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_guardcell_f_to_c_set(
     .                                  mype,pe,lb,iblock,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag,
     .                                  nlayersx,nlayersy,nlayersz)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : amr_1blk_cc_cp_remote,
     &                                amr_1blk_fc_cp_remote,
     &                                amr_1blk_ec_cp_remote,
     &                                amr_1blk_nc_cp_remote

!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks required to fill guard cells on block (pe,lb), assuming that 
! exchange is only required to the current block, from neighbor blocks
! at a higher refinement level.
! The actual exchanges are performed with calls to the routines 
! amr_1blk_cc_cp_remote and amr_1blk_fc_cp_remote.
!
!
! Written :     Peter MacNeice          January 2003
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


      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx,nlayersy,nlayersz

!------------------------------------
! local arrays

        integer,save :: fremote_blk(2,2),fremote_pe(2,2)
        integer,save :: remote_blk,remote_pe
        integer,save :: remote_type
        integer,save :: cchild(2,mchild)

        integer :: nguard0,nguard02,nlayers0,ng0,jfl,jfu,jface
        integer :: nlayers0x,nlayers0y,nlayers0z
        integer :: ilays,jlays,klays,id,jd,kd,is,js,ks
        integer :: ip1,jp1,kp1,ip2,jp2,kp2,ii,jj,kk,nn
        integer :: jblk
        integer :: icoord_loc
        integer :: ip3, jp3, kp3, ip4, jp4, kp4
        integer :: ip5, jp5, kp5, ip6, jp6, kp6
        integer :: ibnd, jbnd, kbnd
        integer :: nblks_ind, idf, jdf, kdf

        logical :: lfound

!------------------------------------

      icoord_loc = 0

      if(iopt.eq.1) then
           nguard0 = nguard
           nguard02 = nguard0 * 2
           nlayers0 = nguard  * 2 
      elseif(iopt.ge.2) then 
           nguard0 = nguard_work
           nguard02 = nguard0  * 2
           nlayers0 = nlayers  * 2 
      endif

!      nlayers0x = min(nlayersx*2 ,nguard02)
!      nlayers0y = min(nlayersy*2 ,nguard02)
!      nlayers0z = min(nlayersz*2 ,nguard02)
      nlayers0x = min(nlayersx*2 + 1,nguard02)
      nlayers0y = min(nlayersy*2 + 1,nguard02)
      nlayers0z = min(nlayersz*2 + 1,nguard02)


      ng0 = nguard0*npgs


!
! error trapping!
      if(nlayers0.gt.2*nguard0) then
         write(*,*) ' nlayers0 ',nlayers0,' nguard0 ',nguard0,
     .              ' iopt ',iopt
         write(*,*) 'amr_1blk_guardcell_f_to_c : Too many guardcell ',
     .              'layers requested to be filled'
         call amr_abort()
      endif
      if( (nlayers0x.gt.nxb) .or. 
     .    ((nlayers0y.gt.nyb).and.ndim.ge.2) .or.
     .    ((nlayers0z.gt.nzb).and.ndim.eq.3) ) then
         write(*,*) ' nlayers0x,y,z, ',nlayers0x,nlayers0y,nlayers0z,
     .              ' iopt ',iopt
         write(*,*) 'amr_1blk_guardcell_f_to_c : Too many ',
     .              'layers requested to be filled'
         call amr_abort()
      endif


      jfl = 1
      jfu = nfaces
      if(icoord_loc.gt.0) then
        jfl = 1 + 2*(icoord_loc-1)
        jfu = jfl + 1
      endif


!
! cycle through block faces
      do jface = jfl,jfu

        if(jface.eq.1) then
          remote_type = surrblks(3,1,2,2)
        elseif(jface.eq.2) then
          remote_type = surrblks(3,3,2,2)
        elseif(jface.eq.3) then
          remote_type = surrblks(3,2,1,2)
        elseif(jface.eq.4) then
          remote_type = surrblks(3,2,3,2)
        elseif(jface.eq.5) then
          remote_type = surrblks(3,2,2,1)
        elseif(jface.eq.6) then
          remote_type = surrblks(3,2,2,3)
        endif

! Default array index limits

! Range - source indeces are initially computed as though there
! are no permanent guardcells.
        ilays = nxb
        jlays = nyb*k2d
        klays = nzb*k3d 
! Starting indeces on destination working block
        id = 1 + nguard02
        jd = 1 + nguard02*k2d
        kd = 1 + nguard02*k3d
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
          id   = 1 + nguard02 - nlayers0x
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
          remote_type = surrblks(3,3,2,2)
          id   = 1 + nxb + nguard02
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
          remote_type = surrblks(3,2,1,2)
          jd   = 1 + nguard02 - nlayers0y
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
          remote_type = surrblks(3,2,3,2)
          jd   = 1 + nyb + nguard02
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
          remote_type = surrblks(3,2,2,1)
          kd   = 1 + nguard02 - nlayers0z
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
          remote_type = surrblks(3,2,2,3)
          kd   = 1 + nzb + nguard02
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

        if(remote_blk.gt.0.and.remote_type.eq.2) then

!
! Locate the fine blocks neighboring this face
#ifdef MPI_USED

          if(remote_pe.ne.mype) then

          lfound = .false.
          jblk = strt_buffer-1
          do while(jblk.lt.last_buffer.and.(.not.lfound))
           jblk = jblk+1
           if(remote_blk.eq.laddress(1,jblk).and.
     .        remote_pe.eq.laddress(2,jblk)) then
             lfound = .true.
             cchild(1:2,1:mchild) = child(1:2,1:mchild,jblk)
           endif
          enddo

          if(.not.lfound) then
            write(*,*) 'Error in amr_1blk_guardcell_f_to_c '
            write(*,*) 'Error: looking for blk ',
     .                  remote_blk,remote_pe
            call amr_abort()
          endif

          else
           cchild(1:2,1:mchild) = child(1:2,1:mchild,remote_blk)
          endif

#else /* MPI_USED */
          call shmem_integer_get(cchild(1,1),child(1,1,remote_blk),
     .                           2*mchild,remote_pe)
#endif /* MPI_USED */


          if(jface.eq.1) then
            fremote_blk(1,1)  = cchild(1,2)
            fremote_blk(2,1)  = cchild(1,4)
            fremote_pe(1,1)   = cchild(2,2)
            fremote_pe(2,1)   = cchild(2,4)
            if(ndim.eq.3) then
              fremote_blk(1,2)  = cchild(1,6)
              fremote_blk(2,2)  = cchild(1,8)
              fremote_pe(1,2)   = cchild(2,6)
              fremote_pe(2,2)   = cchild(2,8)
            endif
            nblks_ind = 13
          elseif(jface.eq.2) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_blk(2,1)  = cchild(1,3)
            fremote_pe(1,1)   = cchild(2,1)
            fremote_pe(2,1)   = cchild(2,3)
            if(ndim.eq.3) then
              fremote_blk(1,2)  = cchild(1,5)
              fremote_blk(2,2)  = cchild(1,7)
              fremote_pe(1,2)   = cchild(2,5)
              fremote_pe(2,2)   = cchild(2,7)
            endif
            nblks_ind = 15
          elseif(jface.eq.3) then
            fremote_blk(1,1)  = cchild(1,3)
            fremote_blk(2,1)  = cchild(1,4)
            fremote_pe(1,1)   = cchild(2,3)
            fremote_pe(2,1)   = cchild(2,4)
            if(ndim.eq.3) then
              fremote_blk(1,2)  = cchild(1,7)
              fremote_blk(2,2)  = cchild(1,8)
              fremote_pe(1,2)   = cchild(2,7)
              fremote_pe(2,2)   = cchild(2,8)
            endif
            nblks_ind = 11
          elseif(jface.eq.4) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_blk(2,1)  = cchild(1,2)
            fremote_pe(1,1)   = cchild(2,1)
            fremote_pe(2,1)   = cchild(2,2)
            if(ndim.eq.3) then
              fremote_blk(1,2)  = cchild(1,5)
              fremote_blk(2,2)  = cchild(1,6)
              fremote_pe(1,2)   = cchild(2,5)
              fremote_pe(2,2)   = cchild(2,6)
            endif
            nblks_ind = 17
          elseif(jface.eq.5) then
            fremote_blk(1,1)  = cchild(1,5)
            fremote_blk(2,1)  = cchild(1,6)
            fremote_pe(1,1)   = cchild(2,5)
            fremote_pe(2,1)   = cchild(2,6)
            fremote_blk(1,2)  = cchild(1,7)
            fremote_blk(2,2)  = cchild(1,8)
            fremote_pe(1,2)   = cchild(2,7)
            fremote_pe(2,2)   = cchild(2,8)
            nblks_ind = 5
          elseif(jface.eq.6) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_blk(2,1)  = cchild(1,2)
            fremote_pe(1,1)   = cchild(2,1)
            fremote_pe(2,1)   = cchild(2,2)
            fremote_blk(1,2)  = cchild(1,3)
            fremote_blk(2,2)  = cchild(1,4)
            fremote_pe(1,2)   = cchild(2,3)
            fremote_pe(2,2)   = cchild(2,4)
            nblks_ind = 23
          endif

          call set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)


! is,js,ks are the same as we would use for the srl routine, except that
! the number of layers we need is doubled.




          do jj = 1,ndim-1                   ! only gt 1 for 3D
          do ii = 1,2

! Set offsets to base values for id,jd,kd for each fine neighbor block
          if(jface.eq.1.or.jface.eq.2) then
            idf = id 
            if(jface.eq.2) idf = id + nxb
            jdf = jd 
            kdf = kd 
            if(ii.eq.1.and.jj.eq.1) then
            elseif(ii.eq.2.and.jj.eq.1) then
              jdf = jd + nyb
            elseif(ii.eq.1.and.jj.eq.2) then
              kdf = kd + nzb
            elseif(ii.eq.2.and.jj.eq.2) then
              jdf = jd + nyb
              kdf = kd + nzb
            endif
          elseif(jface.eq.3.or.jface.eq.4) then
            idf = id 
            jdf = jd 
            if(jface.eq.4) jdf = jd + nyb
            kdf = kd 
            if(ii.eq.1.and.jj.eq.1) then
            elseif(ii.eq.2.and.jj.eq.1) then
              idf = id + nxb
            elseif(ii.eq.1.and.jj.eq.2) then
              kdf = kd + nzb
            elseif(ii.eq.2.and.jj.eq.2) then
              idf = id + nxb
              kdf = kd + nzb
            endif
          elseif(jface.eq.5.or.jface.eq.6) then
            idf = id 
            jdf = jd 
            kdf = kd 
            if(jface.eq.6) kdf = kd + nzb
            if(ii.eq.1.and.jj.eq.1) then
            elseif(ii.eq.2.and.jj.eq.1) then
              idf = id + nxb
            elseif(ii.eq.1.and.jj.eq.2) then
              jdf = jd + nyb
            elseif(ii.eq.2.and.jj.eq.2) then
              idf = id + nxb
              jdf = jd + nyb
            endif
          endif


          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,
     .                 fremote_pe(ii,jj),fremote_blk(ii,jj),
     .                 iblock,iopt,
     .                 idf,jdf,kdf,is,js,ks,ilays,jlays,klays,
     .                 nblks_ind)

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,
     .                 fremote_pe(ii,jj),fremote_blk(ii,jj),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 ip2,jp2,kp2,jface,
     .                 nblks_ind)


          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,
     .                 fremote_pe(ii,jj),fremote_blk(ii,jj),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 ip3,jp3,kp3,jface,
     .                 nblks_ind)


          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,
     .                 fremote_pe(ii,jj),fremote_blk(ii,jj),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 nblks_ind)

           enddo
           enddo

        endif

      enddo


!------------------------------------

      if(ldiag) then

!------------------------------------
      if(ndim.ge.2) then
      if(icoord_loc.ne.3) then
!
! Now fill from edges along the z axis.

! Loop over the 4 corners
      do jj = 1,3,2
      do ii = 1,3,2

!
!
! Reset default index ranges
          klays = nzb*k3d 
          kd = 1 + nguard02*k3d
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
          remote_type = surrblks(3,ii,jj,2)


          ilays = nlayers0x
          jlays = nlayers0y

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x) + (ii-2)*gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard02
          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)*gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard02

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

        if(remote_blk.gt.0.and.remote_type.eq.2) then


!
! Locate the fine blocks neighboring this face
#ifdef MPI_USED
          if(remote_pe.ne.mype) then

          lfound = .false.
          jblk = strt_buffer-1
          do while(jblk.lt.last_buffer.and.(.not.lfound))
           jblk = jblk+1
           if(remote_blk.eq.laddress(1,jblk).and.
     .        remote_pe.eq.laddress(2,jblk)) then
             lfound = .true.
             cchild(1:2,1:mchild) = child(1:2,1:mchild,jblk)
           endif
          enddo

          if(.not.lfound) then
            write(*,*) 'Error in amr_1blk_guardcell_f_to_c '
            call amr_abort()
          endif

          else
           cchild(1:2,1:mchild) = child(1:2,1:mchild,remote_blk)
          endif

#else /* MPI_USED */
          call shmem_integer_get(cchild(1,1),child(1,1,remote_blk),
     .                           2*mchild,remote_pe)
#endif /* MPI_USED */


          if(ii.eq.1.and.jj.eq.1) then
            fremote_blk(1,1)  = cchild(1,4)
            fremote_pe(1,1)   = cchild(2,4)
            if(ndim.eq.3) then
              fremote_blk(2,1)  = cchild(1,8)
              fremote_pe(2,1)   = cchild(2,8)
            endif
            nblks_ind = 10
          elseif(ii.eq.3.and.jj.eq.1) then
            fremote_blk(1,1)  = cchild(1,3)
            fremote_pe(1,1)   = cchild(2,3)
            if(ndim.eq.3) then
              fremote_blk(2,1)  = cchild(1,7)
              fremote_pe(2,1)   = cchild(2,7)
            endif
            nblks_ind = 12
          elseif(ii.eq.1.and.jj.eq.3) then
            fremote_blk(1,1)  = cchild(1,2)
            fremote_pe(1,1)   = cchild(2,2)
            if(ndim.eq.3) then
              fremote_blk(2,1)  = cchild(1,6)
              fremote_pe(2,1)   = cchild(2,6)
            endif
            nblks_ind = 16
          elseif(ii.eq.3.and.jj.eq.3) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_pe(1,1)   = cchild(2,1)
            if(ndim.eq.3) then
              fremote_blk(2,1)  = cchild(1,5)
              fremote_pe(2,1)   = cchild(2,5)
            endif
            nblks_ind = 18
          endif

          call set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)


!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        do nn = 1,ndim-1


! Set offsets to base values for id,jd,kd for each fine neighbor block
          idf = id
          jdf = jd
          kdf = kd
          if(ii.eq.1.and.jj.eq.1) then
            if(nn.eq.1) then
            else
              kdf = kd + nzb
            endif
          elseif(ii.eq.3.and.jj.eq.1) then
            idf = id + nxb
            if(nn.eq.1) then
            else
              kdf = kd + nzb
            endif
          elseif(ii.eq.1.and.jj.eq.3) then
            jdf = jd + nyb
            if(nn.eq.1) then
            else
              kdf = kd + nzb
            endif
          elseif(ii.eq.3.and.jj.eq.3) then
            idf = id + nxb
            jdf = jd + nyb
            if(nn.eq.1) then
            else
              kdf = kd + nzb
            endif
          endif


          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,iopt,
     .                 idf,jdf,kdf,is,js,ks,ilays,jlays,klays,
     .                 nblks_ind)


          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblks_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,0,0,0,
     .                 0,0,0,0,
     .                 nblks_ind)

          if(ii.eq.3) ip1 = 1
          if(jj.eq.3) jp1 = 1
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 nblks_ind)

          if(ii.eq.3) ip1 = 0           !reset to avoid errors during nn loop
          if(jj.eq.3) jp1 = 0


        enddo

        endif

      enddo
      enddo

      endif                                 ! end of if test icoord_loc=3
      endif                                 ! end of if test for 2D corners


      if(ndim.eq.3) then
      if(icoord_loc.ne.2) then
!
! Now fill from edges along the y axis.

! Loop over the 4 corners
      do kk = 1,3,2
      do ii = 1,3,2
    
!
! Reset default index ranges
          jlays = nyb*k2d 
          jd = 1 + nguard02*k2d
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
          remote_type = surrblks(3,ii,2,kk)
          ilays = nlayers0x
          klays = nlayers0z

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x)  + (ii-2)*gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard02
          ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z)  + (kk-2)*gc_off_z
          kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard02

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



        if(remote_blk.gt.0.and.remote_type.eq.2) then
!
! Locate the fine blocks neighboring this face
#ifdef MPI_USED
          if(remote_pe.ne.mype) then

          lfound = .false.
          jblk = strt_buffer-1
          do while(jblk.lt.last_buffer.and.(.not.lfound))
           jblk = jblk+1
           if(remote_blk.eq.laddress(1,jblk).and.
     .        remote_pe.eq.laddress(2,jblk)) then
             lfound = .true.
             cchild(1:2,1:mchild) = child(1:2,1:mchild,jblk)
           endif
          enddo

          if(.not.lfound) then
            write(*,*) 'Error in amr_1blk_guardcell_f_to_c '
            call amr_abort()
          endif

          else
           cchild(1:2,1:mchild) = child(1:2,1:mchild,remote_blk)
          endif
#else /* MPI_USED */
          call shmem_integer_get(cchild(1,1),child(1,1,remote_blk),
     .                           2*mchild,remote_pe)
#endif /* MPI_USED */


          if(ii.eq.1.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,6)
            fremote_pe(1,1)   = cchild(2,6)
            fremote_blk(2,1)  = cchild(1,8)
            fremote_pe(2,1)   = cchild(2,8)
            nblks_ind = 4
          elseif(ii.eq.3.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,5)
            fremote_pe(1,1)   = cchild(2,5)
            fremote_blk(2,1)  = cchild(1,7)
            fremote_pe(2,1)   = cchild(2,7)
            nblks_ind = 6
          elseif(ii.eq.1.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,2)
            fremote_pe(1,1)   = cchild(2,2)
            fremote_blk(2,1)  = cchild(1,4)
            fremote_pe(2,1)   = cchild(2,4)
            nblks_ind = 22
          elseif(ii.eq.3.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_pe(1,1)   = cchild(2,1)
            fremote_blk(2,1)  = cchild(1,3)
            fremote_pe(2,1)   = cchild(2,3)
            nblks_ind = 24
          endif

          call set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)

!
!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        do nn = 1,ndim-1

! Set offsets to base values for id,jd,kd for each fine neighbor block
          idf = id
          jdf = jd
          kdf = kd
          if(ii.eq.1.and.kk.eq.1) then
            if(nn.eq.1) then
            else
              jdf = jd + nyb
            endif
          elseif(ii.eq.3.and.kk.eq.1) then
            idf = id + nxb
            if(nn.eq.1) then
            else
              jdf = jd + nyb
            endif
          elseif(ii.eq.1.and.kk.eq.3) then
            kdf = kd + nzb
            if(nn.eq.1) then
            else
              jdf = jd + nyb
            endif
          elseif(ii.eq.3.and.kk.eq.3) then
            idf = id + nxb
            kdf = kd + nzb
            if(nn.eq.1) then
            else
              jdf = jd + nyb
            endif
          endif


          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,iopt,
     .                 idf,jdf,kdf,is,js,ks,ilays,jlays,klays,
     .                 nblks_ind)


          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblks_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblks_ind)

          if(ii.eq.3) ip1 = 1
          if(kk.eq.3) kp1 = 1
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 nblks_ind)
          if(ii.eq.3) ip1 = 0           !reset to avoid errors during nn loop
          if(kk.eq.3) kp1 = 0


        enddo

        endif

      enddo
      enddo

      endif                                 ! end of if test icoord_loc=2

!
! Now fill from edges along the x axis.
      if(icoord_loc.ne.1) then

! Loop over the 4 corners
      do kk = 1,3,2
      do jj = 1,3,2
    
!
! Reset default index ranges
          ilays = nxb
          id = 1 + nguard0*2
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
          remote_type = surrblks(3,2,jj,kk)
          jlays = nlayers0y
          klays = nlayers0z

          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)*gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard02
          ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z) + (kk-2)*gc_off_z
          kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard02

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



        if(remote_blk.gt.0.and.remote_type.eq.2) then


!
! Locate the fine blocks neighboring this face
#ifdef MPI_USED
          if(remote_pe.ne.mype) then

          lfound = .false.
          jblk = strt_buffer-1
          do while(jblk.lt.last_buffer.and.(.not.lfound))
           jblk = jblk+1
           if(remote_blk.eq.laddress(1,jblk).and.
     .        remote_pe.eq.laddress(2,jblk)) then
             lfound = .true.
             cchild(1:2,1:mchild) = child(1:2,1:mchild,jblk)
           endif
          enddo

          if(.not.lfound) then
            write(*,*) 'Error in amr_1blk_guardcell_f_to_c '
            call amr_abort()
          endif

          else
           cchild(1:2,1:mchild) = child(1:2,1:mchild,remote_blk)
          endif

#else /* MPI_USED */
          call shmem_integer_get(cchild(1,1),child(1,1,remote_blk),
     .                           2*mchild,remote_pe)
#endif /* MPI_USED */


          if(jj.eq.1.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,7)
            fremote_pe(1,1)   = cchild(2,7)
            fremote_blk(2,1)  = cchild(1,8)
            fremote_pe(2,1)   = cchild(2,8)
            nblks_ind = 2
          elseif(jj.eq.3.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,5)
            fremote_pe(1,1)   = cchild(2,5)
            fremote_blk(2,1)  = cchild(1,6)
            fremote_pe(2,1)   = cchild(2,6)
            nblks_ind = 8
          elseif(jj.eq.1.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,3)
            fremote_pe(1,1)   = cchild(2,3)
            fremote_blk(2,1)  = cchild(1,4)
            fremote_pe(2,1)   = cchild(2,4)
            nblks_ind = 20
          elseif(jj.eq.3.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_pe(1,1)   = cchild(2,1)
            fremote_blk(2,1)  = cchild(1,2)
            fremote_pe(2,1)   = cchild(2,2)
            nblks_ind = 26
          endif

          call set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)



!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

        do nn = 1,ndim-1

! Set offsets to base values for id,jd,kd for each fine neighbor block
          idf = id
          jdf = jd
          kdf = kd
          if(jj.eq.1.and.kk.eq.1) then
            if(nn.eq.1) then
            else
              idf = id + nxb
            endif
          elseif(jj.eq.3.and.kk.eq.1) then
            jdf = jd + nyb
            if(nn.eq.1) then
            else
              idf = id + nxb
            endif
          elseif(jj.eq.1.and.kk.eq.3) then
            kdf = kd + nzb
            if(nn.eq.1) then
            else
              idf = id + nxb
            endif
          elseif(jj.eq.3.and.kk.eq.3) then
            jdf = jd + nyb
            kdf = kd + nzb
            if(nn.eq.1) then
            else
              idf = id + nxb
            endif
          endif


          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,iopt,
     .                 idf,jdf,kdf,is,js,ks,ilays,jlays,klays,
     .                 nblks_ind)

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblks_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblks_ind)

          if(jj.eq.3) jp1 = 1
          if(kk.eq.3) kp1 = 1
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,
     .                 fremote_pe(nn,1),fremote_blk(nn,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 nblks_ind)
          if(jj.eq.3) jp1 = 0           !reset to avoid errors during nn loop
          if(kk.eq.3) kp1 = 0


        enddo

        endif

      enddo
      enddo

      endif                                 ! end of if test icoord_loc=1

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
          remote_type = surrblks(3,ii,jj,kk)

          ilays = nlayers0x
          jlays = nlayers0y
          klays = nlayers0z

          is = (ii/2) + (1-ii/2)*(nxb+1-nlayers0x) + (ii-2)* gc_off_x
          id = (ii/2)*nxb + (1-ii/2)*(-nlayers0x) + 1 + nguard02
          js = (jj/2) + (1-jj/2)*(nyb+1-nlayers0y) + (jj-2)* gc_off_y
          jd = (jj/2)*nyb + (1-jj/2)*(-nlayers0y) + 1 + nguard02
          ks = (kk/2) + (1-kk/2)*(nzb+1-nlayers0z) + (kk-2)* gc_off_z
          kd = (kk/2)*nzb + (1-kk/2)*(-nlayers0z) + 1 + nguard02

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


        if(remote_blk.gt.0.and.remote_type.eq.2) then


!
! Locate the fine blocks neighboring this face
#ifdef MPI_USED
          if(remote_pe.ne.mype) then

          lfound = .false.
          jblk = strt_buffer-1
          do while(jblk.lt.last_buffer.and.(.not.lfound))
           jblk = jblk+1
           if(remote_blk.eq.laddress(1,jblk).and.
     .        remote_pe.eq.laddress(2,jblk)) then
             lfound = .true.
             cchild(1:2,1:mchild) = child(1:2,1:mchild,jblk)
           endif
          enddo

          if(.not.lfound) then
            write(*,*) 'Error in amr_1blk_guardcell_f_to_c '
            call amr_abort()
          endif

          else
           cchild(1:2,1:mchild) = child(1:2,1:mchild,remote_blk)
          endif

#else /* MPI_USED */
          call shmem_integer_get(cchild(1,1),child(1,1,remote_blk),
     .                           2*mchild,remote_pe)
#endif /* MPI_USED */

          if(ii.eq.1.and.jj.eq.1.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,8)
            fremote_pe(1,1)   = cchild(2,8)
            nblks_ind = 1
          elseif(ii.eq.3.and.jj.eq.1.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,7)
            fremote_pe(1,1)   = cchild(2,7)
            nblks_ind = 3
          elseif(ii.eq.1.and.jj.eq.3.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,6)
            fremote_pe(1,1)   = cchild(2,6)
            nblks_ind = 7
          elseif(ii.eq.3.and.jj.eq.3.and.kk.eq.1) then
            fremote_blk(1,1)  = cchild(1,5)
            fremote_pe(1,1)   = cchild(2,5)
            nblks_ind = 9
          elseif(ii.eq.1.and.jj.eq.1.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,4)
            fremote_pe(1,1)   = cchild(2,4)
            nblks_ind = 19
          elseif(ii.eq.3.and.jj.eq.1.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,3)
            fremote_pe(1,1)   = cchild(2,3)
            nblks_ind = 21
          elseif(ii.eq.1.and.jj.eq.3.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,2)
            fremote_pe(1,1)   = cchild(2,2)
            nblks_ind = 25
          elseif(ii.eq.3.and.jj.eq.3.and.kk.eq.3) then
            fremote_blk(1,1)  = cchild(1,1)
            fremote_pe(1,1)   = cchild(2,1)
            nblks_ind = 27
          endif

          call set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)


!
! If a neighbor exists at this blocks refinement level then fill guardcells
! from its data.

          idf = id
          jdf = jd
          kdf = kd
          if(ii.eq.3) idf = id + nxb
          if(jj.eq.3) jdf = jd + nyb
          if(kk.eq.3) kdf = kd + nzb


          if(lcc) call amr_1blk_cc_cp_remote(
     .                 mype,
     .                 fremote_pe(1,1),fremote_blk(1,1),
     .                 iblock,iopt,
     .                 idf,jdf,kdf,is,js,ks,ilays,jlays,klays,
     .                 nblks_ind)

          if(lfc) call amr_1blk_fc_cp_remote( 
     .                 mype,
     .                 fremote_pe(1,1),fremote_blk(1,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,ip1,jp1,kp1,
     .                 ip2,jp2,kp2,0,
     .                 nblks_ind)

          if(lec) call amr_1blk_ec_cp_remote( 
     .                 mype,
     .                 fremote_pe(1,1),fremote_blk(1,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,
     .                 0,0,0,0,0,0,0,
     .                 nblks_ind)

          if(ii.eq.3) ip1 = 1
          if(jj.eq.3) jp1 = 1
          if(kk.eq.3) kp1 = 1
          if(lnc) call amr_1blk_nc_cp_remote( 
     .                 mype,
     .                 fremote_pe(1,1),fremote_blk(1,1),
     .                 iblock,
     .                 idf,jdf,kdf,is,js,ks,
     .                 ilays,jlays,klays,
     .                 ip1,jp1,kp1,ip3,jp3,kp3,
     .                 nblks_ind)
          if(ii.eq.3) ip1 = 0           !reset to avoid errors during nn loop
          if(jj.eq.3) jp1 = 0
          if(kk.eq.3) kp1 = 0

        endif


      enddo
      enddo
      enddo

      endif                                 ! end of loop over 3D corners

!------------------------------------

      endif                                 ! end of ldiag if test

!------------------------------------

#ifdef DEBUG
      if(lcc.and.iopt.eq.1.and.lb.eq.12) then
!        write(*,*) 'Block ',lb,' unk1_fl'
!        do j= 1,(nguard02+nyb)*2
        do k = 7,8
        do j = 7,8
!         write(*,*) j,k,unk1_fl(1,:,j,k)
        enddo
        enddo
!        write(*,*) 'Block ',lb,' facevarx1_fl'
        do j= 1,(nguard02+nyb)*2
!         write(*,*) j,facevarx1_fl(1,:,j,1)
        enddo
!        write(*,*) 'Block ',lb,' unk_n1_fl'
!        do j= 1,(nguard02+nyb)*2+1
        do k=1,17
        do j=1,17
!         write(*,*) j,k,unk_n1_fl(1,:,j,k)
        enddo
        enddo
!        write(*,*) 'Block ',lb,' unk_e_x1_fl'
        do k=1,17
        do j=1,17
!         write(*,*) j,k,unk_e_x1_fl(1,:,j,k)
        enddo
        enddo
      endif
#endif /* DEBUG */


      return
      end subroutine amr_1blk_guardcell_f_to_c_set





      subroutine set_f2c_indexes(lcc,lfc,lec,lnc,nblks_ind,iopt)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit none

      logical, intent(in) :: lcc,lfc,lec,lnc
      integer, intent(in) :: nblks_ind,iopt

          if(lcc) then
          if(iopt.eq.1) then
            f2c_ind_unk(1,1,nblks_ind) = iu_bnd1
     .                        + nxb + nguard
            f2c_ind_unk(1,2,nblks_ind) = ju_bnd1
     .                        + (nyb + nguard)*k2d
            f2c_ind_unk(1,3,nblks_ind) = ku_bnd1
     .                        + (nzb + nguard)*k3d
            f2c_ind_unk(2,1,nblks_ind) = il_bnd1
            f2c_ind_unk(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unk(2,3,nblks_ind) = kl_bnd1
          else
            f2c_ind_work(1,1,nblks_ind) = iuw1
     .                        + nxb + nguard_work
            f2c_ind_work(1,2,nblks_ind) = juw1
     .                        + (nyb + nguard_work)*k2d
            f2c_ind_work(1,3,nblks_ind) = kuw1
     .                        + (nzb + nguard_work)*k3d
            f2c_ind_work(2,1,nblks_ind) = ilw1
            f2c_ind_work(2,2,nblks_ind) = jlw1
            f2c_ind_work(2,3,nblks_ind) = klw1
          endif
          endif
          if(lfc) then
            f2c_ind_facex(1,1,nblks_ind) = iu_bnd1 + 1
     .                        + nxb + nguard
            f2c_ind_facex(1,2,nblks_ind) = ju_bnd1
     .                        + (nyb + nguard)*k2d
            f2c_ind_facex(1,3,nblks_ind) = ku_bnd1
     .                        + (nzb + nguard)*k3d
            f2c_ind_facex(2,1,nblks_ind) = il_bnd1
            f2c_ind_facex(2,2,nblks_ind) = jl_bnd1
            f2c_ind_facex(2,3,nblks_ind) = kl_bnd1
            f2c_ind_facey(1,1,nblks_ind) = iu_bnd1 
     .                        + nxb + nguard
            f2c_ind_facey(1,2,nblks_ind) = ju_bnd1 + k2d
     .                        + (nyb + nguard)*k2d
            f2c_ind_facey(1,3,nblks_ind) = ku_bnd1
     .                        + (nzb + nguard)*k3d
            f2c_ind_facey(2,1,nblks_ind) = il_bnd1
            f2c_ind_facey(2,2,nblks_ind) = jl_bnd1
            f2c_ind_facey(2,3,nblks_ind) = kl_bnd1
            f2c_ind_facez(1,1,nblks_ind) = iu_bnd1 
     .                        + nxb + nguard
            f2c_ind_facez(1,2,nblks_ind) = ju_bnd1
     .                        + (nyb + nguard)*k2d
            f2c_ind_facez(1,3,nblks_ind) = ku_bnd1 + k3d
     .                        + (nzb + nguard)*k3d
            f2c_ind_facez(2,1,nblks_ind) = il_bnd1
            f2c_ind_facez(2,2,nblks_ind) = jl_bnd1
            f2c_ind_facez(2,3,nblks_ind) = kl_bnd1
          endif
          if(lec) then
            f2c_ind_unkex(1,1,nblks_ind) = iu_bnd1
     .                        + nxb + nguard
            f2c_ind_unkex(1,2,nblks_ind) = ju_bnd1 + k2d
     .                        + (nyb + nguard)*k2d
            f2c_ind_unkex(1,3,nblks_ind) = ku_bnd1 + k3d
     .                        + (nzb + nguard)*k3d
            f2c_ind_unkex(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkex(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkex(2,3,nblks_ind) = kl_bnd1
            f2c_ind_unkey(1,1,nblks_ind) = iu_bnd1 + 1
     .                        + nxb + nguard
            f2c_ind_unkey(1,2,nblks_ind) = ju_bnd1
     .                        + (nyb + nguard)*k2d
            f2c_ind_unkey(1,3,nblks_ind) = ku_bnd1 + k3d
     .                        + (nzb + nguard)*k3d
            f2c_ind_unkey(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkey(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkey(2,3,nblks_ind) = kl_bnd1
            f2c_ind_unkez(1,1,nblks_ind) = iu_bnd1 + 1
     .                        + nxb + nguard
            f2c_ind_unkez(1,2,nblks_ind) = ju_bnd1 + k2d
     .                        + (nyb + nguard)*k2d
            f2c_ind_unkez(1,3,nblks_ind) = ku_bnd1
     .                        + (nzb + nguard)*k3d
            f2c_ind_unkez(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkez(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkez(2,3,nblks_ind) = kl_bnd1
          endif
          if(lnc) then
            f2c_ind_unkn(1,1,nblks_ind) = iu_bnd1 + 1
     .                        + nxb + nguard
            f2c_ind_unkn(1,2,nblks_ind) = ju_bnd1 + k2d
     .                        + (nyb + nguard)*k2d
            f2c_ind_unkn(1,3,nblks_ind) = ku_bnd1 + k3d
     .                        + (nzb + nguard)*k3d
            f2c_ind_unkn(2,1,nblks_ind) = il_bnd1
            f2c_ind_unkn(2,2,nblks_ind) = jl_bnd1
            f2c_ind_unkn(2,3,nblks_ind) = kl_bnd1
          endif

          return

      end subroutine set_f2c_indexes
