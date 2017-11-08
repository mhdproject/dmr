#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#include "jet_conditions.fh"
      subroutine amr_1blk_bcset(mype,ibc,lb,pe,
     .    idest,iopt,iface0,jface0,kface0,surrblks)


!------------------------------------------------------------------------
!
! This routine sets guardcell values at external boundaries in the case
! where a single block is having its guardcells filled.
!
! It can be assumed in writing this routine, that all guardcells for this
! block which are not across an external boundary have already been
! properly filled.
!
! This routine provides code to set boundary conditions for a selection
! of common boundary conditions.
!
!   ibc = -21                       zero gradient
!   ibc = -22                       reflecting
!
!
! It is important to point out that any vector components stored
! in either unk or unk_n will have different reflection properties
! across different boundary planes.
! Note, in this example the following case is assumed:
!           unk(1,....)  is a scalar
!           unk(2,....)  is the x component of a vector
!           unk(3,....)  is the y component of a vector
!           unk(4,....)  is the z component of a vector
!           unk(5,....)  is a scalar
!
!           unk_n        all components are scalar.
!
! This is significant because it means that when ibc=-22, unk(2,...)
! changes it sign across x boundary planes but not across y or z
! boundary planes. Similarly unk(3,...) changes it sign across y 
! boundary planes but not across x or z boundary planes.
!
!------------------------------------------------------------------------
!
!
! Arguments:
!      mype             local processor
!      ibc              the integer specifying the particular boundary
!                        condition to be imposed
!      lb               block number of selected block
!      pe               processor on which block lb is located
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      iface            identifies location of cells to be filled with
!                        respect to the x boundaries, eg
!                        iface = -1      left x boundary guardcells
!                        iface =  0      x coords of guardcells located 
!                                        between left and right x boundary
!                        iface = +1      right x boundary guardcells
!
! Written :     Peter MacNeice          August 2001
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use athena_specific_parameters
      use physcons
      use strongfields

      implicit none

#include "amr_shmem.fh"

#ifdef MPI_USED
      include "mpif.h"
      integer,external :: shmem_my_pe
      integer          :: ierrorcode,ierr
#endif /* MPI_USED */

      real, save :: bbox(2,3)

      integer, intent(in) :: mype,ibc,lb,pe
      integer, intent(in) :: idest,iopt,iface0,jface0,kface0
      integer, intent(in) :: surrblks(:,:,:,:)

      integer :: iface, jface, kface,
     .     il_extra, iu_extra, jl_extra, ju_extra, kl_extra, ku_extra,
     .     nguard0, ilays, jlays, klays,  ip1, jp1, kp1, 
     .     id, jd, kd, il0, jl0, kl0,
     .     ndel
      
      real :: xmax, xmin, ymax, ymin, zmax, zmin, delx, dely, delz,
     .     xpos, ypos, zpos, rpos, distance, p, eth, rpos2,
     .     rposb, rpos2b, distanceb, distance1, distance2,
     .     maxv, x
      real :: variable_input_factor
      real :: jet_1_x_dist
      real :: jet_1_z_dist
      real :: jet_2_x_dist
      real :: jet_2_z_dist
      real :: jet1_x_orig
      real :: jet1_z_orig          
      real :: jet2_x_orig          
      real :: jet2_z_orig


          


      real ccoord(3),csize(3)
      save    ccoord,csize

      integer il,jl,kl,id1,jd1,kd1,i,j,k
      integer ks,js,is,lbw,ivar
      real :: fact

!-----x-direction from left supersonic inflow
!
      real :: b2,bv,bx,by,bz
      real :: rl,vxl,vyl,vzl,pl,el,bxl,byl,bzl

      common/work_bc_control/ ibc_work
      integer :: ibc_work

      integer :: lb0, ilb


      common/left_state/rl,vxl,vyl,vzl,pl,bxl,byl,bzl
      real :: rr,vx,vy,vz,bx1,by1,bz1,bx0,by0,bz0,v2,ek,bk,et,px,py,pz
      real :: time,t0,vl
      common/evolution/time
 
      real :: jet_den, jet_p, jet_vx,jet_vy,jet_vz,jet_radius, jet_beta,
     .     ambient_den, ambient_vx, ambient_vy, ambient_p
      common/jet_state/jet_den, jet_p, jet_vx,jet_vy,jet_vz,jet_radius, jet_beta,
     .     ambient_den, ambient_vx, ambient_vy, ambient_p

      if(mype.eq.2.and.lb.eq.5) then
c         write(*,*) mype, time,'amr_1blk in',unk1(1,11,5,1,idest)
      endif

!---------------------------------------------------------------------------
!      if(lb.eq.3) then
!     do k=2-k3d,2+k3d
!     do j=2-k2d,2+k2d
!     write(*,10) j,k, (surrblks(1,i,j,k),
!    & surrblks(2,i,j,k),i=1,3)
!0    format(" bcset: ",i3,2x,i3,3(2x,'(',i3,'/',i3,')'))
!     enddo
!     enddo

!      endif

      call strongfields_wb(pe,lb)




!---------------------------------------------------------------------------
! Do not modify this section 
      iface  = iface0
      jface  = jface0
      kface  = kface0

      if(iopt.eq.1) then
        nguard0 = nguard
      elseif(iopt.gt.1) then
        nguard0 = nguard_work
      endif

      ilays = nguard0
      jlays = nguard0*k2d
      klays = nguard0*k3d
      if(iface.eq.0) ilays = nxb
      if(jface.eq.0) jlays = nyb
      if(kface.eq.0) klays = nzb

      ip1 = 0
      jp1 = 0
      kp1 = 0
      if(iface.eq.1) ip1 = 1
      if(jface.eq.1) jp1 = k2d
      if(kface.eq.1) kp1 = k3d

      
      id = 1
      if(iface.eq.0)  id = 1+nguard0
      if(iface.eq.+1) id = nxb+1+nguard0
      jd = 1
      if(jface.eq.0)  jd = 1+nguard0*k2d
      if(jface.eq.+1) jd = nyb+(1+nguard0)*k2d
      kd = 1
      if(kface.eq.0)  kd = 1+nguard0*k3d
      if(kface.eq.+1) kd = nzb+(1+nguard0)*k3d

!
! Adjust index ranges
      il = ilays-1
      jl = (jlays-1)*k2d
      kl = (klays-1)*k3d


!
! Now reset iface,jface,kface so that blocks next to a boundary
! which treat their diagonal and corner guardcells correctly.
!     iface0 = iface
!     jface0 = jface
!     kface0 = kface
      if(iface.ne.0.and.surrblks(1,iface+2,2,2).gt.-20) iface=0
      if(jface.ne.0.and.surrblks(1,2,jface+2,2).gt.-20) jface=0
      if(kface.ne.0.and.surrblks(1,2,2,kface+2).gt.-20) kface=0




!---------------------------------------------------------------------------
! Section to be modified by user


! Which boundary condition has been specified?
       if(ibc.eq.-22) then
!
!-------------------------
! reflecting - all variables

        lbw = idest

        if(iopt.eq.1) then

!
! Do cell-face-centered data
        if(nfacevar.gt.0) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d

! set the following index range in facevarx1
!
!         facevarx1(:,id1:id1+il,jd:jd+jl,kd:kd+kl,idest) =
!     .            ????

          il0 = 0
          if(iface0.eq.0) il0 = 1

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id1 ,id1+il+il0
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+2
                     fact = -1.
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1)+gc_off_x)
                     fact = -1.
                 endif
                facevarx1(:nfacevar,i,j,k,lbw) = 
     .                    facevarx1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


! set the following index range in facevary1
!
!              facevary1(:,id:id+il,jd1:jd1+jl,kd:kd+kl,idest) =
!     .            ????

          jl0 = 0
          if(jface0.eq.0) jl0 = k2d

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd1,jd1+jl+jl0
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   fact = -1.
                   js = (2*nguard-j+gc_off_y+1)*k2d+1
               elseif(jface.eq.1) then
                   fact = -1.
                   js = (nyb+(nguard+1)*k2d)
     .                 -(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevary1(:nfacevar,i,j,k,lbw) =
     .                     facevary1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo

! set the following index range in facevarz1
!
!              facevarz1(:,id:id+il,jd:jd+jl,kd1:kd1+kl,idest) =
!     .            ????

          kl0 = 0
          if(kface0.eq.0) kl0 = k3d

          fact = 1.
          do k = kd1,kd1+kl+kl0
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               fact = -1.
               ks = (2*nguard-k+gc_off_z+1)*k3d+1
            elseif(kface.eq.1) then
               fact = -1.
               ks = (nzb+(nguard+1)*k3d)
     .             -(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevarz1(:nfacevar,i,j,k,lbw) =
     .                     facevarz1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


        endif                            ! end of nfacevar if test


!
! Now do cell centered data


        if(nvar.gt.0) then

! set the following index range in unk1
!
!          unk1(:,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard-(i-(nxb+nguard+1)+gc_off_x)
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                do ivar = 1,nvar
                    fact= 1.0
                  if(ivar.eq.2.and.abs(iface).ne.0)then
                    fact=-1.0
                  elseif(ivar.eq.3.and.abs(jface).ne.0)then
                    fact=-1.0
#if N_DIM == 3
                  elseif(ivar.eq.4.and.abs(kface).ne.0)then
                    fact=-1.0
#endif
                  endif

                  unk1(ivar,i,j,k,lbw) = fact*unk1(ivar,is,js,ks,lbw)


                enddo
              enddo
            enddo
          enddo

        endif                     ! end of nvar if test

!
! Now do cell corner data

        if(nvarcorn.gt.0) then

! set the appropriate index range in unk_n1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvarcorn
                   fact= 1.0
                   unk_n1(ivar,i,j,k,lbw) = 
     .                     fact*unk_n1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          endif                    ! end of nvarcorn if test

!
! Now do cell edge centered data

#ifdef N_DIM > 1

        if(nvaredge.gt.0) then
! First unk_e_x1
!
          il_extra = 0
          iu_extra = 0
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = nguard+(nguard+1-i)
                     fact = -1.
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard))
                     fact = -1.
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvaredge
                   unk_e_x1(ivar,i,j,k,lbw) = 
     .                     fact*unk_e_x1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

! Now unk_e_y1
!
          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = 1+(nguard-(j-nguard))*k2d
                  fact = -1.
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard))*k2d
                  fact = -1.
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (1+nguard)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (1+(nxb+nguard))-(i-(nxb+nguard+1))
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvaredge
                   unk_e_y1(ivar,i,j,k,lbw) = 
     .                     fact*unk_e_y1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

#if N_DIM == 3
! finally unk_e_z1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard-k)*k3d
                fact = -1.
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard))*k3d
                fact = -1.
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

! apply sign change to vector component perpendicular to the appropriate
! boundaries 
                 do ivar = 1,nvaredge
                   unk_e_z1(ivar,i,j,k,lbw) = 
     .                     fact*unk_e_z1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo
#endif /*  N_DIM == 3 */
#endif /*  N_DIM > 1 */

         endif                    ! end of nvaredge if test


        elseif(iopt.ge.2) then

! set the following index range in work1
!
!          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard_work-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard_work*k3d-
     .               (k-(nzb+nguard_work+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard_work-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard_work*k2d-
     .                   (j-(nyb+nguard_work+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard_work-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard_work-
     .                     (i-(nxb+nguard_work+1)+gc_off_x)
                 endif

                work1(i,j,k,lbw) = work1(is,js,ks,lbw)


              enddo
            enddo
          enddo


        endif


!-------------------------


      elseif(ibc.eq.-21) then


!-------------------------
! zero gradient - all variables

        lbw = idest

        if(iopt.eq.1) then

!
! Do cell-face-centered data
        if(nfacevar.gt.0) then

          id1 = id + ip1
          jd1 = jd + jp1*k2d
          kd1 = kd + kp1*k3d

! set the following index range in facevarx1
!
!         facevarx1(:,id1:id1+il,jd:jd+jl,kd:kd+kl,idest) =
!     .            ????

          il0 = 0
          if(iface0.eq.0) il0 = 1

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id1 ,id1+il+il0
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+2
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1)+gc_off_x)
                 endif
                facevarx1(:nfacevar,i,j,k,lbw) = 
     .                    facevarx1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


! set the following index range in facevary1
!
!              facevary1(:,id:id+il,jd1:jd1+jl,kd:kd+kl,idest) =
!     .            ????

          jl0 = 0
          if(jface0.eq.0) jl0 = k2d

          fact = 1.
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd1,jd1+jl+jl0
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y+1)*k2d+1
               elseif(jface.eq.1) then
                   js = (nyb+(nguard+1)*k2d)
     .                 -(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevary1(:nfacevar,i,j,k,lbw) =
     .                     facevary1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo

! set the following index range in facevarz1
!
!              facevarz1(:,id:id+il,jd:jd+jl,kd1:kd1+kl,idest) =
!     .            ????

          kl0 = 0
          if(kface0.eq.0) kl0 = k3d

          fact = 1.
          do k = kd1,kd1+kl+kl0
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z+1)*k3d+1
            elseif(kface.eq.1) then
               ks = (nzb+(nguard+1)*k3d)
     .             -(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = (nxb+nguard)-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                facevarz1(:nfacevar,i,j,k,lbw) =
     .                     facevarz1(:nfacevar,is,js,ks,lbw)*fact
              enddo
            enddo
          enddo


        endif                            ! end of nfacevar if test


!
! Now do cell centered data


        if(nvar.gt.0) then

! set the following index range in unk1
!
!          unk1(:,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard*k3d-(k-(nzb+nguard+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard*k2d-(j-(nyb+nguard+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard-(i-(nxb+nguard+1)+gc_off_x)
                 endif

                do ivar = 1,nvar
                  fact= 1.0
                  unk1(ivar,i,j,k,lbw) = fact*unk1(ivar,is,js,ks,lbw)
!                 if(lb.eq.3.and.ivar.eq.2.and.j.eq.4)then
!                   write(*,*)'i,j,fact,is,js,unk1 ',i,j,fact,is,js,
!    &          unk1(ivar,is,js,ks,lbw)
!                 endif

                enddo
              enddo
            enddo
          enddo

        endif                     ! end of nvar if test

!
! Now do cell corner data

        if(nvarcorn.gt.0) then

! set the appropriate index range in unk_n1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

                 do ivar = 1,nvarcorn
                   fact= 1.0
                   unk_n1(ivar,i,j,k,lbw) = 
     .                     fact*unk_n1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

          endif                    ! end of nvarcorn if test

!
! Now do cell edge centered data

#ifdef N_DIM > 1

        if(nvaredge.gt.0) then
! First unk_e_x1
!
          il_extra = 0
          iu_extra = 0
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = nguard+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard))
                 endif

                 do ivar = 1,nvaredge
                   unk_e_x1(ivar,i,j,k,lbw) = 
     .                     fact*unk_e_x1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

! Now unk_e_y1
!
          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          kl_extra = 0
          ku_extra = 0
          if(kface0.eq.0) then
            ku_extra = k3d
          elseif(kface0.eq.+1) then
            kl_extra = k3d
            ku_extra = k3d
          endif

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard+1-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard+1))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = 1+(nguard-(j-nguard))*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (1+nguard)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (1+(nxb+nguard))-(i-(nxb+nguard+1))
                 endif

                 do ivar = 1,nvaredge
                   unk_e_y1(ivar,i,j,k,lbw) = 
     .                     fact*unk_e_y1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo

#if N_DIM == 3
! finally unk_e_z1
!

          il_extra = 0
          iu_extra = 0
          if(iface0.eq.0) then
            iu_extra = 1
          elseif(iface0.eq.+1) then
            il_extra = 1
            iu_extra = 1
          endif
          jl_extra = 0
          ju_extra = 0
          if(jface0.eq.0) then
            ju_extra = k2d
          elseif(jface0.eq.+1) then
            jl_extra = k2d
            ju_extra = k2d
          endif
          kl_extra = 0
          ku_extra = 0

          fact = 1.

          lbw = idest

          do k = kd+kl_extra,kd+kl+ku_extra
            if(kface.eq.0 ) then
                ks = k 
            elseif(kface.eq.-1) then
                ks = (1+nguard*k3d)+(nguard-k)*k3d
            elseif(kface.eq.1) then
                ks = (1+(nzb+nguard)*k3d)-(k-(nzb+nguard))*k3d
            endif

            do j = jd+jl_extra,jd+jl+ju_extra
              if(jface.eq.0 ) then
                  js = j 
              elseif(jface.eq.-1) then
                  js = (1+nguard*k2d)+(nguard+1-j)*k2d
              elseif(jface.eq.1) then
                  js = (1+(nyb+nguard)*k2d)-(j-(nyb+nguard+1))*k2d
              endif

              do i = id+il_extra,id+il+iu_extra
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (nguard+1)+(nguard+1-i)
                 elseif(iface.eq.1) then
                     is = (nxb+nguard+1)-(i-(nxb+nguard+1))
                 endif

                 do ivar = 1,nvaredge
                   unk_e_z1(ivar,i,j,k,lbw) = 
     .                     fact*unk_e_z1(ivar,is,js,ks,lbw)
                 enddo

              enddo
            enddo
          enddo
#endif /*  N_DIM == 3 */
#endif /*  N_DIM > 1 */

         endif                    ! end of nvaredge if test


        elseif(iopt.ge.2) then

! set the following index range in work1
!
!          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

          lbw = idest
          do k = kd,kd+kl
            if(kface.eq.0 ) then
               ks = k
            elseif(kface.eq.-1) then
               ks = (2*nguard_work-k+gc_off_z)*k3d+1
            elseif(kface.eq.1) then
               ks = nzb+nguard_work*k3d-
     .               (k-(nzb+nguard_work+1)+gc_off_z)*k3d
            endif
            do j = jd,jd+jl
               if(jface.eq.0 ) then
                   js = j 
               elseif(jface.eq.-1) then
                   js = (2*nguard_work-j+gc_off_y)*k2d+1
               elseif(jface.eq.1) then
                   js = nyb+nguard_work*k2d-
     .                   (j-(nyb+nguard_work+1)+gc_off_y)*k2d
               endif
              do i = id,id+il
                 if(iface.eq.0 ) then
                     is = i 
                 elseif(iface.eq.-1) then
                     is = (2*nguard_work-i+gc_off_x)+1
                 elseif(iface.eq.1) then
                     is = nxb+nguard_work-
     .                     (i-(nxb+nguard_work+1)+gc_off_x)
                 endif

                work1(i,j,k,lbw) = work1(is,js,ks,lbw)


              enddo
            enddo
          enddo


        endif


!
!-------------------------

      elseif(ibc.eq.-30) then
!
!-------------------------
! supersonic inflow from left x-boundary
c          write(*,*)'entering ibc.eq.-30' 
          lbw = idest
#ifdef MAGNETOSPHERE
          t0  = 50./0.4 ! time scale for ramping up inflow velocity
#else
          t0  = 100./0.4 ! time scale for ramping up inflow velocity
#endif
!pmn      vl  = vxl*(1.-exp(-time/t0))
          vl  = vxl
          b2  = bxl**2+byl**2+bzl**2 
          bv  = (bxl*vl+byl*vyl+bzl*vzl)
          el  =  0.5*rl*vl**2+pl*gammam1i+pie8i*b2
c          write(*,*)'bxl,byl,bzl=',bxl,byl,bzl
c          write(*,*)'vxl,vyl,vzl=',vxl,vyl,vzl
!
! Now do cell centered data
          lbw = idest
          do k = kl_bnd1,ku_bnd1
           do j = jl_bnd1,ju_bnd1
             do i = 1,nguard
!
!------first convert to primitives
!
                rr                 = rl
!pmn            vx                 = vxl*(1.-exp(-time/t0))
                vx                 = vxl
                vy                 = vyl
                vz                 = vzl
                bx0                = b_sf(1,i,j,k)
                by0                = b_sf(2,i,j,k)
                bz0                = b_sf(3,i,j,k)
                bx1                = bxl-bx0
                by1                = byl-by0
                bz1                = bzl-bz0
                bx                 = bx1
                by                 = by1
                bz                 = bz1
                v2                 = vx*vx+vy*vy+vz*vz
                bv                 = (bx*vx+by*vy+bz*vz)
                b2                 = bx*bx+by*by+bz*bz
                ek                 = 0.5*rl*vx*vx
                bk                 = (bxl*bxl+byl*byl+bzl*bzl)*pie8i
                et                 = ek + pl*gammam1i+bk 
                px                 = rr*vx + afact*(b2*vx - bv*bx)
                py                 = rr*vy + afact*(b2*vy - bv*by)
                pz                 = rr*vz + afact*(b2*vz - bv*bz)
                px                 = rl*vx
                py                 = rl*vy
                pz                 = rl*vz

              if(iopt.eq.1) then
                unk1( 1,i,j,k,lbw) = rl
                unk1( 2,i,j,k,lbw) = px  
                unk1( 3,i,j,k,lbw) = py 
                unk1( 4,i,j,k,lbw) = pz 
                unk1( 5,i,j,k,lbw) = et
                unk1( 6,i,j,k,lbw) = bxl-b_sf(1,i,j,k)
                unk1( 7,i,j,k,lbw) = byl-b_sf(2,i,j,k)
                unk1( 8,i,j,k,lbw) = bzl-b_sf(3,i,j,k)
c                unk1( 6,i,j,k,lbw) = bxl
c                unk1( 7,i,j,k,lbw) = byl
c                unk1( 8,i,j,k,lbw) = bzl
              elseif(iopt.gt.1) then
                if(ibc_work.eq.1) then
                  work1(i,j,k,iopt-1) = rl
                endif
              endif

#ifdef DEBUGG
            if(i.eq.nguard) then
              write(*,*)'i,j,nguard        =',i,j,nguard
              write(*,*)'u(1,ix,iy,iz) = ',unk1(1,i,j,k,lbw)
              write(*,*)'u(2,ix,iy,iz) = ',unk1(2,i,j,k,lbw)
              write(*,*)'u(3,ix,iy,iz) = ',unk1(3,i,j,k,lbw)
              write(*,*)'u(4,ix,iy,iz) = ',unk1(4,i,j,k,lbw)
              write(*,*)'u(5,ix,iy,iz) = ',unk1(5,i,j,k,lbw)
              write(*,*)'u(6,ix,iy,iz) = ',unk1(6,i,j,k,lbw)
              write(*,*)'u(8,ix,iy,iz) = ',unk1(8,i,j,k,lbw)
              write(*,*)'u(7,ix,iy,iz) = ',unk1(7,i,j,k,lbw)
              write(*,*)'bxl,b_sf(1,ix,iy,iz) = ',bxl,bx0
              write(*,*)'byl,b_sf(2,ix,iy,iz) = ',byl,by0
              write(*,*)'bzl,b_sf(3,ix,iy,iz) = ',bzl,bz0
c              pause
            endif
#endif /*DEBUGG*/

              enddo
             enddo
           enddo
#ifdef XWIND
#ifdef NOLEFTSTATE
      do k=kl_bnd,ku_bnd
       do j=jl_bnd,ju_bnd
        do i=1,nguard
        facevarx1(1,i,j,k,lbw)=0.0
        enddo   
       enddo   
      enddo
!
!     y faces
!         
      do k=kl_bnd,ku_bnd
       do j=jl_bnd,ju_bnd+k2d
        do i=1,nguard
        facevary1(1,i,j,k,lbw)=0.0
        enddo   
       enddo   
      enddo
!
!     z faces
!         
      do k=kl_bnd,ku_bnd+k3d
       do j=jl_bnd,ju_bnd
        do i=1,nguard
        facevarz1(1,i,j,k,lbw)=0.0
        enddo   
       enddo   
      enddo
#else
      do k=kl_bnd1,ku_bnd1
       do j=jl_bnd1,ju_bnd1
        do i=1,nguard
         facevarx1(1,i,j,k,lbw)=bxl-bxf_sf(1,i,j,k)
        enddo   
       enddo   
      enddo
!
!     y faces
!         
      do k=kl_bnd1,ku_bnd1
       do j=jl_bnd1,ju_bnd1+k2d
        do i=1,nguard
        facevary1(1,i,j,k,lbw)=byl-byf_sf(1,i,j,k)
c        write(*,*)'byl,byf_sf(1,i,j,k)=',byl,byf_sf(1,i,j,k)
        enddo   
       enddo   
      enddo
!
!     z faces
!         
      do k=kl_bnd1,ku_bnd1+k3d
       do j=jl_bnd1,ju_bnd1
        do i=1,nguard
        facevarz1(1,i,j,k,lbw)=bzl-bzf_sf(1,i,j,k)
        enddo   
       enddo   
      enddo
#endif
#endif

!
!-------------------------

      elseif(ibc.eq.-40) then


cc lower y boundary - up to 1/6 is postshock - the rest is reflecting
!
!-------------------------

         lb0=lb
#ifdef MPI_USED
      if(pe.ne.mype) then

         lb0=0;  ilb=strt_buffer;

         do while( (ilb.ge.strt_buffer.and.
     &              ilb.le.last_buffer).and.
     &              lb0.eq.0) 
            if(lb.eq.laddress(1,ilb).and.pe.eq.laddress(2,ilb)) then
                lb0=ilb
             endif
             ilb=ilb+1
          end do

         if(lb0.eq.0) then
            write(*,*) 'Error: amr_1blk_bcset : could not find block'
            call amr_close
            stop
         endif
      endif
#endif
      call shmem_real_get(bbox(1,1),bnd_box(1,1,lb0),6,mype)


      xmin  = bbox(1,1)
      xmax  = bbox(2,1)
      delx  = (xmax - xmin)/float(nxb)

      ymin  = bbox(1,2)
      ymax  = bbox(2,2)
      dely  = (ymax - ymin)/float(nyb)

#if N_DIM == 3
      zmin  = bbox(1,3)
      zmax  = bbox(2,3)
      delz  = (zmax - zmin)/float(nzb)
#else
      zmax  = 0.
      zmin  = 0.
      delz  = 0.
#endif /*NDIM3*/




      ambient_p=AMBIENT_P
      ambient_den=AMBIENT_DEN
      ambient_vx=AMBIENT_VX
      ambient_vy=AMBIENT_VY

!

ciic      variable_input_factor = (1- (0.15*sin(0.785*time)))
      variable_input_factor = 1.

      jet_p = JET_P
      jet_vx = JET_VX
      jet_vy = JET_VY
      jet_vz = JET_VY
      jet_den = JET_DEN


! Now do cell centered data
ciic
      jet1_x_orig  =  JET_X_0 
      maxv=0.0
      x=0.
      if(time.gt.x) then
         jet1_x_orig  =  JET_X_0  -0.1*(time-x)*maxv
      endif
c      jet1_x_orig  =  JET_X_0  - 0.95 * time * (delx/2 ** ( lrefine_max - lrefine(lb)))/1.e-3  
c      write (*,*) time , delx , delx/ 2 ** (lrefine_max - lrefine (lb))
c      jet1_x_orig  =  JET_X_0  * cos ( 0.785 * time/40)
c      jet1_z_orig  =  JET_X_0  * sin ( 0.785 * time/40)

      jet2_x_orig  = -JET_X_0  
c      jet2_x_orig  = -JET_X_0  * cos ( 0.785 * time/40)
c      jet2_z_orig  = -JET_X_0  * sin ( 0.785 * time/40)

      jet1_z_orig  =  JET_X_0 * 0

      jet2_z_orig  = -JET_X_0 * 0


      if(iopt.eq.1) then
! set the following index range in unk1
!          unk1(:,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 

          lbw = idest
          do k = kd,kd+kl
             do i = id,id+il
                do j = jd,jd+jl
                   js = 1+2*nguard -j
                   xpos = xmin + delx*(i-nguard-0.5)
#if N_DIM == 3
                   zpos = zmin + delz*(k-nguard-0.5)
#else 
                   zpos = 0.
#endif
               jet_1_x_dist = sqrt ((xpos -jet1_x_orig  )**2)
               jet_2_x_dist = sqrt ((xpos -jet2_x_orig  )**2)
               jet_1_z_dist = sqrt ((zpos -jet1_z_orig  )**2)
               jet_2_z_dist = sqrt ((zpos -jet2_z_orig  )**2)

                   rpos2 = ((xpos-jet1_x_orig   )**2) + ((zpos-jet1_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance1=rpos

                   rpos2 = ((xpos-jet2_x_orig   )**2) + ((zpos-jet2_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance2=rpos


       
ciic distanceb????


cc lower y boundary - up to 1/6 is postshock - the rest is reflecting

              if(xpos .lt. ( JET_XMAX /24) )then

                      
                      rr=ambient_den
                      p=ambient_p
                      ek=0.5*rr*(  ambient_vx*ambient_vx + ambient_vy*ambient_vy)
                      eth=p/(gamma-1.)
                      
                      unk1(:,i,j,k,lbw) = 0.
                     
                      unk1( 1,i,j,k,lbw) = rr
                      unk1( 2,i,j,k,lbw) = rr *ambient_vx
                      unk1( 3,i,j,k,lbw) = rr *ambient_vy
                      unk1( 4,i,j,k,lbw) = 0.0
                      unk1( 5,i,j,k,lbw) = ek+eth ! no fields so far
                      
                   else         ! reflecting (zero inflow)
                      

                      unk1(:,i,j,k,lbw) =      unk1(:,i,js,k,lbw)
                      unk1(3,i,j,k,lbw) =     -unk1(3,i,js,k,lbw)
ceec                      unk1(3,i,j,k,lbw) =     max(0.0,0.5*unk1(3,i,js,k,lbw))


                      
                   endif
                enddo   
             enddo   
          enddo
          
       elseif(iopt.gt.1) then


! set the following index range in work1
!
!          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

c         write(*,*) 'iopt,ibc,id,id+il',iopt,ibc,id,id+il
c         write(*,*) '                 ',iopt,ibc,jd,jd+jl

         lbw = idest
         do k = kd,kd+kl
            do j = jd,jd+jl
               do i = id,id+il
                  xpos = xmin + delx*(i-nguard-0.5)
#if N_DIM == 3
                  zpos = zmin + delz*(k-nguard-0.5)
#else  
                  zpos = 0.
#endif

                   rpos2 = ((xpos-jet1_x_orig   )**2) + ((zpos-jet1_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance1=rpos

                   rpos2 = ((xpos-jet2_x_orig   )**2) + ((zpos-jet2_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance2=rpos
ciic distanceb????
c                   if((distance1 .lt. jet_radius).or.(distance2 .lt. jet_radius))then ! jet aperture
                   if((distance1 .lt. jet_radius))then ! jet aperture
                     rr=jet_den
                     work1(i,j,k,lbw) = rr
                  else
                     js = 1+2*nguard_work -j
                     work1(i,j,k,lbw) = work1(i,js,k,lbw)
                  endif
               enddo
            enddo
         enddo
      endif
                

!
!     y faces
!         
      do k = 1,nz
         do i = 1,nx
            do j = 1,nguard
               js = 2+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-0.5)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-0.5)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               if(distance .lt. jet_radius)then ! jet aperture
                  facevary1(1,i,j,k,lbw)=0.0
               else
                  facevary1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo


!
!     x faces
!         
      do k = 1,nz
         do i = 1,nx+1
            do j = 1,nguard
               js = 1+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-1.0)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-0.5)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               if(distance .lt. jet_radius)then ! jet aperture
                  facevarx1(1,i,j,k,lbw)=0.0
               else
                  facevarx1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo

!
!     z faces
!         
#if N_DIM == 3
      do k = 1,nz+1
         do i = 1,nx
            do j = 1,nguard
               js = 1+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-0.5)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-1.0)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               if(distance .lt. jet_radius)then ! jet aperture
                  facevarz1(1,i,j,k,lbw)=0.0
               else
                  facevarz1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo
#endif

!
!-------------------------
      elseif(ibc.eq.-50) then
!
!-------------------------
! the moving upper boundary - time dependent
         lb0=lb
#ifdef MPI_USED
      if(pe.ne.mype) then

         lb0=0;  ilb=strt_buffer;

         do while( (ilb.ge.strt_buffer.and.
     &              ilb.le.last_buffer).and.
     &              lb0.eq.0) 
            if(lb.eq.laddress(1,ilb).and.pe.eq.laddress(2,ilb)) then
                lb0=ilb
             endif
             ilb=ilb+1
          end do

         if(lb0.eq.0) then
            write(*,*) 'Error: amr_1blk_bcset : could not find block'
            call amr_close
            stop
         endif
      endif
#endif
      call shmem_real_get(bbox(1,1),bnd_box(1,1,lb0),6,mype)


      xmin  = bbox(1,1)
      xmax  = bbox(2,1)
      delx  = (xmax - xmin)/float(nxb)

      ymin  = bbox(1,2)
      ymax  = bbox(2,2)
      dely  = (ymax - ymin)/float(nyb)

#if N_DIM == 3
      zmin  = bbox(1,3)
      zmax  = bbox(2,3)
      delz  = (zmax - zmin)/float(nzb)
#else
      zmax  = 0.
      zmin  = 0.
      delz  = 0.
#endif /*NDIM3*/

      ambient_p=AMBIENT_P
      ambient_den=AMBIENT_DEN
      ambient_vx=AMBIENT_VX
      ambient_vy=AMBIENT_VY
ciic      variable_input_factor = (1- (0.15*sin(0.785*time)))
      variable_input_factor = 1.

      jet_p = JET_P
      jet_vx = JET_VX
      jet_vy = JET_VY
      jet_vz = JET_VZ
      jet_den = JET_DEN


! Now do cell centered data
ciic
      jet1_x_orig  =  JET_X_0 
      maxv=0.0
      x=0.
      if(time.gt.x) then
         jet1_x_orig  =  JET_X_0  -0.1*(time-x)*maxv
      endif

      jet2_x_orig  = -JET_X_0  

      jet1_z_orig  =  JET_X_0 * 0

      jet2_z_orig  = -JET_X_0 * 0


      if(iopt.eq.1) then
! set the following index range in unk1
!          unk1(:,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 

          lbw = idest
          do k = kd,kd+kl
             do i = id,id+il
                do j = jd,jd+jl
                   js = 1+2*nguard -j
                   xpos = xmin + delx*(i-nguard-0.5)
#if N_DIM == 3
                   zpos = zmin + delz*(k-nguard-0.5)
#else 
                   zpos = 0.
#endif
               jet_1_x_dist = sqrt ((xpos -jet1_x_orig  )**2)
               jet_2_x_dist = sqrt ((xpos -jet2_x_orig  )**2)
               jet_1_z_dist = sqrt ((zpos -jet1_z_orig  )**2)
               jet_2_z_dist = sqrt ((zpos -jet2_z_orig  )**2)

                   rpos2 = ((xpos-jet1_x_orig   )**2) + ((zpos-jet1_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance1=rpos

                   rpos2 = ((xpos-jet2_x_orig   )**2) + ((zpos-jet2_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance2=rpos


       
ciic distanceb????

              if(xpos .lt. (JET_XMAX/24 + JET_YMAX/sqrt(3.) + (10 * time /0.866)) )then

                      
                      rr=ambient_den
                      p=ambient_p
                      ek=0.5*rr*(  ambient_vx*ambient_vx + ambient_vy*ambient_vy)
                      eth=p/(gamma-1.)
                      
                      unk1(:,i,j,k,lbw) = 0.
                      unk1( 1,i,j,k,lbw) = rr
                      unk1( 2,i,j,k,lbw) = rr *ambient_vx
                      unk1( 3,i,j,k,lbw) = rr *ambient_vy
                      unk1( 4,i,j,k,lbw) = 0.0
                      unk1( 5,i,j,k,lbw) = ek+eth ! no fields so far
                      
                   else         ! reflecting (zero inflow)
                      

                      rr=jet_den
                      p=jet_p
                      ek=0.5*rr*(  jet_vx*jet_vx + jet_vy*jet_vy)
                      eth=p/(gamma-1.)
                      
                      unk1(:,i,j,k,lbw) = 0.
                      unk1( 1,i,j,k,lbw) = rr
                      unk1( 2,i,j,k,lbw) = rr *jet_vx
                      unk1( 3,i,j,k,lbw) = rr *jet_vy
                      unk1( 4,i,j,k,lbw) = 0.0
                      unk1( 5,i,j,k,lbw) = ek+eth ! no fields so far

                      
                   endif
                enddo   
             enddo   
          enddo
          
       elseif(iopt.gt.1) then


! set the following index range in work1
!
!          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

c         write(*,*) 'iopt,ibc,id,id+il',iopt,ibc,id,id+il
c         write(*,*) '                 ',iopt,ibc,jd,jd+jl

         lbw = idest
         do k = kd,kd+kl
            do j = jd,jd+jl
               do i = id,id+il
                  xpos = xmin + delx*(i-nguard-0.5)
#if N_DIM == 3
                  zpos = zmin + delz*(k-nguard-0.5)
#else  
                  zpos = 0.
#endif

                   rpos2 = ((xpos-jet1_x_orig   )**2) + ((zpos-jet1_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance1=rpos

                   rpos2 = ((xpos-jet2_x_orig   )**2) + ((zpos-jet2_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance2=rpos
ciic distanceb????
c                   if((distance1 .lt. jet_radius).or.(distance2 .lt. jet_radius))then ! jet aperture
c              if(xpos .lt. (16.666667 + 57.735 + (10 * time /0.866)) )then
              if(xpos .lt. (JET_XMAX/24 + JET_YMAX/sqrt(3.) + (10 * time /0.866)) )then
c                   if((distance1 .lt. jet_radius))then ! jet aperture
                     rr=jet_den
                     work1(i,j,k,lbw) = rr
                  else
                     js = 1+2*nguard_work -j
                     work1(i,j,k,lbw) = work1(i,js,k,lbw)
                  endif
               enddo
            enddo
         enddo
      endif
                

!
!     y faces
!         
      do k = 1,nz
         do i = 1,nx
            do j = 1,nguard
               js = 2+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-0.5)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-0.5)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
c               if(distance .lt. jet_radius)then ! jet aperture
c              if(xpos .lt. (16.666667 + 57.735 + (10 * time /0.866)) )then
              if(xpos .lt. (JET_XMAX/24 + JET_YMAX/sqrt(3.) + (10 * time /0.866)) )then
                  facevary1(1,i,j,k,lbw)=0.0
               else
                  facevary1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo


!
!     x faces
!         
      do k = 1,nz
         do i = 1,nx+1
            do j = 1,nguard
               js = 1+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-1.0)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-0.5)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
c               if(distance .lt. jet_radius)then ! jet aperture
c              if(xpos .lt. (16.666667 + 57.735 + (10 * time /0.866)) )then
              if(xpos .lt. (JET_XMAX/24 + JET_YMAX/sqrt(3.) + (10 * time /0.866)) )then
                  facevarx1(1,i,j,k,lbw)=0.0
               else
                  facevarx1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo

!
!     z faces
!         
#if N_DIM == 3
      do k = 1,nz+1
         do i = 1,nx
            do j = 1,nguard
               js = 1+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-0.5)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-1.0)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
c               if(distance .lt. jet_radius)then ! jet aperture
c              if(xpos .lt. (16.666667 + 57.735 + (10 * time /0.866)) )then
              if(xpos .lt. (JET_XMAX/24 + JET_YMAX/sqrt(3.) + (10 * time /0.866)) )then
                  facevarz1(1,i,j,k,lbw)=0.0
               else
                  facevarz1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo
#endif

!
!-------------------------

!-------------------------
      elseif(ibc.eq.-38) then
! lower x boundary - left boundary - set to post shock values 

! sets it to ambient values = postshock values 
!
!-------------------------
         lb0=lb
#ifdef MPI_USED
      if(pe.ne.mype) then

         lb0=0;  ilb=strt_buffer;

         do while( (ilb.ge.strt_buffer.and.
     &              ilb.le.last_buffer).and.
     &              lb0.eq.0) 
            if(lb.eq.laddress(1,ilb).and.pe.eq.laddress(2,ilb)) then
                lb0=ilb
             endif
             ilb=ilb+1
          end do

         if(lb0.eq.0) then
            write(*,*) 'Error: amr_1blk_bcset : could not find block'
            call amr_close
            stop
         endif
      endif
#endif
      call shmem_real_get(bbox(1,1),bnd_box(1,1,lb0),6,mype)


      xmin  = bbox(1,1)
      xmax  = bbox(2,1)
      delx  = (xmax - xmin)/float(nxb)

      ymin  = bbox(1,2)
      ymax  = bbox(2,2)
      dely  = (ymax - ymin)/float(nyb)

#if N_DIM == 3
      zmin  = bbox(1,3)
      zmax  = bbox(2,3)
      delz  = (zmax - zmin)/float(nzb)
#else
      zmax  = 0.
      zmin  = 0.
      delz  = 0.
#endif /*NDIM3*/

      ambient_p=AMBIENT_P
      ambient_den=AMBIENT_DEN
      ambient_vx=AMBIENT_VX
      ambient_vy=AMBIENT_VY
ciic      variable_input_factor = (1- (0.15*sin(0.785*time)))
      variable_input_factor = 1.

      jet_p = JET_P
      jet_vx = JET_VX
      jet_vy = JET_VY
      jet_vz = JET_VZ
      jet_den = JET_DEN


! Now do cell centered data
ciic
      jet1_x_orig  =  JET_X_0 
      maxv=0.0
      x=0.
      if(time.gt.x) then
         jet1_x_orig  =  JET_X_0  -0.1*(time-x)*maxv
      endif

      jet2_x_orig  = -JET_X_0  

      jet1_z_orig  =  JET_X_0 * 0

      jet2_z_orig  = -JET_X_0 * 0


      if(iopt.eq.1) then
! set the following index range in unk1
!          unk1(:,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 

          lbw = idest
          do k = kd,kd+kl
             do i = id,id+il
                do j = jd,jd+jl
                   js = 1+2*nguard -j
                   xpos = xmin + delx*(i-nguard-0.5)
#if N_DIM == 3
                   zpos = zmin + delz*(k-nguard-0.5)
#else 
                   zpos = 0.
#endif
               jet_1_x_dist = sqrt ((xpos -jet1_x_orig  )**2)
               jet_2_x_dist = sqrt ((xpos -jet2_x_orig  )**2)
               jet_1_z_dist = sqrt ((zpos -jet1_z_orig  )**2)
               jet_2_z_dist = sqrt ((zpos -jet2_z_orig  )**2)

                   rpos2 = ((xpos-jet1_x_orig   )**2) + ((zpos-jet1_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance1=rpos

                   rpos2 = ((xpos-jet2_x_orig   )**2) + ((zpos-jet2_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance2=rpos


       
ciic distanceb????


                      
                      rr=ambient_den
                      p=ambient_p
                      ek=0.5*rr*(  ambient_vx*ambient_vx + ambient_vy*ambient_vy)
                      eth=p/(gamma-1.)
                      
                      unk1(:,i,j,k,lbw) = 0.
                      unk1( 1,i,j,k,lbw) = rr
                      unk1( 2,i,j,k,lbw) = rr *ambient_vx
                      unk1( 3,i,j,k,lbw) = rr *ambient_vy
                      unk1( 4,i,j,k,lbw) = 0.0
                      unk1( 5,i,j,k,lbw) = ek+eth ! no fields so far
                      
                      


                      
                enddo   
             enddo   
          enddo
          
       elseif(iopt.gt.1) then


! set the following index range in work1
!
!          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
!     .             ????

c         write(*,*) 'iopt,ibc,id,id+il',iopt,ibc,id,id+il
c         write(*,*) '                 ',iopt,ibc,jd,jd+jl

         lbw = idest
         do k = kd,kd+kl
            do j = jd,jd+jl
               do i = id,id+il
                  xpos = xmin + delx*(i-nguard-0.5)
#if N_DIM == 3
                  zpos = zmin + delz*(k-nguard-0.5)
#else  
                  zpos = 0.
#endif

                   rpos2 = ((xpos-jet1_x_orig   )**2) + ((zpos-jet1_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance1=rpos

                   rpos2 = ((xpos-jet2_x_orig   )**2) + ((zpos-jet2_z_orig  )**2)
                   rpos = sqrt(rpos2)
                   distance2=rpos
ciic distanceb????
c                   if((distance1 .lt. jet_radius).or.(distance2 .lt. jet_radius))then ! jet aperture
                   if((distance1 .lt. jet_radius))then ! jet aperture
                     rr=jet_den
                     work1(i,j,k,lbw) = rr
                  else
                     js = 1+2*nguard_work -j
                     work1(i,j,k,lbw) = work1(i,js,k,lbw)
                  endif
               enddo
            enddo
         enddo
      endif
                

!
!     y faces
!         
      do k = 1,nz
         do i = 1,nx
            do j = 1,nguard
               js = 2+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-0.5)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-0.5)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               if(distance .lt. jet_radius)then ! jet aperture
                  facevary1(1,i,j,k,lbw)=0.0
               else
                  facevary1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo


!
!     x faces
!         
      do k = 1,nz
         do i = 1,nx+1
            do j = 1,nguard
               js = 1+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-1.0)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-0.5)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               if(distance .lt. jet_radius)then ! jet aperture
                  facevarx1(1,i,j,k,lbw)=0.0
               else
                  facevarx1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo

!
!     z faces
!         
#if N_DIM == 3
      do k = 1,nz+1
         do i = 1,nx
            do j = 1,nguard
               js = 1+2*nguard -j
   
               xpos = xmin + delx*(i-nguard-0.5)

#if N_DIM == 3
               zpos = zmin + delz*(k-nguard-1.0)
#else
               zpos = 0.
#endif
               rpos2 = (xpos**2) + (zpos**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               if(distance .lt. jet_radius)then ! jet aperture
                  facevarz1(1,i,j,k,lbw)=0.0
               else
                  facevarz1(1,i,j,k,lbw)=0.0
               endif
            enddo   
         enddo   
      enddo
#endif

      endif                     ! end of test of bc flag


! End of Section to be modified by user
!---------------------------------------------------------------------------

      if(mype.eq.2.and.lb.eq.5) then
c         write(*,*) mype, time,'amr_1blk out',unk1(1,11,5,1,idest)
      endif

      return
      end subroutine amr_1blk_bcset
