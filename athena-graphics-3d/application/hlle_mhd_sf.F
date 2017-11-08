#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#undef LAPIDIUS
       subroutine hlle_mhd(n1,n2,vml,vmr,uc,bn,bf,divv,flux,
     &                     vfout,idir,lb)
       use physcons
       use paramesh_dimensions
       use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */
       implicit none
c       integer, parameter :: nmax  =max(iu_bnd1,ju_bnd1,ku_bnd1)
       integer                    ,intent(in   )::n1,n2,idir,lb
       real, dimension(  8,nmax+1),intent(in   )::vml,vmr
       real, dimension(  8,nmax  ),intent(in   )::uc
       real, dimension(    nmax  ),intent(in   )::bn,divv 
       real, dimension( 11,nmax  ),intent(  out)::flux
       real, dimension(  3,nmax+1),intent(in   )::bf
       real, dimension(  8,nmax  ),intent(  out)::vfout


       real, dimension(10,nmax  )  :: fl
       real, dimension( 7,nmax  )  :: vface,vbar,vl,vr


       integer ::  i,neq
       real fr1,fr2,fr3,fr4,fr5,fr6,fr7,fr8,fr9,fr10
       real fl1,fl2,fl3,fl4,fl5,fl6,fl7,fl8,fl9,fl10
       real rhr,rri,vxr,vyr,vzr,pr,bxr,byr,bzr,etr
       real rhl,rli,vxl,vyl,vzl,pl,bxl,byl,bzl,etl
       real a2l,a2r,cax2r,cax2l,ca2r,ca2l,cfastr,cfastl
       real diff,ekl,ekr,vmaxr,vmaxl,cmax2
       real v2r,v2l,hl,hr,hh,rdum,vh2,bh2
       real bx2,b2r,b2l,b2,sgnr,sgnl,ccl,ccr
       real bv,rvxr,rvyr,rvzr,rvxl,rvyl,rvzl,rbi
       real dif,difi
       real evminl,lamlm,evminr,lamlp
       real rh,uh,vh,wh,buh,bvh,bwh,ph
       real bu2,rhi,a2h,ca2h,cax2h,cfasth,evmin,evmax
       real arg
       real rf,mfx,mfy,mfz,ef,bfx,bfy,bfz,bl,br
       real bxs,bys,bzs,bxt,bytl,bztl,bytr,bztr,b2h,bvth,bwth
       real bdotbs,bsq,bdotu,byt,bzt
       real :: vmax,r_p,ar_p,phi,mu_p,xi_p,coef
       real :: vx,vy,vz,pf,vsq
       real :: bx,by,bz,bdotb0l,bdotb0r
       real :: a2f,ca2f,cax2f,cfastf,nu,div1
       vfout = 0.0
       flux  = 0.0
       vbar  = 0.0
       do i=n1,n2
        vl(1:5,i) = vml(1:5,i)
        vr(1:5,i) = vmr(1:5,i)
        vl(6:7,i) = vml(7:8,i)
        vr(6:7,i) = vmr(7:8,i)
       enddo

       do i=n1,n2
          bx       = bn(  i)
          bxs      = bf(1,i)
          bys      = bf(2,i)
          bzs      = bf(3,i)
          bxt      = bxs + bx
          bx2      = bxt*bxt
c
c    left state
c
          rhl      = vl(1,i)
c         if(lb.eq.73.and.idir.eq.3)then
c          write(*,*)'i,rhl,rhr=',i,rhl,vr(1,i)
c          write(*,*)'i,pl,pr=',i,vl(5,i),vr(5,i)
c         endif
          rli      = 1.0/rhl
          vxl      = vl(2,i)
          vyl      = vl(3,i)
          vzl      = vl(4,i)
          pl       = max(vl(5,i),pmin)
          byl      = vl(6,i)
          bzl      = vl(7,i)
          bytl     = bys+byl
          bztl     = bzs+bzl
          bdotbs   = bx*bxs + byl*bys + bzl*bzs
          bdotu    = bx*vxl + byl*vyl + bzl*vzl
          b2       = bx*bx+byl*byl+bzl*bzl
          b2l      = bx2+bytl*bytl+bztl*bztl
          a2l      = gamma*pl*rli
          rbi      = 1.0/rhl
c          rbi      = 1.0/(rhl+afact*b2l)
          ca2l     = pie4i*b2l*rbi
          cax2l    = pie4i*bx2*rbi
          cmax2    = a2l + ca2l
          arg      = max(0.0,cmax2*cmax2-4.0*a2l*cax2l)
          cfastl   = (0.5*(cmax2+sqrt(arg)))
          cfastl   = sqrt(cfastl)
          v2l      = vxl*vxl+vyl*vyl+vzl*vzl
          vmaxl    = sqrt(v2l)+sqrt(cmax2)
          ekl      = 0.5*rhl*v2l
          etl      = ekl+pl*gammam1i+b2*pie8i
          hl       = (ekl+gamma*pl*gammam1i+b2*pie4i)*rli
c
          bv       = bxt*vxl+bytl*vyl+bztl*vzl
          rvxl     = rhl*vxl
          rvyl     = rhl*vyl
          rvzl     = rhl*vzl
c          rvxl     = (rvxl + afact*(b2l*vxl-bv*bxt ))
c          rvyl     = (rvyl + afact*(b2l*vyl-bv*bytl))
c          rvzl     = (rvzl + afact*(b2l*vzl-bv*bztl))
c
          fl1      = rhl*vxl
          fl2      = fl1*vxl+pl
          fl3      = fl1*vyl 
          fl4      = fl1*vzl 
          fl5      = vxl*(etl+pl)
          fl5      = vxl*(etl+pl+pie8i*b2+pie4i*bdotbs)
     $             -pie4i*bdotu*bxt           

          fl6      = vxl*bytl-bxt*vyl
          fl7      = vxl*bztl-bxt*vzl
#ifdef LORENTZ_FORCE
          fl8      = pie8i*(byl*byl+bzl*bzl-bx*bx)
     $              +pie4i*(bdotbs-2.0*bxs*bx)
          fl9      =-pie4i*(bxt*byl+bx*bys)
          fl10     =-pie4i*(bxt*bzl+bx*bzs)
#else
          fl8      = 0.0
          fl9      = 0.0
          fl10     = 0.0
#endif
          evminl   = vxl + cfastl
          lamlm    = vxl - cfastl
c
c      right state
c
          rhr      = vr(1,i)
          rri      = 1.0/rhr
          vxr      = vr(2,i)
          vyr      = vr(3,i)
          vzr      = vr(4,i)
          pr       = max(vr(5,i),pmin)
          byr      = vr(6,i)
          bzr      = vr(7,i)
          bytr     = bys+byr
          bztr     = bzs+bzr
          bdotbs   = bx*bxs + byr*bys+bzr*bzs
          bdotu    = bx*vxr + byr*vyr+bzr*vzr
          b2       = bx*bx+byr*byr+bzr*bzr
          b2r      = bx2+bytr*bytr+bztr*bztr
          rbi      = 1.0/rhr
c          rbi      = 1.0/(rhr+afact*b2r)
          a2r      = gamma*pr*rri
          ca2r     = pie4i*b2r*rbi
          cax2r    = pie4i*bx2*rbi
          cmax2    = a2r + ca2r
          arg      = max(0.0,cmax2*cmax2-4.0*a2r*cax2r)
          cfastr   = (0.5*(cmax2+sqrt(arg)))
          cfastr   = sqrt(cfastr)
          v2r      = vxr*vxr+vyr*vyr+vzr*vzr
          vmaxr    = sqrt(v2r)+sqrt(cmax2)
          ekr      = 0.5*rhr*v2r
          etr      = ekr+pr*gammam1i+b2*pie8i
          hr       = (ekr+gamma*pr*gammam1i+b2*pie4i)*rri
!
          bv       = bxt*vxr+bytr*vyr+bztr*vzr
          rvxr     = rhr*vxr 
          rvyr     = rhr*vyr 
          rvzr     = rhr*vzr 
!
          fr1      = rhr*vxr
          fr2      = fr1*vxr+pr
          fr3      = fr1*vyr
          fr4      = fr1*vzr
          fr5      = vxr*(etr+pr)
          fr5      = vxr*(etr+pr+pie8i*b2+pie4i*bdotbs)
     $             -pie4i*bdotu*bxt           
          fr6      = vxr*bytr-bxt*vyr
          fr7      = vxr*bztr-bxt*vzr
#ifdef LORENTZ_FORCE
          fr8      = pie8i*(byr*byr+bzr*bzr-bx*bx)
     $              +pie4i*(bdotbs-2.0*bxs*bx)
          fr9      =-pie4i*(bxt*byr+bx*bys)
          fr10     =-pie4i*(bxt*bzr+bx*bzs)
#else
          fr8      = 0.0
          fr9      = 0.0
          fr10     = 0.0
#endif
          evminr   = vxr - cfastr
          lamlp    = vxr + cfastr
!
!
        bdotb0l     = bn(i)*bf(1,i)+vl(6,i)*bf(2,i)+vl(7,i)*bf(3,i)
        ekl         = 0.5*vl(1,i)*(vl(2,i)**2+vl(3,i)**2+vl(4,i)**2)
        b2l         = (bn(i)**2+vl(6,i)**2+vl(7,i)**2)*pie8i
        etl         = ekl+b2l+vl(5,i)*gammam1i
        rli         = 1.0/vl(1,i)

        bdotb0r     = bn(i)*bf(1,i)+vr(6,i)*bf(2,i)+vr(7,i)*bf(3,i)
        ekr         = 0.5*vr(1,i)*(vr(2,i)**2+vr(3,i)**2+vr(4,i)**2)
        b2r         = (bn(i)**2+vr(6,i)**2+vr(7,i)**2)*pie8i
        etr         = ekr+b2r+vr(5,i)*gammam1i
        rri         = 1.0/vr(1,i)

#ifdef ROE_AVERAGE
!
!-------Original Roe Averaging !!
!
          sqrhl     = sqrt(vl(1,i))
          sqrhr     = sqrt(vr(1,i))
          rh        = sqrhr*sqrhl
          rhi         = 1.0/(sqrhr+sqrhl)
          uh          = (sqrhl*vl(2,i)+sqrhr*vr(2,i))*rhi
          vh          = (sqrhl*vl(3,i)+sqrhr*vr(3,i))*rhi
          wh          = (sqrhl*vl(4,i)+sqrhr*vr(4,i))*rhi
          vh2         = (uh**2+vh**2+wh**2)
          hr          = (etr+vr(5,i)+pie8i*(b2r+2.0*bdotb0r))*rri
          hl          = (etl+vl(5,i)+pie8i*(b2l+2.0*bdotb0l))*rli
          hh          = (sqrhl*hl+sqrhr*hr)*rhi
          ph          = (sqrhl*vl(5,i)+sqrhr*vr(5,i))*rhi
          bvh         = (sqrhl*vl(6,i)+sqrhr*vr(6,i))*rhi
          bwh         = (sqrhl*vl(7,i)+sqrhr*vr(7,i))*rhi
          bh2         = bn(i)**2+bvh**2+bwh**2
          bh2         = bh2+2.0*(bf(1,i)*bn(i)+bf(2,i)*bvh+bf(3,i)*bwh)
          a2h         = gammam1*(hh-0.5*vh2+pie4i*bh2/rh)
          ph          = a2h*rh/gamma
          bvh         = bvh+bf(2,i)
          bwh         = bwh+bf(3,i)
          ph       = rh*a2h/gamma
#else
          rh       = 0.5*(rhl  + rhr)
          rhi      = 1.0/rh
          uh       = 0.5*(vxl  + vxr)
          vh       = 0.5*(vyl  + vyr)
          wh       = 0.5*(vzl  + vzr)
          buh      = bx
          bvh      = 0.5*(byl  + byr)
          bwh      = 0.5*(bzl  + bzr)
          ph       = 0.5*(pl   + pr  )
          a2h      = gamma*ph*rhi
#endif
          vbar(1,i)= rh
          vbar(2,i)= uh
          vbar(3,i)= vh
          vbar(4,i)= wh
          vbar(5,i)= ph
          vbar(6,i)= bvh
          vbar(7,i)= bwh
          bvth     = bys+bvh
          bwth     = bzs+bwh
          b2h      = bx2+bvth*bvth+bwth*bwth
          rbi      = 1.0/rh
c          rbi      = 1.0/(rh+afact*b2h)
          ca2h     = pie4i*b2h*rbi
          cax2h    = pie4i*bx2*rbi
          cmax2    = a2h + ca2h
          cfasth   = (0.5*(cmax2+
     $         sqrt(max(0.,cmax2*cmax2-4.0*a2h*cax2h))))
          cfasth   = sqrt(cfasth)
          evmin    = uh - cfasth
          evmax    = uh + cfasth
          evmin    = min(evmin,lamlm,0.)
          evmax    = max(evmax,lamlp,0.)
!
!
!         conserved quantites at the interface
!
          difi = 1.0 /(evmax - evmin)
          rf   =  ((evmax*rhr  - evmin*rhl )-(fr1-fl1))*difi
          mfx  =  ((evmax*rvxr - evmin*rvxl)-(fr2-fl2))*difi
          mfy  =  ((evmax*rvyr - evmin*rvyl)-(fr3-fl3))*difi
          mfz  =  ((evmax*rvzr - evmin*rvzl)-(fr4-fl4))*difi
          ef   =  ((evmax*etr  - evmin*etl )-(fr5-fl5 ))*difi
          bfx  =  bx
          bfy  =  ((evmax*byr  - evmin*byl )-(fr6-fl6))*difi
          bfz  =  ((evmax*bzr  - evmin*bzl )-(fr7-fl7))*difi
!
          dif     = evmin*evmax
!
         fl( 1,i) = difi*((evmax*fl1-evmin*fr1)+dif*(rhr -rhl ))
         fl( 2,i) = difi*((evmax*fl2-evmin*fr2)+dif*(rvxr-rvxl))
         fl( 3,i) = difi*((evmax*fl3-evmin*fr3)+dif*(rvyr-rvyl))
         fl( 4,i) = difi*((evmax*fl4-evmin*fr4)+dif*(rvzr-rvzl))
         fl( 5,i) = difi*((evmax*fl5-evmin*fr5)+dif*(etr -etl ))
         fl( 6,i) = difi*((evmax*fl6-evmin*fr6)+dif*(byr -byl ))
         fl( 7,i) = difi*((evmax*fl7-evmin*fr7)+dif*(bzr -bzl ))

         fl( 8,i) = difi*(evmax*fl8 -evmin*fr8 )
         fl( 9,i) = difi*(evmax*fl9 -evmin*fr9 )
         fl(10,i) = difi*(evmax*fl10-evmin*fr10) 

!=====================================================================
!
          bxt    = bf(1,i)  + bfx
          byt    = bf(2,i)  + bfy
          bzt    = bf(3,i)  + bfz
          b2     = bxt*bxt+byt*byt+bzt*bzt
          bsq    = bfx*bfx+bfy*bfy+bfz*bfz
          rbi    = 1.0/rf
c          rbi    = 1.0/(rf+afact*b2)
          bv     = (bxt*mfx+byt*mfy+bzt*mfz)*rbi
          vx     = mfx*rbi
          vy     = mfy*rbi
          vz     = mfz*rbi
c          vx     = (mfx + afact*bv*bxt)*rbi
c          vy     = (mfy + afact*bv*byt)*rbi
c          vz     = (mfz + afact*bv*bzt)*rbi
          vsq    = vx*vx+vy*vy+vz*vz 
          pf     = gammam1*(ef-0.5*rf*vsq)
          a2f    = gamma*pf/rf
          ca2f   = pie4i*b2*rbi
          cax2f  = pie4i*bxt*bxt*rbi
          cmax2  = a2f + ca2f
          cfastf = (0.5*(cmax2+
     $         sqrt(max(0.,cmax2*cmax2-4.0*a2f*cax2f))))
          cfastf   = sqrt(cfastf)
!
          vface(1,i) = rf
          vface(2,i) = vx
          vface(3,i) = vy
          vface(4,i) = vz
          vface(5,i) = pf
          vface(6,i) = bfy
          vface(7,i) = bfz
      enddo

#ifdef PRIMITIVE_FLUXESX

      do i=n1,n2
       bdotb0   = bn(i)*bf(1,i)+vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)
       bdotul   = bn(i)*vf(2,i)+vf(6,i)*vf(3,i)+vf(7,i)*vf(4,i)
       b2f      = (vf(6,i)**2+vf(7,i)**2-bn(i)**2)*pie8i
       ek       = 0.5*vf(1,i)*(vf(2,i)**2+vf(3,i)**2+vf(4,i)**2)
       eb       = (vf(6,i)**2+vf(7,i)**2+bn(i)**2)
       eb       = eb*pie8i
       et       = ek+eb+vf(5,i)*gammam1i
       bxt      = bf(1,i)  + bn(i)
       byt      = bf(2,i)  + vf(6,i)
       bzt      = bf(3,i)  + vf(7,i)
       fl(1,i)  = vf(1,i)*vf(2,i)
       fl(2,i)  = fl(1,i)*vf(2,i) + vf(5,i) 
       fl(3,i)  = fl(1,i)*vf(3,i) 
       fl(4,i)  = fl(1,i)*vf(4,i) 
       fl(5,i)  = (et+vf(5,i)+eb+pie4i*bdotb0)*vf(2,i)
     &               -bdotul*bxt*pie4i
       fl(6,i)  = vf(2,i)*byt-vf(3,i)*bxt
       fl(7,i)  = vf(2,i)*bzt-vf(4,i)*bxt

       fl( 8,i) = b2f + pie4i*(bdotb0-2.0*bf(1,i)*bn(i))
       fl( 9,i) =-(bxt*vf(6,i)+bn(i)*bf(2,i))*pie4i
       fl(10,i) =-(bxt*vf(7,i)+bn(i)*bf(3,i))*pie4i
      enddo
#endif PRIMITIVE_FLUXESX

#ifdef LAPIDIUS
!
!-----add artificial viscosity
!
!-----Lapidius AV
      nu = 0.1
      do i=n1,n2
       div1     = nu*min( divv(i),0.)
       fl(:5,i) = fl(:5,i) + div1*(uc(:5,i)-uc(:5,i-1))
       fl( 6,i) = fl( 6,i) + div1*(uc( 7,i)-uc( 7,i-1))
       fl( 7,i) = fl( 7,i) + div1*(uc( 8,i)-uc( 8,i-1))
      enddo
#endif /*LAPIDIUS*/
#ifdef CFD
       fl(6:10,n1:n2)=0.0
#endif /*CFD*/
!
!-----rotate fluxes here
!
      if(idir.eq.1)then

       do i = n1,n2
        flux( 1,i) = fl( 1,i)
        flux( 2,i) = fl( 2,i) 
        flux( 3,i) = fl( 3,i)
        flux( 4,i) = fl( 4,i)
        flux( 5,i) = fl( 5,i)
        flux( 6,i) = 0.0
        flux( 7,i) = fl( 6,i) 
        flux( 8,i) = fl( 7,i)
        flux( 9,i) = fl( 8,i) 
        flux(10,i) = fl( 9,i)
        flux(11,i) = fl(10,i)
!
!------primitives at faces
!
        vfout( :5,i) = vface( :5,i)
        vfout(  6,i) = bn(i)
        vfout(7:8,i) = vface(6:7,i)
!
       enddo

      elseif(idir.eq.2)then

       do i = n1,n2

        flux( 1,i) = fl(1,i)
        flux( 3,i) = fl(2,i) 
        flux( 4,i) = fl(3,i)
        flux( 2,i) = fl(4,i)
        flux( 5,i) = fl(5,i)
        flux( 7,i) = 0.0 
        flux( 8,i) = fl( 6,i)
        flux( 6,i) = fl( 7,i) 
        flux(10,i) = fl( 8,i) 
        flux(11,i) = fl( 9,i)
        flux( 9,i) = fl(10,i)
!
!------primitives at faces
!
        vfout(1,i) = vface(1,i)
        vfout(3,i) = vface(2,i)
        vfout(4,i) = vface(3,i)
        vfout(2,i) = vface(4,i)
        vfout(5,i) = vface(5,i)
        vfout(7,i) = bn(i)
        vfout(8,i) = vface(6,i)
        vfout(6,i) = vface(7,i)
!
       enddo

      else

       do i = n1,n2

        flux( 1,i) = fl( 1,i)
        flux( 4,i) = fl( 2,i)
        flux( 2,i) = fl( 3,i)
        flux( 3,i) = fl( 4,i)
        flux( 5,i) = fl( 5,i)
        flux( 8,i) = 0.0
        flux( 6,i) = fl( 6,i) 
        flux( 7,i) = fl( 7,i) 
        flux(11,i) = fl( 8,i) 
        flux( 9,i) = fl( 9,i)
        flux(10,i) = fl(10,i)
!
!------primitives at faces
!
        vfout(1,i) = vface(1,i)
        vfout(4,i) = vface(2,i)
        vfout(2,i) = vface(3,i)
        vfout(3,i) = vface(4,i)
        vfout(5,i) = vface(5,i)
        vfout(8,i) = bn(i)
        vfout(6,i) = vface(6,i)
        vfout(7,i) = vface(7,i)
       enddo


      endif
c

       return 
       end
