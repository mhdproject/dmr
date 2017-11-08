#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#undef PRIMITIVE_FLUXES
#undef HARTEN_ENROPY_FIX
#undef  STRONG_RAREFACTION_WAVE_FIX
#undef HLLEM_ENTROPY_FIX
#undef ROE_AVERAGE
      subroutine eigenvalues(n1,n2,vv,bx,eigenval,cfast)
      use paramesh_dimensions
      use physcons
      integer, parameter :: nmax  =max(iu_bnd1,ju_bnd1,ku_bnd1)
      integer i,j,k,n1,n2
      real vv(7,nmax),bx(nmax),cfast(nmax)
      real eigenval(7,nmax)
      do i=n1,n2
        rho    = max(rmin,vv(1,i))
        rhoi   = 1.0/rho
        p      = vv(5,i)
        u      = vv(2,i)
        v      = vv(3,i)
        w      = vv(4,i)
        by     = vv(6,i)
        bz     = vv(7,i)
        q2     = (u**2+v**2+w**2)
        bperp  = by**2+bz**2
        bx2    = bx(i)*bx(i)
        b2     = (bx2+bperp)
        p      = max(pmin,p)
        a2     = gamma*p*rhoi
        a2i    = 1.0/a2
        a      = sqrt(a2)
        ca2    = pie4i*b2*rhoi
        cax2   = pie4i*bx2*rhoi
        cax    = sqrt(cax2)
        cmax2  = a2 + ca2
        arg    = max(0.0,cmax2*cmax2-4.0*a2*cax2)
        arg    = sqrt(arg)
        cf22   = 0.5*(cmax2+arg)
        cf     = sqrt(cf22)
        cs2    = 0.5*(cmax2-arg)
        cs     = sqrt(max(cs2,0.0))
c
        cfast(i) = cf
c
c----compute eigenvalues
c
        eigenval(1,i) = u-cf
        eigenval(2,i) = u-cax
        eigenval(3,i) = u-cs
        eigenval(4,i) = u
        eigenval(5,i) = u+cs
        eigenval(6,i) = u+cax
        eigenval(7,i) = u+cf

      enddo
      return
      end
c
c--------------------------------------------------------------------------
c
      subroutine mhd_eigenvectors ( rhe ,lhe ,rhevc,lhevc,eigenval,
     &                              vv  ,bx  ,wt   ,dv   ,
     &                              ivec,n1,n2,cfast,speeds,dvy,dvz,
     &                              idir)
      use paramesh_dimensions
      use physcons
      integer :: ivec,n1,n2,idir
      integer, parameter :: nmax  =max(iu_bnd1,ju_bnd1,ku_bnd1)
      real,dimension(  4,nmax) :: speeds
      real,dimension(  7,nmax) :: dv,vv,wt,eigenval
      real,dimension(7,7,nmax) :: rhe  ,lhe  
      real,dimension(7,7,nmax) :: rhevc,lhevc
      real,dimension(7,7,nmax) :: M    ,MI
      real,dimension(    nmax) :: bx,cfast,dvy,dvz
      real, parameter :: small_speed = 0.001
      real            :: sqrt2i 
      sqrt2i = 1.0/sqrt(2.0)


      if(ivec.eq.1)then
c
c-----compute primitive eigenvectors
c

      do i=n1,n2
        rho    = max(rmin,vv(1,i))
        rhoi   = 1.0/rho
        p      = vv(5,i)
        u      = vv(2,i)
        v      = vv(3,i)
        w      = vv(4,i)
        by     = vv(6,i)
        bz     = vv(7,i)
        q2     = (u**2+v**2+w**2)
        bperp  = by**2+bz**2
        bx2    = bx(i)*bx(i)
        b2     = (bx2+bperp)
        p      = max(pmin,p)
        aa2    = gamma*p*rhoi
        a2i    = 1.0/aa2
        a      = sqrt(aa2)
        ca2    = pie4i*b2*rhoi
        cax2   = pie4i*bx2*rhoi
        cax    = sqrt(cax2)
        cmax2  = aa2 + ca2
        
        arg    = max(0.0,cmax2*cmax2-4.0*aa2*cax2)
        arg    = sqrt(arg)
        cf22   = 0.5*(cmax2+arg)
        cf     = sqrt(cf22)
        cs2    = 0.5*(cmax2-arg)
        cs     = sqrt(max(cs2,0.0))
        sgnBx  = sign(1.0,bx(i)) 
c 
      
        bxd  = bx(i)
        byd  = by
        bzd  = bz

        s4prho = sqrt(pie4*rho)
        s4prhoi= 1.0/s4prho
        vay    = by*s4prho
        vaz    = bz*s4prho
        vaperp = sqrt(bperp)*s4prhoi
        eps    = small_speed*a
        if(vaperp.gt.eps)then
          vpi   = 1.0/(vaperp*s4prho)
          betay = by*vpi
          betaz = bz*vpi
        else 
           vt=sqrt(dvy(i)*dvy(i)+dvz(i)*dvz(i))
           if(.true. .and. vt.gt.eps) then
              betay=sgnBx*dvy(i)/vt
              betaz=sgnBx*dvz(i)/vt
           else
              betay = 0.
              betaz = 0.
           endif
        endif
c
        
        if(cf.ge.(cs+eps))then
c
c-----far from umbilic point
c
          tempaa = 1.0/(cf22-cs2)
          alphaf = max((aa2 -cs2),0.0)
          alphaf = sqrt(tempaa*alphaf)
          alphas = max((cf22 -aa2),0.0)
          alphas = sqrt(tempaa*alphas)
        elseif(abs(cax-a)>smin*a)then
c
c-----near umbilic point
c
          tempaa = vaperp/(cax - a)
          tempaa = 0.5*atan(tempaa)
          alphaf = abs(sin(tempaa))
          alphas = abs(cos(tempaa))
        else
c
c-----at umbilic point
c
          tempaa = 0.25*pie*sign(1.0,(cax - a))
          alphaf = abs(sin(tempaa))
          alphas = abs(cos(tempaa))
        endif
c
c----compute eigenvalues
c
        eigenval(1,i) = u-cf
        eigenval(2,i) = u-cax
        eigenval(3,i) = u-cs
        eigenval(4,i) = u
        eigenval(5,i) = u+cs
        eigenval(6,i) = u+cax
        eigenval(7,i) = u+cf
c
        speeds(1,i)   = a
        speeds(2,i)   = cs
        speeds(3,i)   = cax
        speeds(4,i)   = cf
        cfast(i)      = cf
c
c----compute right hand primitive eigenvectors
c
c
c----compute right hand primitive eigenvectors
c
        coef0      = s4prho*sgnBx
        coef1      = rho*alphas
        coef2      = rho*alphaf
        coef3      = a*s4prhoi*alphas
        coef4      = 0.5*a2i
        coef5      = a*s4prhoi*alphaf
        coef6      = alphas*cs*sgnBx
        coef7      = a*s4prho*alphas
        coef8      = a*s4prho*alphaf
        coef9      = alphaf*cf*sgnBx
        coefa      = 0.5*s4prhoi*sgnBx
        coefb      = pie4i*rhoi*coef4
        coefc      = coef4*alphas*a*s4prhoi
        coefd      = coef4*alphaf*a*s4prhoi
c
c------fast modes
c
        a1         = coef4*alphaf*cf*dv(2,i)
        a2         = coef4*alphaf*dv(5,i)*rhoi
        a3         = coef4*coef3*(betay*dv(6,i)+betaz*dv(7,i))
        a4         = coef4*coef6*(betay*dv(3,i)+betaz*dv(4,i))
        wt(1,i)    =-a1+a2+a3+a4
        wt(7,i)    = a1+a2+a3-a4
c
c------Alfven modes
c
        a1         = 0.5*betaz*dv(3,i)
        a2         = 0.5*betay*dv(4,i)
        a3         = betaz*coefa*dv(6,i)
        a4         = betay*coefa*dv(7,i)
        wt(2,i)    =-a1+a2-a3+a4
        wt(6,i)    = a1-a2-a3+a4
c
c-----Slow modes
c
        a1         = coef4*alphas*cs*dv(2,i)
        a2         = coef4*alphas*dv(5,i)*rhoi
        a3         = coef4*coef5*(betay*dv(6,i)+betaz*dv(7,i))
        a4         = coef4*coef9*(betay*dv(3,i)+betaz*dv(4,i))
        wt(3,i)    =-a1+a2-a3-a4
        wt(5,i)    = a1+a2-a3+a4
c
c-----Entropy wave
c
        wt(4,i)    = dv(1,i)-a2i*dv(5,i)

c
c------Alfven modes
c
        rhe(1,2,i) = 0.0
        rhe(1,6,i) = 0.0

        rhe(2,2,i) = 0.0
        rhe(2,6,i) = 0.0

        rhe(3,2,i) =-betaz
        rhe(3,6,i) = betaz

        rhe(4,2,i) = betay
        rhe(4,6,i) =-betay

        rhe(5,2,i) = 0.0
        rhe(5,6,i) = 0.0

        rhe(6,2,i) =-s4prho*sgnBx*betaz
        rhe(6,6,i) = rhe(6,2,i)

        rhe(7,2,i) = s4prho*sgnBx*betay
        rhe(7,6,i) = rhe(7,2,i)
c
        lhe(2,1,i) = 0.0
        lhe(6,1,i) = 0.0

        lhe(2,2,i) = 0.0
        lhe(6,2,i) = 0.0

        lhe(2,3,i) =-0.5*betaz
        lhe(6,3,i) = 0.5*betaz

        lhe(2,4,i) = 0.5*betay
        lhe(6,4,i) =-0.5*betay

        lhe(2,5,i) = 0.0
        lhe(6,5,i) = 0.0

        lhe(2,6,i) =-0.5*s4prhoi*sgnBx*betaz
        lhe(6,6,i) = lhe(2,6,i)

        lhe(2,7,i) = 0.5*s4prhoi*sgnBx*betay
        lhe(6,7,i) = lhe(2,7,i)
c
c------fast modes
c
        rhe(1,1,i) = rho*alphaf
        rhe(1,7,i) = rhe(1,1,i)

        rhe(2,1,i) =-alphaf*cf
        rhe(2,7,i) =-rhe(2,1,i)

        rhe(3,1,i) = alphas*cs*sgnBx*betay
        rhe(3,7,i) =-rhe(3,1,i)

        rhe(4,1,i) = alphas*cs*sgnBx*betaz
        rhe(4,7,i) =-rhe(4,1,i)

        rhe(5,1,i) = rho*alphaf*aa2 
        rhe(5,7,i) = rhe(5,1,i)

        rhe(6,1,i) = a*s4prho*alphas*betay
        rhe(6,7,i) = rhe(6,1,i)

        rhe(7,1,i) = a*s4prho*alphas*betaz
        rhe(7,7,i) = rhe(7,1,i)
c
        lhe(1,1,i) = 0.0
        lhe(7,1,i) = 0.0

        lhe(1,2,i) =-0.5*a2i*alphaf*cf
        lhe(7,2,i) =-lhe(1,2,i)

        lhe(1,3,i) = 0.5*a2i*alphas*cs*sgnBx*betay
        lhe(7,3,i) =-lhe(1,3,i)

        lhe(1,4,i) = 0.5*a2i*alphas*cs*sgnBx*betaz
        lhe(7,4,i) =-lhe(1,4,i)

        lhe(1,5,i) = 0.5*a2i*alphaf*rhoi 
        lhe(7,5,i) = lhe(1,5,i)

        lhe(1,6,i) = 0.5*a2i*a*s4prhoi*alphas*betay  ! <<==
        lhe(7,6,i) = lhe(1,6,i)

        lhe(1,7,i) = 0.5*a2i*a*s4prhoi*alphas*betaz
        lhe(7,7,i) = lhe(1,7,i)

c
c-----Slow mode
c
        rhe(1,3,i) = rho*alphas
        rhe(1,5,i) = rhe(1,3,i)

        rhe(2,3,i) =-alphas*cs
        rhe(2,5,i) =-rhe(2,3,i)

        rhe(3,3,i) =-alphaf*cf*sgnBx*betay
        rhe(3,5,i) =-rhe(3,3,i)

        rhe(4,3,i) =-alphaf*cf*sgnBx*betaz
        rhe(4,5,i) =-rhe(4,3,i)

        rhe(5,3,i) = rho*alphas*aa2
        rhe(5,5,i) = rhe(5,3,i)

        rhe(6,3,i) =-alphaf*a*s4prho*betay
        rhe(6,5,i) = rhe(6,3,i)

        rhe(7,3,i) =-alphaf*a*s4prho*betaz
        rhe(7,5,i) = rhe(7,3,i)
c
        lhe(3,1,i) = 0.0
        lhe(5,1,i) = 0.0

        lhe(3,2,i) =-0.5*a2i*alphas*cs
        lhe(5,2,i) =-lhe(3,2,i)

        lhe(3,3,i) =-0.5*a2i*alphaf*cf*sgnBx*betay
        lhe(5,3,i) =-lhe(3,3,i)

        lhe(3,4,i) =-0.5*a2i*alphaf*cf*sgnBx*betaz
        lhe(5,4,i) =-lhe(3,4,i)

        lhe(3,5,i) = 0.5*a2i*alphas*rhoi
        lhe(5,5,i) = lhe(3,5,i)

        lhe(3,6,i) =-0.5*a2i*alphaf*a*s4prhoi*betay
        lhe(5,6,i) = lhe(3,6,i)

        lhe(3,7,i) =-0.5*a2i*alphaf*a*s4prhoi*betaz
        lhe(5,7,i) = lhe(3,7,i)
c
c-----Entropy wave
c
        rhe(1,4,i) = 1.0
        rhe(2,4,i) = 0.0
        rhe(3,4,i) = 0.0
        rhe(4,4,i) = 0.0
        rhe(5,4,i) = 0.0
        rhe(6,4,i) = 0.0
        rhe(7,4,i) = 0.0
c
        lhe(4,1,i) = 1.0
        lhe(4,2,i) = 0.0
        lhe(4,3,i) = 0.0
        lhe(4,4,i) = 0.0
        lhe(4,5,i) =-a2i
        lhe(4,6,i) = 0.0
        lhe(4,7,i) = 0.0
c
      enddo
      return

      else
c
c-----compute conservative eigenvectors
c
      do i=n1,n2
        rho    = max(rmin,vv(1,i))
        rhoi   = 1.0/rho
        p      = vv(5,i)
        u      = vv(2,i)
        v      = vv(3,i)
        w      = vv(4,i)
        by     = vv(6,i)
        bz     = vv(7,i)
        q2     = (u**2+v**2+w**2)
        bperp  = by**2+bz**2
        bx2    = bx(i)*bx(i)
        b2     = (bx2+bperp)
        p      = max(pmin,p)
        aa2    = gamma*p*rhoi
        a2i    = 1.0/aa2
        a      = sqrt(aa2)
        ca2    = pie4i*b2*rhoi
        cax2   = pie4i*bx2*rhoi
        cax    = sqrt(cax2)
        cmax2  = aa2 + ca2
        
        arg    = max(0.0,cmax2*cmax2-4.0*aa2*cax2)
        arg    = sqrt(arg)
        cf22   = 0.5*(cmax2+arg)
        cf     = sqrt(cf22)
        cs2    = 0.5*(cmax2-arg)
        cs     = sqrt(max(cs2,0.0))
        sgnBx  = sign(1.0,bx(i)) 
c 
      
        bxd  = bx(i)
        byd  = by
        bzd  = bz

        s4prho = sqrt(pie4*rho)
        s4prhoi= 1.0/s4prho
        vay    = by*s4prho
        vaz    = bz*s4prho
        vaperp = sqrt(bperp)*s4prhoi
        eps    = small_speed*a
        if(vaperp.gt.eps)then
          vpi   = 1.0/(vaperp*s4prho)
          betay = by*vpi
          betaz = bz*vpi
        else 
          vt=sqrt(dvy(i)*dvy(i)+dvz(i)*dvz(i))
          if(.true. .and. vt.gt.eps) then
             betay=sgnBx*dvy(i)/vt
             betaz=sgnBx*dvz(i)/vt
          else
             betay = 0.
             betaz = 0.
           endif
        endif
c
        if(cf.ge.(cs+eps))then
c
c-----far from umbilic point
c
          tempaa = 1.0/(cf22-cs2)
          alphaf = max((aa2 -cs2),0.0)
          alphaf = sqrt(tempaa*alphaf)
          alphas = max((cf22 -aa2),0.0)
          alphas = sqrt(tempaa*alphas)
        elseif(abs(cax-a)>smin*a)then
c
c-----near umbilic point
c
          tempaa = vaperp/(cax - a)
          tempaa = 0.5*atan(tempaa)
          alphaf = abs(sin(tempaa))
          alphas = abs(cos(tempaa))
        else
c
c-----at umbilic point
c
          tempaa = 0.25*pie*sign(1.0,(cax - a))
          alphaf = abs(sin(tempaa))
          alphas = abs(cos(tempaa))
        endif
c
c----compute eigenvalues
c
        eigenval(1,i) = u-cf
        eigenval(2,i) = u-cax
        eigenval(3,i) = u-cs
        eigenval(4,i) = u
        eigenval(5,i) = u+cs
        eigenval(6,i) = u+cax
        eigenval(7,i) = u+cf
c
        speeds(1,i)   = a
        speeds(2,i)   = cs
        speeds(3,i)   = cax
        speeds(4,i)   = cf
        cfast(i)      = cf
c
c----compute right hand primitive eigenvectors
c
c
c----compute right hand primitive eigenvectors
c
        coef0      = s4prho*sgnBx
        coef1      = rho*alphas
        coef2      = rho*alphaf
        coef3      = a*s4prhoi*alphas
        coef4      = 0.5*a2i
        coef5      = a*s4prhoi*alphaf
        coef6      = alphas*cs*sgnBx
        coef7      = a*s4prho*alphas
        coef8      = a*s4prho*alphaf
        coef9      = alphaf*cf*sgnBx
        coefa      = 0.5*s4prhoi*sgnBx
        coefb      = pie4i*rhoi*coef4
        coefc      = coef4*alphas*a*s4prhoi
        coefd      = coef4*alphaf*a*s4prhoi
c
c------fast modes
c
        a1         = coef4*alphaf*cf*dv(2,i)
        a2         = coef4*alphaf*dv(5,i)*rhoi
        a3         = coef4*coef3*(betay*dv(6,i)+betaz*dv(7,i))
        a4         = coef4*coef6*(betay*dv(3,i)+betaz*dv(4,i))
        wt(1,i)    =-a1+a2+a3+a4
        wt(7,i)    = a1+a2+a3-a4
c
c------Alfven modes
c
        a1         = 0.5*betaz*dv(3,i)
        a2         = 0.5*betay*dv(4,i)
        a3         = betaz*coefa*dv(6,i)
        a4         = betay*coefa*dv(7,i)
        wt(2,i)    =-a1+a2-a3+a4
        wt(6,i)    = a1-a2-a3+a4
c
c-----Slow modes
c
        a1         = coef4*alphas*cs*dv(2,i)
        a2         = coef4*alphas*dv(5,i)*rhoi
        a3         = coef4*coef5*(betay*dv(6,i)+betaz*dv(7,i))
        a4         = coef4*coef9*(betay*dv(3,i)+betaz*dv(4,i))
        wt(3,i)    =-a1+a2-a3-a4
        wt(5,i)    = a1+a2-a3+a4
c
c-----Entropy wave
c
        wt(4,i)    = dv(1,i)-a2i*dv(5,i)
c
c------Alfven modes
c
        rhe(1,2,i) = 0.0
        rhe(1,6,i) = 0.0

        rhe(2,2,i) = 0.0
        rhe(2,6,i) = 0.0

        rhe(3,2,i) =-betaz
        rhe(3,6,i) = betaz

        rhe(4,2,i) = betay
        rhe(4,6,i) =-betay

        rhe(5,2,i) = 0.0
        rhe(5,6,i) = 0.0

        rhe(6,2,i) =-s4prho*sgnBx*betaz
        rhe(6,6,i) = rhe(6,2,i)

        rhe(7,2,i) = s4prho*sgnBx*betay
        rhe(7,6,i) = rhe(7,2,i)
c
        lhe(2,1,i) = 0.0
        lhe(6,1,i) = 0.0

        lhe(2,2,i) = 0.0
        lhe(6,2,i) = 0.0

        lhe(2,3,i) =-0.5*betaz
        lhe(6,3,i) = 0.5*betaz

        lhe(2,4,i) = 0.5*betay
        lhe(6,4,i) =-0.5*betay

        lhe(2,5,i) = 0.0
        lhe(6,5,i) = 0.0

        lhe(2,6,i) =-0.5*s4prhoi*sgnBx*betaz
        lhe(6,6,i) = lhe(2,6,i)

        lhe(2,7,i) = 0.5*s4prhoi*sgnBx*betay
        lhe(6,7,i) = lhe(2,7,i)
c
c------fast modes
c
        rhe(1,1,i) = rho*alphaf
        rhe(1,7,i) = rhe(1,1,i)

        rhe(2,1,i) =-alphaf*cf
        rhe(2,7,i) =-rhe(2,1,i)

        rhe(3,1,i) = alphas*cs*sgnBx*betay
        rhe(3,7,i) =-rhe(3,1,i)

        rhe(4,1,i) = alphas*cs*sgnBx*betaz
        rhe(4,7,i) =-rhe(4,1,i)

        rhe(5,1,i) = rho*alphaf*aa2 
        rhe(5,7,i) = rhe(5,1,i)

        rhe(6,1,i) = a*s4prho*alphas*betay
        rhe(6,7,i) = rhe(6,1,i)

        rhe(7,1,i) = a*s4prho*alphas*betaz
        rhe(7,7,i) = rhe(7,1,i)
c
        lhe(1,1,i) = 0.0
        lhe(7,1,i) = 0.0

        lhe(1,2,i) =-0.5*a2i*alphaf*cf
        lhe(7,2,i) =-lhe(1,2,i)

        lhe(1,3,i) = 0.5*a2i*alphas*cs*sgnBx*betay
        lhe(7,3,i) =-lhe(1,3,i)

        lhe(1,4,i) = 0.5*a2i*alphas*cs*sgnBx*betaz
        lhe(7,4,i) =-lhe(1,4,i)

        lhe(1,5,i) = 0.5*a2i*alphaf*rhoi 
        lhe(7,5,i) = lhe(1,5,i)

        lhe(1,6,i) = 0.5*a2i*a*s4prhoi*alphas*betay  ! <<==
        lhe(7,6,i) = lhe(1,6,i)

        lhe(1,7,i) = 0.5*a2i*a*s4prhoi*alphas*betaz
        lhe(7,7,i) = lhe(1,7,i)

c
c-----Slow mode
c
        rhe(1,3,i) = rho*alphas
        rhe(1,5,i) = rhe(1,3,i)

        rhe(2,3,i) =-alphas*cs
        rhe(2,5,i) =-rhe(2,3,i)

        rhe(3,3,i) =-alphaf*cf*sgnBx*betay
        rhe(3,5,i) =-rhe(3,3,i)

        rhe(4,3,i) =-alphaf*cf*sgnBx*betaz
        rhe(4,5,i) =-rhe(4,3,i)

        rhe(5,3,i) = rho*alphas*aa2
        rhe(5,5,i) = rhe(5,3,i)

        rhe(6,3,i) =-alphaf*a*s4prho*betay
        rhe(6,5,i) = rhe(6,3,i)

        rhe(7,3,i) =-alphaf*a*s4prho*betaz
        rhe(7,5,i) = rhe(7,3,i)
c
        lhe(3,1,i) = 0.0
        lhe(5,1,i) = 0.0

        lhe(3,2,i) =-0.5*a2i*alphas*cs
        lhe(5,2,i) =-lhe(3,2,i)

        lhe(3,3,i) =-0.5*a2i*alphaf*cf*sgnBx*betay
        lhe(5,3,i) =-lhe(3,3,i)

        lhe(3,4,i) =-0.5*a2i*alphaf*cf*sgnBx*betaz
        lhe(5,4,i) =-lhe(3,4,i)

        lhe(3,5,i) = 0.5*a2i*alphas*rhoi
        lhe(5,5,i) = lhe(3,5,i)

        lhe(3,6,i) =-0.5*a2i*alphaf*a*s4prhoi*betay
        lhe(5,6,i) = lhe(3,6,i)

        lhe(3,7,i) =-0.5*a2i*alphaf*a*s4prhoi*betaz
        lhe(5,7,i) = lhe(3,7,i)
c
c-----Entropy wave
c
        rhe(1,4,i) = 1.0
        rhe(2,4,i) = 0.0
        rhe(3,4,i) = 0.0
        rhe(4,4,i) = 0.0
        rhe(5,4,i) = 0.0
        rhe(6,4,i) = 0.0
        rhe(7,4,i) = 0.0
c
        lhe(4,1,i) = 1.0
        lhe(4,2,i) = 0.0
        lhe(4,3,i) = 0.0
        lhe(4,4,i) = 0.0
        lhe(4,5,i) =-a2i
        lhe(4,6,i) = 0.0
        lhe(4,7,i) = 0.0
c
c
c-----Transform primitive eigenvectors to conservative eigen vectors
c     See Hirsch Vol. II p.147 eqns.(16.2.47)
c

        M(1,1,i) = 1.0
c       M(1,2,i) = 0.0
c       M(1,3,i) = 0.0
c       M(1,4,i) = 0.0
c       M(1,5,i) = 0.0
c       M(1,6,i) = 0.0
c       M(1,7,i) = 0.0


        M(2,1,i) = u
        M(2,2,i) = rho
c       M(2,3,i) = 0.0
c       M(2,4,i) = 0.0
c       M(2,5,i) = 0.0
c       M(2,6,i) = 0.0
c       M(2,7,i) = 0.0


        M(3,1,i) = v
c       M(3,2,i) = 0.0
        M(3,3,i) = rho
c       M(3,4,i) = 0.0
c       M(3,5,i) = 0.0
c       M(3,6,i) = 0.0
c       M(3,7,i) = 0.0


        M(4,1,i) = w
c       M(4,2,i) = 0.0
c       M(4,3,i) = 0.0
        M(4,4,i) = rho
c       M(4,5,i) = 0.0
c       M(4,6,i) = 0.0
c       M(4,7,i) = 0.0


        M(5,1,i) = 0.5*q2
        M(5,2,i) = rho*u
        M(5,3,i) = rho*v
        M(5,4,i) = rho*w
        M(5,5,i) = gammam1i
        M(5,6,i) = by*pie4i
        M(5,7,i) = bz*pie4i


c       M(6,1,i) = 0.0
c       M(6,2,i) = 0.0
c       M(6,3,i) = 0.0
c       M(6,4,i) = 0.0
c       M(6,5,i) = 0.0
        M(6,6,i) = 1.0
c       M(6,7,i) = 0.0


c       M(7,1,i) = 0.0
c       M(7,2,i) = 0.0
c       M(7,3,i) = 0.0
c       M(7,4,i) = 0.0
c       M(7,5,i) = 0.0
c       M(7,6,i) = 0.0
        M(7,7,i) = 1.0
c

        MI(1,1,i) = 1.0
        MI(2,1,i) =-rhoi*u
        MI(3,1,i) =-rhoi*v
        MI(4,1,i) =-rhoi*w
        MI(5,1,i) = 0.5*gammam1*q2
c       MI(6,1,i) = 0.0
c       MI(7,1,i) = 0.0


c       MI(1,2,i) = 0.0
        MI(2,2,i) = rhoi
c       MI(3,2,i) = 0.0
c       MI(4,2,i) = 0.0
        MI(5,2,i) =-gammam1*u
c       MI(6,2,i) = 0.0
c       MI(7,2,i) = 0.0


c       MI(1,3,i) = 0.0
c       MI(2,3,i) = 0.0
        MI(3,3,i) = rhoi
c       MI(4,3,i) = 0.0
        MI(5,3,i) =-gammam1*v
c       MI(6,3,i) = 0.0
c       MI(7,3,i) = 0.0


c       MI(1,4,i) = 0.0
c       MI(2,4,i) = 0.0
c       MI(3,4,i) = 0.0
        MI(4,4,i) = rhoi
        MI(5,4,i) =-gammam1*w
c       MI(6,4,i) = 0.0
c       MI(7,4,i) = 0.0


c       MI(1,5,i) = 0.0
c       MI(2,5,i) = 0.0
c       MI(3,5,i) = 0.0
c       MI(4,5,i) = 0.0
        MI(5,5,i) = gammam1
c       MI(6,5,i) = 0.0
c       MI(7,5,i) = 0.0


c       MI(1,6,i) = 0.0
c       MI(2,6,i) = 0.0
c       MI(3,6,i) = 0.0
c       MI(4,6,i) = 0.0
        MI(5,6,i) =-gammam1*by*pie4i
        MI(6,6,i) = 1.0
c       MI(7,6,i) = 0.0


c       MI(1,7,i) = 0.0
c       MI(2,7,i) = 0.0
c       MI(3,7,i) = 0.0
c       MI(4,7,i) = 0.0
        MI(5,7,i) =-gammam1*bz*pie4i
c       MI(6,7,i) = 0.0
        MI(7,7,i) = 1.0
      enddo
c
c-----compute conserved eigenvectors
c
      do i=n1,n2
        do k=1,7
           rhevc(1,k,i) = M(1,1,i)*rhe(1,k,i)
           rhevc(2,k,i) = M(2,1,i)*rhe(1,k,i) + M(2,2,i)*rhe(2,k,i)
           rhevc(3,k,i) = M(3,1,i)*rhe(1,k,i) + M(3,3,i)*rhe(3,k,i)
           rhevc(4,k,i) = M(4,1,i)*rhe(1,k,i) + M(4,4,i)*rhe(4,k,i)
           rhevc(5,k,i) = M(5,1,i)*rhe(1,k,i) + M(5,2,i)*rhe(2,k,i)
     &                  + M(5,3,i)*rhe(3,k,i) + M(5,4,i)*rhe(4,k,i)
     &                  + M(5,5,i)*rhe(5,k,i) + M(5,6,i)*rhe(6,k,i)
     &                  + M(5,7,i)*rhe(7,k,i)
           rhevc(6,k,i) = M(6,6,i)*rhe(6,k,i)
           rhevc(7,k,i) = M(7,7,i)*rhe(7,k,i)


           lhevc(k,1,i) = lhe(k,1,i)*MI(1,1,i) + lhe(k,2,i)*MI(2,1,i)
     &                  + lhe(k,3,i)*MI(3,1,i) + lhe(k,4,i)*MI(4,1,i)
     &                  + lhe(k,5,i)*MI(5,1,i)
           lhevc(k,2,i) = lhe(k,2,i)*MI(2,2,i) + lhe(k,5,i)*MI(5,2,i)
           lhevc(k,3,i) = lhe(k,3,i)*MI(3,3,i) + lhe(k,5,i)*MI(5,3,i)
           lhevc(k,4,i) = lhe(k,4,i)*MI(4,4,i) + lhe(k,5,i)*MI(5,4,i)
           lhevc(k,5,i) = lhe(k,5,i)*MI(5,5,i)
           lhevc(k,6,i) = lhe(k,5,i)*MI(5,6,i) + lhe(k,6,i)*MI(6,6,i)
           lhevc(k,7,i) = lhe(k,5,i)*MI(5,7,i) + lhe(k,7,i)*MI(7,7,i)
        enddo   
      enddo
      endif
      return
      end
c=======================================================================
c
c=======================================================================
      subroutine roe_mhd(n1,n2,vml,vmr,uc,bn,bf,divv,flux,
     &                      vfout,idir)
      use paramesh_dimensions
      use physcons
      use grid
      use athena_specific_parameters
c      integer, parameter :: nmax  =max(iu_bnd1,ju_bnd1,ku_bnd1)
      integer                    ,intent(in   )::n1,n2,idir
      real, dimension(  8,nmax+1),intent(in   )::vml,vmr
      real, dimension(  8,nmax  ),intent(in   )::uc
      real, dimension(    nmax  ),intent(in   )::bn,divv 
      real, dimension( 11,nmax  ),intent(  out)::flux
      real, dimension(  3,nmax+1),intent(in   )::bf
      real, dimension(  8,nmax  ),intent(  out)::vfout
      integer  :: i,j,k,l,ientr,ii,jj,kk,ll,kn,km
      real, dimension(  7,nmax  ) :: vl,vr,vdl,vdr,ul,ur,eigenwt,dfl,dul
      real, dimension(  7,nmax  ) :: dvl
      real, dimension(    nmax  ) :: dvy, dvz
      real, dimension(7,7,nmax  ) :: rhe,lhe,rhec,lhec
      real, dimension(  7,nmax  ) :: eigenval,sign_eval,eps_zone
      real, dimension( 10,nmax  ) :: fl,fr
      real, dimension(  7,nmax  ) :: eval_l,eval_r,vbar,delu,delv,vf,vdum
      real, dimension(  7,nmax  ) :: vbar_l,vbar_r,wt,wtt,uf,frd,fld
      real, dimension(    nmax  ) :: cfastl,cfastr,cfast,bnt
      real, dimension(  4,nmax  ) :: speeds
      real, dimension( 10,nmax  ) :: fx
      real, dimension(7,7       ) :: unit,rhed,lhed
      real, dimension(    nmax  ) :: einfeldt_fix,indx_ein
      integer :: numein
      real   , dimension(  nmax)::a,cf,cs,va,vax,alpha,blpha
      real                      ::sumr,suml,sum1,nu,div1
      real                      ::bminus,bplus,denom,coef0,coef1,coef2
      real                      ::et,eb,ek,b2f,bdotb0,bdotul, vt
      real                      ::bdotb0l,bdotb0r
      real                      ::sdif,fmax,bm,bp
      real,parameter            ::const1=0.25
      logical, dimension(  nmax)::lrare
      if(idir.eq.1)then
        delta_x=delx
      elseif(idir.eq.2)then
        delta_x=dely
      elseif(idir.eq.3)then
        delta_x=delz
      endif
       do i=n1,n2
        vl(1:5,i) = vml(1:5,i)
        vr(1:5,i) = vmr(1:5,i)
        vl(6:7,i) = vml(7:8,i)
        vr(6:7,i) = vmr(7:8,i)

        vdl(1:5,i) = vl(1:5,i) 
        vdl(6:7,i) = vl(6:7,i) + bf(2:3,i)
        vdr(1:5,i) = vr(1:5,i) 
        vdr(6:7,i) = vr(6:7,i) + bf(2:3,i)
        bnt(    i) = bn(i)     + bf(  1,i)

        vl(1,i)   = max(rmin,vl(1,i))
        vr(1,i)   = max(rmin,vr(1,i))
        vl(5,i)   = max(pmin,vl(5,i))
        vr(5,i)   = max(pmin,vr(5,i))
       enddo

       do i=n1,n2
        vl(1:5,i) = vml(1:5,i)
        vr(1:5,i) = vmr(1:5,i)
        vl(6:7,i) = vml(7:8,i)
        vr(6:7,i) = vmr(7:8,i)

        vdl(1:5,i) = vl(1:5,i) 
        vdl(6:7,i) = vl(6:7,i) + bf(2:3,i)
        vdr(1:5,i) = vr(1:5,i) 
        vdr(6:7,i) = vr(6:7,i) + bf(2:3,i)
        bnt(    i) = bn(i)     + bf(  1,i)

       enddo

c
c-----averages
c
      do i=n1,n2
        bdotb0l     = bn(i)*bf(1,i)+vl(6,i)*bf(2,i)+vl(7,i)*bf(3,i)
        ekl         = 0.5*vl(1,i)*(vl(2,i)**2+vl(3,i)**2+vl(4,i)**2)
        b2l         = (bn(i)**2+vl(6,i)**2+vl(7,i)**2)*pie8i
        ul(1,i)     = vl(1,i)
        ul(2,i)     = vl(1,i)*vl(2,i)
        ul(3,i)     = vl(1,i)*vl(3,i)
        ul(4,i)     = vl(1,i)*vl(4,i)
        ul(5,i)     = ekl+b2l+vl(5,i)*gammam1i
        ul(6,i)     = vl(6,i)
        ul(7,i)     = vl(7,i)
        etl         = ul(5,i)
        rli         = 1.0/vl(1,i)

        bdotb0r     = bn(i)*bf(1,i)+vr(6,i)*bf(2,i)+vr(7,i)*bf(3,i)
        ekr         = 0.5*vr(1,i)*(vr(2,i)**2+vr(3,i)**2+vr(4,i)**2)
        b2r         = (bn(i)**2+vr(6,i)**2+vr(7,i)**2)*pie8i
        ur(1,i)     = vr(1,i)
        ur(2,i)     = vr(1,i)*vr(2,i)
        ur(3,i)     = vr(1,i)*vr(3,i)
        ur(4,i)     = vr(1,i)*vr(4,i)
        ur(5,i)     = ekr+b2r+vr(5,i)*gammam1i
        ur(6,i)     = vr(6,i)
        ur(7,i)     = vr(7,i)
        etr         = ur(5,i)
        rri         = 1.0/vr(1,i)

#ifdef ROE_AVERAGE
!
!-------Original Roe Averaging !!
!

        sqrhl       = sqrt(vl(1,i))
        sqrhr       = sqrt(vr(1,i))
        rh          = sqrhr*sqrhl
        rhi         = 1.0/(sqrhr+sqrhl)
        uh          = (sqrhl*vl(2,i)+sqrhr*vr(2,i))*rhi
        vh          = (sqrhl*vl(3,i)+sqrhr*vr(3,i))*rhi
        wh          = (sqrhl*vl(4,i)+sqrhr*vr(4,i))*rhi
        vh2         = (uh**2+vh**2+wh**2)
        hr          = (etr+vr(5,i)+pie8i*(b2r+2.0*bdotb0r))*rri
        hl          = (etl+vl(5,i)+pie8i*(b2l+2.0*bdotb0l))*rli
        hh          = (sqrhl*hl+sqrhr*hr)*rhi
c       ph          = (sqrhl*vl(5,i)+sqrhr*vr(5,i))*rhi
        bvh         = (sqrhl*vl(6,i)+sqrhr*vr(6,i))*rhi
        bwh         = (sqrhl*vl(7,i)+sqrhr*vr(7,i))*rhi
        bh2         = bn(i)**2+bvh**2+bwh**2
        bh2         = bh2+2.0*(bf(1,i)*bn(i)+bf(2,i)*bvh+bf(3,i)*bwh)
        a2h         = gammam1*(hh-0.5*vh2+pie4i*bh2/rh)
        ph          = a2h*rh/gamma
        bvh         = bvh+bf(2,i)
        bwh         = bwh+bf(3,i)

#else
        rh          = 0.5*(vl(1,i)+vr(1,i))
        uh          = 0.5*(vl(2,i)+vr(2,i))
        vh          = 0.5*(vl(3,i)+vr(3,i))
        wh          = 0.5*(vl(4,i)+vr(4,i))
        ph          = 0.5*(vl(5,i)+vr(5,i))
        bvh         = 0.5*(vl(6,i)+vr(6,i))
        bwh         = 0.5*(vl(7,i)+vr(7,i))
        bvh         = bvh+bf(2,i)
        bwh         = bwh+bf(3,i)
#endif        
        vbar(  1,i) = rh
        vbar(  2,i) = uh
        vbar(  3,i) = vh
        vbar(  4,i) = wh
        vbar(  5,i) = ph
        vbar(  6,i) = bvh
        vbar(  7,i) = bwh
      enddo

c
c-----compute delta v
c
      do i=n1,n2
       delu(:7,i) = ur(:7,i)-ul(:7,i)
       delv(:7,i) = vr(:7,i)-vl(:7,i)
       dvy(i)     = vr(3,i)-vl(3,i)
       dvz(i)     = vr(4,i)-vl(4,i)
      enddo
      vf(:7,:)  = 0.0

c
c-----get the primitive eigenvectors and eigenweights for the average state
c
      call eigenvalues(n1,n2,vdl ,bnt,eval_l,cfastl)
      call eigenvalues(n1,n2,vdr ,bnt,eval_r,cfastr)

#ifndef PRIMITIVE_FLUXES
      call mhd_eigenvectors(rhe,lhe ,rhec,lhec,eigenval,vbar ,bnt   ,
     &                      wt ,delv,2   ,n1  ,n2      ,cfast,speeds,
     &                      dvy, dvz,idir)
#else
      call mhd_eigenvectors(rhe,lhe ,rhec,lhec,eigenval,vbar ,bnt   ,
     &                      wt ,delv,1   ,n1  ,n2      ,cfast,speeds,
     &                      dvy, dvz,idir)
#endif
#ifdef HARTEN_ENROPY_FIX
c
c----sonic point fix
c
      do i=n1,n2
       do k=1, 7
         eps_zone(k,i) = max(smin,(eigenval(k,i)-  eval_l(k,i)),
     &                            (  eval_r(k,i)-eigenval(k,i)))
         sign_eval(k,i)= sign(1.0,eigenval(k,i))
         eigenval(k,i) = abs(eigenval(k,i) )

         if(eigenval(k,i).le.eps_zone(k,i))
     &      eigenval(k,i) = 0.5*(eigenval(k,i)**2/eps_zone(k,i) 
     &                                          + eps_zone(k,i))
          wt(k,i)       = 0.5*eigenval(k,i)*wt(k,i)
          eigenval(k,i) = eigenval(k,i)*sign_eval(k,i) 
        enddo
      enddo
#endif /*HARTEN_ENROPY_FIX*/
 
#ifdef HLLEM_ENTROPY_FIX
c
c-----make Roe's method both positively conserving and entropy satisfying
c
      do i=n1,n2
       bplus  = max(eigenval(7,i),eval_r(7,i),0.0)
       bminus = min(eigenval(1,i),eval_l(1,i),0.0)
       denom  = cfast(i)+0.5*abs(max(eigenval(7,i),eval_r(7,i))+
     &                           min(eigenval(1,i),eval_l(1,i)))
       delta  = cfast(i)/denom
       coef0  = 1.0/(bplus-bminus)
       coef1  = (bplus+bminus)*coef0
       coef2  = (bplus*bminus)*coef0*2.0
       eigenval(1,i) = coef1*eigenval(1,i) - coef2
       eigenval(2,i) = coef1*eigenval(2,i) - coef2*(1.0-delta)
       eigenval(3,i) = coef1*eigenval(3,i) - coef2
       eigenval(4,i) = coef1*eigenval(4,i) - coef2*(1.0-delta)
       eigenval(5,i) = coef1*eigenval(5,i) - coef2
       eigenval(6,i) = coef1*eigenval(6,i) - coef2*(1.0-delta)
       eigenval(7,i) = coef1*eigenval(7,i) - coef2
      enddo
#endif /*HLLEM_ENTROPY_FIX*/

c
c-----compute resolved state primitives
c

      vbar_r(:,:) = vr(:,:)
      vbar_l(:,:) = vl(:,:)
      do i=n1,n2
       do j=1,7
        do k=1,7
         if(eigenval(k,i).lt.-smin)then
          vbar_l(j,i)=vbar_l(j,i)+wt(k,i)*rhe(j,k,i)
         endif
         kk = 8-k
         if(eigenval(kk,i).gt.smin)then
          vbar_r(j,i)=vbar_r(j,i)-wt(kk,i)*rhe(j,kk,i)
         endif
        enddo
       enddo
      enddo
      do i=n1,n2
        if(eigenval(4,i)<-smin)then
         do k=1,7
          vf(k,i) = vbar_r(k,i)
         enddo
        elseif(eigenval(4,i)>smin)then
         do k=1,7
          vf(k,i) = vbar_l(k,i)
         enddo
        else
         do k=1,7
          vf(k,i) = 0.5*(vbar_l(k,i)+vbar_r(k,i))
         enddo
        endif
      enddo

#ifndef PRIMITIVE_FLUXES
c
c----construct left and right conserved variables and fluxes
c
      do i=n1,n2
c
c----left
c
        ul(1,i)     = vl(1,i)
        ul(2,i)     = ul(1,i)*vl(2,i)
        ul(3,i)     = ul(1,i)*vl(3,i)
        ul(4,i)     = ul(1,i)*vl(4,i)
        ek          = 0.5*vl(1,i)*(vl(2,i)**2+vl(3,i)**2+vl(4,i)**2)
        eb          = (bn(i)**2+vl(6,i)**2+vl(7,i)**2)*pie8i
        ul(5,i)     = ek+eb+vl(5,i)*gammam1i
        ul(6,i)     = vl(6,i)
        ul(7,i)     = vl(7,i)

        bdotul      = bn(i)*vl(2,i)+vl(6,i)*vl(3,i)+vl(7,i)*vl(4,i)
        bdotb0      = bn(i)*bf(1,i)+vl(6,i)*bf(2,i)+vl(7,i)*bf(3,i)
        bylt        = bf(2,i)+vl(6,i)
        bzlt        = bf(3,i)+vl(7,i)
        b2l         = (vl(6,i)**2+vl(7,i)**2-bn(i)**2)*pie8i

        fl( 1,i) = ul(2,i)
        fl( 2,i) = ul(2,i)*vl(2,i) + vl(5,i) 
        fl( 3,i) = ul(2,i)*vl(3,i) 
        fl( 4,i) = ul(2,i)*vl(4,i) 
        fl( 5,i) = (ul(5,i)+vl(5,i)+eb+pie4i*bdotb0)*vl(2,i)
     &               -bdotul*bnt(i)*pie4i
        fl( 6,i) = vl(2,i)*bylt-vl(3,i)*bnt(i)
        fl( 7,i) = vl(2,i)*bzlt-vl(4,i)*bnt(i)
        fl( 8,i) = b2l + pie4i*(bdotb0-2.0*bf(1,i)*bn(i))
        fl( 9,i) =-(bnt(i)*vl(6,i)+bn(i)*bf(2,i))*pie4i
        fl(10,i) =-(bnt(i)*vl(7,i)+bn(i)*bf(3,i))*pie4i
c
c----right
c
        ur(1,i)     = vr(1,i)
        ur(2,i)     = ur(1,i)*vr(2,i)
        ur(3,i)     = ur(1,i)*vr(3,i)
        ur(4,i)     = ur(1,i)*vr(4,i)
        ek          = 0.5*vr(1,i)*(vr(2,i)**2+vr(3,i)**2+vr(4,i)**2)
        eb          = (bn(i)**2+vr(6,i)**2+vr(7,i)**2)*pie8i
        ur(5,i)     = ek+eb+vr(5,i)*gammam1i
        ur(6,i)     = vr(6,i)
        ur(7,i)     = vr(7,i)

        bdotur      = bn(i)*vr(2,i)+vr(6,i)*vr(3,i)+vr(7,i)*vr(4,i)
        bdotb0      = bn(i)*bf(1,i)+vr(6,i)*bf(2,i)+vr(7,i)*bf(3,i)
        byrt        = bf(2,i)+vr(6,i)
        bzrt        = bf(3,i)+vr(7,i)
        b2r         = (vr(6,i)**2+vr(7,i)**2-bn(i)**2)*pie8i

        fr( 1,i) = ur(2,i)
        fr( 2,i) = ur(2,i)*vr(2,i) + vr(5,i) 
        fr( 3,i) = ur(2,i)*vr(3,i) 
        fr( 4,i) = ur(2,i)*vr(4,i) 
        fr( 5,i) = (ur(5,i)+vr(5,i)+eb+pie4i*bdotb0)*vr(2,i)
     &               -bdotur*bnt(i)*pie4i
        fr( 6,i) = vr(2,i)*byrt-vr(3,i)*bnt(i)
        fr( 7,i) = vr(2,i)*bzrt-vr(4,i)*bnt(i)
        fr( 8,i) = b2r + pie4i*(bdotb0-2.0*bf(1,i)*bn(i))
        fr( 9,i) = -(bnt(i)*vr(6,i)+bn(i)*bf(2,i))*pie4i
        fr(10,i) = -(bnt(i)*vr(7,i)+bn(i)*bf(3,i))*pie4i

        do j=1,7
         fld(j,i) = fl(j,i)
         frd(j,i) = fr(j,i)
      enddo

        do j=2,4
         fld(j,i) = fld(j,i)+fl(j+6,i)
         frd(j,i) = frd(j,i)+fr(j+6,i)
       enddo

      enddo

      do i=n1,n2
       do k=1,10
        fx(k,i)       = 0.5*(fl(k,i)+fr(k,i))
       enddo
      enddo
c
c-----Construct fluxes from resolved state
c
      do i=n1,n2
       do j=1,7
        do k=1,7
         fx(j,i)=fx(j,i)-0.5*abs(eigenval(k,i))*wt(k,i)*rhec(j,k,i)
         enddo
         enddo
       enddo

#else

      do i=n1,n2
       bdotb0   = bn(i)*bf(1,i)+vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)
       bdotul   = bn(i)*vf(2,i)+vf(6,i)*vf(3,i)+vf(7,i)*vf(4,i)
       b2f      = (vf(6,i)**2+vf(7,i)**2-bn(i)**2)*pie8i
       ek       = 0.5*vf(1,i)*(vf(2,i)**2+vf(3,i)**2+vf(4,i)**2)
       eb       = (vf(6,i)**2+vf(7,i)**2+bn(i)**2)
       eb       = eb*pie8i
       et       = ek+eb+vf(5,i)*gammam1i
       byt      = bf(2,i)  + vf(6,i)
       bzt      = bf(3,i)  + vf(7,i)
       fx(1,i)  = vf(1,i)*vf(2,i)
       fx(2,i)  = fx(1,i)*vf(2,i) + vf(5,i) 
       fx(3,i)  = fx(1,i)*vf(3,i) 
       fx(4,i)  = fx(1,i)*vf(4,i) 
       fx(5,i)  = (et+vf(5,i)+eb+pie4i*bdotb0)*vf(2,i)
     &               -bdotul*bnt(i)*pie4i
       fx(6,i)  = vf(2,i)*byt-vf(3,i)*bnt(i)
       fx(7,i)  = vf(2,i)*bzt-vf(4,i)*bnt(i)

       fx( 8,i) = b2f + pie4i*(vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)-bf(1,i)*bn(i))
       fx( 9,i) =-(bnt(i)*vf(6,i)+bn(i)*bf(2,i))*pie4i
       fx(10,i) =-(bnt(i)*vf(7,i)+bn(i)*bf(3,i))*pie4i

      enddo
#endif /*PRIMITIVE_FLUXES*/

#ifdef LAPIDIUS
c-----Lapidius AV

      nu = 0.1
      do i=n1,n2
       div1     = nu*min( divv(i),0.)
       fx(:5,i) = fx(:5,i) + div1*(uc(:5,i)-uc(:5,i-1))
       fx( 6,i) = fx( 6,i) + div1*(uc( 7,i)-uc( 7,i-1))
       fx( 7,i) = fx( 7,i) + div1*(uc( 8,i)-uc( 8,i-1))
      enddo
#endif /*LAPIDIUS*/

#ifdef STRONG_RAREFACTION_WAVE_FIX
      lrare(:)=.false.
      do i=n1,n2
       sdif = vr(2,i)-vl(2,i)
       fmax = eigenval(7,i) - eigenval(4,i)
       if(sdif>=const1*fmax )	then
	lrare(i)= .true.
       endif
      enddo
      do i=n1,n2
       if(lrare(i))then
c        write(*,*)'lrare(i)=',lrare(i)
        bm   = min(eval_l(1,i),eval_l(4,i),eigenval(1,i),0.0)
        bp   = min(eval_r(7,i),eval_r(4,i),eigenval(7,i),0.0)
        difi = 1.0/(bp-bm)
        dif  = bp*bm
        do j=1,7
         uf(j,i)=difi*((bp*ur(j,i)-bm*ul(j,i))-    (frd(j,i)-fld(j,i)))
         fx(j,i)=difi*((bp*fr(j,i)-bm*fl(j,i))+dif*( ur(j,i)- ul(j,i)))
        enddo 
       endif
      enddo


#endif /*STRONG_RAREFACTION_WAVE_FIX*/
#ifdef INTERIOR_BOUNDARYXX
C     hardwall Boundary
      
      if(ifaces_1d.gt.0) then
      do ii =1, ifaces_1d
       ip = iface_1d(1,ii)
       f  = iface_1d(2,ii)
       if(f.eq.2)then         ! boundary face to your right
        i        = ip+1
#ifdef NO_SLIP
        bdotb0   = bn(i)*bf(1,i)+vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)
        bdotul   = bn(i)*vf(2,i)+vf(6,i)*vf(3,i)+vf(7,i)*vf(4,i)
        b2f      = (vf(6,i)**2+vf(7,i)**2-bn(i)**2)*pie8i
        ek       = 0.0
        eb       = (vf(6,i)**2+vf(7,i)**2+bn(i)**2)
        eb       = eb*pie8i
        byt      = bf(2,i)  + vf(6,i)
        bzt      = bf(3,i)  + vf(7,i)

        fx(1,i)  = 0.0
        fx(2,i)  = vf(5,i) 
        fx(3,i)  = 0.0 
        fx(4,i)  = 0.0 
        fx(5,i)  = 0.0
        fx(6,i)  = 0.0
        fx(7,i)  = 0.0

        fx( 8,i) = b2f + pie4i*(vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)-bf(1,i)*bn(i))
        fx( 9,i) =-(bnt(i)*vf(6,i)+bn(i)*bf(2,i))*pie4i
        fx(10,i) =-(bnt(i)*vf(7,i)+bn(i)*bf(3,i))*pie4i
#else
        bdotb0   = bn(i)*bf(1,i)+vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)
        bdotul   =               vf(6,i)*vf(3,i)+vf(7,i)*vf(4,i)
        b2f      = (vf(6,i)**2+vf(7,i)**2-bn(i)**2)*pie8i
        ek       = 0.0
        eb       = (vf(6,i)**2+vf(7,i)**2+bn(i)**2)
        eb       = eb*pie8i
        byt      = bf(2,i)  + vf(6,i)
        bzt      = bf(3,i)  + vf(7,i)

        fx(1,i)  = 0.0
        fx(2,i)  = vf(5,i) 
        fx(3,i)  = 0.0 
        fx(4,i)  = 0.0 
        fx(5,i)  = -bdotul*bnt(i)*pie4i
        fx(6,i)  = -vf(3,i)*bnt(i)
        fx(7,i)  = -vf(4,i)*bnt(i)

        fx( 8,i) = b2f + pie4i*(vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)-bf(1,i)*bn(i))
        fx( 9,i) =-(bnt(i)*vf(6,i)+bn(i)*bf(2,i))*pie4i
        fx(10,i) =-(bnt(i)*vf(7,i)+bn(i)*bf(3,i))*pie4i
#endif
       else
        i        = ip
#ifdef NO_SLIP
        bdotb0   = bn(i)*bf(1,i)+vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)
        bdotul   = bn(i)*vf(2,i)+vf(6,i)*vf(3,i)+vf(7,i)*vf(4,i)
        b2f      = (vf(6,i)**2+vf(7,i)**2-bn(i)**2)*pie8i
        ek       = 0.0
        eb       = (vf(6,i)**2+vf(7,i)**2+bn(i)**2)
        eb       = eb*pie8i
        byt      = bf(2,i)  + vf(6,i)
        bzt      = bf(3,i)  + vf(7,i)

        fx(1,i)  = 0.0
        fx(2,i)  = vf(5,i) 
        fx(3,i)  = 0.0 
        fx(4,i)  = 0.0 
        fx(5,i)  = 0.0
        fx(6,i)  = 0.0
        fx(7,i)  = 0.0

        fx( 8,i) = b2f + pie4i*(vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)-bf(1,i)*bn(i))
        fx( 9,i) =-(bnt(i)*vf(6,i)+bn(i)*bf(2,i))*pie4i
        fx(10,i) =-(bnt(i)*vf(7,i)+bn(i)*bf(3,i))*pie4i
#else
        bdotb0   = bn(i)*bf(1,i)+vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)
        bdotul   =               vf(6,i)*vf(3,i)+vf(7,i)*vf(4,i)
        b2f      = (vf(6,i)**2+vf(7,i)**2-bn(i)**2)*pie8i
        ek       = 0.0
        eb       = (vf(6,i)**2+vf(7,i)**2+bn(i)**2)
        eb       = eb*pie8i
        byt      = bf(2,i)  + vf(6,i)
        bzt      = bf(3,i)  + vf(7,i)

        fx(1,i)  = 0.0
        fx(2,i)  = vf(5,i) 
        fx(3,i)  = 0.0 
        fx(4,i)  = 0.0 
        fx(5,i)  = -bdotul*bnt(i)*pie4i
        fx(6,i)  = -vf(3,i)*bnt(i)
        fx(7,i)  = -vf(4,i)*bnt(i)

        fx( 8,i) = b2f + pie4i*(vf(6,i)*bf(2,i)+vf(7,i)*bf(3,i)-bf(1,i)*bn(i))
        fx( 9,i) =-(bnt(i)*vf(6,i)+bn(i)*bf(2,i))*pie4i
        fx(10,i) =-(bnt(i)*vf(7,i)+bn(i)*bf(3,i))*pie4i
#endif

       endif
      enddo
      endif

#endif /* INTERIOR_BOUNDARY */

!
!-----rotate fluxes here
!
      if(idir.eq.1)then

       do i = n1,n2
        flux( 1,i) = fx( 1,i)
        flux( 2,i) = fx( 2,i) 
        flux( 3,i) = fx( 3,i)
        flux( 4,i) = fx( 4,i)
        flux( 5,i) = fx( 5,i)
        flux( 6,i) = 0.0
        flux( 7,i) = fx( 6,i) 
        flux( 8,i) = fx( 7,i)
        flux( 9,i) = fx( 8,i) 
        flux(10,i) = fx( 9,i)
        flux(11,i) = fx(10,i)
!
!------primitives at faces
!
        vfout( :5,i) = vf(:5,i)
        vfout(  6,i) = bn(i)
        vfout(7:8,i) = vf(6:7,i)
        
        
!
       enddo

      elseif(idir.eq.2)then

       do i = n1,n2

        flux( 1,i) = fx(1,i)
        flux( 3,i) = fx(2,i) 
        flux( 4,i) = fx(3,i)
        flux( 2,i) = fx(4,i)
        flux( 5,i) = fx(5,i)
        flux( 7,i) = 0.0 
        flux( 8,i) = fx( 6,i)
        flux( 6,i) = fx( 7,i) 
        flux(10,i) = fx( 8,i) 
        flux(11,i) = fx( 9,i)
        flux( 9,i) = fx(10,i)
!
!------primitives at faces
!
        vfout(1,i) = vf(1,i)
        vfout(3,i) = vf(2,i)
        vfout(4,i) = vf(3,i)
        vfout(2,i) = vf(4,i)
        vfout(5,i) = vf(5,i)
        vfout(7,i) = bn(i)
        vfout(8,i) = vf(6,i)
        vfout(6,i) = vf(7,i)
!
       enddo

      else

       do i = n1,n2

        flux( 1,i) = fx( 1,i)
        flux( 4,i) = fx( 2,i)
        flux( 2,i) = fx( 3,i)
        flux( 3,i) = fx( 4,i)
        flux( 5,i) = fx( 5,i)
        flux( 8,i) = 0.0
        flux( 6,i) = fx( 6,i) 
        flux( 7,i) = fx( 7,i) 
        flux(11,i) = fx( 8,i) 
        flux( 9,i) = fx( 9,i)
        flux(10,i) = fx(10,i)
!
!------primitives at faces
!
        vfout(1,i) = vf(1,i)
        vfout(4,i) = vf(2,i)
        vfout(2,i) = vf(3,i)
        vfout(3,i) = vf(4,i)
        vfout(5,i) = vf(5,i)
        vfout(8,i) = bn(i)
        vfout(6,i) = vf(6,i)
        vfout(7,i) = vf(7,i)

        

       enddo

      endif
c
      return
      end subroutine roe_mhd

