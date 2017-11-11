      subroutine polint(xa,ya,n,x,y,dy)

      integer n,nmax
      real    dy,x,y,xa(n),ya(n)
      parameter(nmax=10)
      integer i,m,ns
      real den,dif,dift,ho,hp,w,c(nmax),d(nmax)

      ns = 1
      dif = abs(x-xa(1))

      do i=1,n
         dift = abs(x-xa(i))
         if(dift.lt.dif) then
           ns = i
           dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
      enddo

      y = ya(ns)
      ns = ns-1

      do m=1,n-1
        do i=1,n-m
           ho = xa(i)-x
           hp = xa(i+m)-x
           w = c(i+1)-d(i)
           den = ho - hp

! error trap : error occurs if two xa values are identical to roundoff
           if(den.eq.0) stop 'failure in polint'

           den = w/den
           d(i) = hp*den
           c(i) = ho*den
        enddo
        if(2*ns.lt.n-m) then
           dy = c(ns+1)
        else
           dy = d(ns)
           ns = ns-1
        endif
        y = y+dy
      enddo

      return
      end
