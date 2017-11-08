      subroutine read_control(infile, maxstp,
     +       lgraph,nstep_graphics,restart,
     +       nstep_dump,nstep_print,
     +       ltecplot,nstep_tecplot,
     +       rmin,pmin,smin,
     +       cfl,dt,deltim,
     +       ctore,ctode,epsil,nuref,
     +       luref, muref,
     +       nmovie,
     +       movie_si,
     +       lmovie,
     +       arfact2)

      character*(*) infile
      logical :: lgraph, restart,lrefine_again
      integer :: nstep_graphics
      integer :: nstep_dump
      integer :: ltecplot
      integer :: nstep_tecplot
      real    :: rmin, pmin, smin
      real    :: cfl
      real    :: dt
      real    :: ctore,ctode,epsil
      integer :: nuref,luref(muref)
      integer :: nmovie, movie_si
      logical :: lmovie
      real    :: arfact2

      open(unit=2,file=infile,status = 'old', err=500)
      read(2,*,err=501)
      read(2,*,err=501) maxstp
      read(2,*,err=501)
      read(2,*,err=501)lgraph,nstep_graphics,restart,
     &     nstep_dump,nstep_print
      read(2,*,err=501)
      read(2,*,err=501)ltecplot,nstep_tecplot
      read(2,*,err=501)
      read(2,*,err=501)rmin,pmin,smin
      read(2,*,err=501)
      read(2,*,err=501)cfl,dt,deltim
      read(2,*,err=501)
      read(2,*,err=501)ctore,ctode,epsil,nuref
      read(2,*,err=501)
      read(2,*,err=501)(luref(i),i=1,nuref)
      read(2,*,err=501)
      read(2,*,err=501)nmovie
      read(2,*,err=501)
      read(2,*,err=501)movie_si
      read(2,*,err=501)
      read(2,*,err=501)lmovie
      read(2,*,err=501)
      read(2,*,err=501)arfact2
      close(2,err=502)
      return

 500  print*,'Unable to open file ', infile
      return

 501  print*,'Read error for file', infile
      return

 502  print*,'Error closing file', infile
      return

      end
