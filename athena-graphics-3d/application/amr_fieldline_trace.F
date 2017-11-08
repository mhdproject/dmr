#include "paramesh_preprocessor.fh"
       subroutine amr_fieldline_trace(loop)
       use physcons
       use paramesh_dimensions
       use physicaldata

! include file defining the tree
       use tree


      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */




#ifdef MPI_USED
       include "mpif.h"
      integer :: tag_offset,max_blks_sent
#endif /* MPI_USED */


! include file required for shmem library.
#include "amr_shmem.fh"

#include "fieldline_tracer.fh"

! local amr variables
        integer nprocs,errcode
        integer shmem_my_pe,shmem_n_pes
        logical lcc,lfc,lzeroflux,newchildt(maxblocks_tr),ldebug


#ifdef TECPLOT_IO
        character (len=6) :: filenumber
#endif /* TECPLOT_IO */

!---------------------------------------------------------------
!
        mype = shmem_my_pe()
        nprocs = shmem_n_pes()

!mk5
!          write(*,*)'amr_fieldline_trace: ', mype

!---------------------------------------------------------------
! exchange guardcell information
        iopt = 1
        nlayers = nguard
#ifdef MPI_USED
!        call mpi_amr_guardcell(mype,iopt,nlayers,nprocs)
        call amr_guardcell(mype,iopt,nlayers)
#else
        call     amr_guardcell(mype,iopt,nlayers)
#endif /*MPI_USED*/

!       write(*,*) 'calling field_at_corners'
        call field_at_cell_corners
!       write(*,*) 'exited field_at_corners'

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


#ifdef TECPLOT_IO
        if(mype.eq.0) then

        write(filenumber,"(I6.6)") loop
        iout = 10
        open(unit=iout,status='new',file='line.dat.'//filenumber)

        write(iout,*) 'TITLE = "Field-lines  time = ',time,'"'
        write(iout,*) 'VARIABLES = "x", "y", "z"'

        endif                                ! mype if test
#endif /* TECPLOT_IO */
!----------------------------------------
        ldebug = .false.
! fieldlines from inner radial boundary
!mk5
!       write(*,*)'field_trace: grid_xmax,grid_xmin ',
!    &grid_xmax,grid_xmin
!       write(*,*)'field_trace: grid_ymax,grid_ymin ',
!    &grid_ymax,grid_ymin
        delta_x = (grid_xmax-grid_xmin)/float(nflines)
        do jline=1,nflines,2
#if N_DIM == 3
        x0(jline) = grid_xmin+delta_x*real(jline-1)
        y0(jline) = 0.
        y0(jline) = 0.
        z0(jline) = grid_zmax
#else
        x0(jline) = grid_xmin+delta_x*real(jline-1)
!dss        x0(jline) = -3822.+100.*real(jline-1)
!dss        x0(jline) = -700.+25.*real(jline-1)
        y0(jline) = 0.
!mk5
!       y0(jline) = grid_ymax
!dss     y0(jline) = 5500.
        z0(jline) = 0.
#endif

! fieldlines in +ve direction
        call fieldlines(x0(jline),y0(jline),z0(jline),
     .                  xf(1,jline),yf(1,jline),zf(1,jline),
     .                  1.,1,npts,mype,ldebug)
        fline_npts(jline) = npts
!mk5
!        if(mype.eq.1) 
!    &   write(*,*)'# of fieldpts.,npts,mype ',npts,mype
#ifdef TECPLOT_IO
        write(iout,*) 'ZONE I=',npts,' F=POINT'
        do i=1,npts
        write(iout,100) xf(i,jline),yf(i,jline),zf(i,jline)
        enddo
#endif /* TECPLOT_IO */

! fieldlines in -ve direction
        x0(jline+1) = x0(jline)
c        y0(jline+1) = y0(jline)
        y0(jline+1) =-y0(jline)
        z0(jline+1) = z0(jline)
        call fieldlines(x0(jline+1),y0(jline+1),z0(jline+1),
     .                  xf(1,jline+1),yf(1,jline+1),zf(1,jline+1),
     .                  -1.,1,npts,mype,ldebug)
        fline_npts(jline+1) = npts
!mk5
!        if(mype.eq.1) 
!    &   write(*,*)'# of fieldpts.,npts,mype ',npts,mype
#ifdef TECPLOT_IO
        write(iout,*) 'ZONE I=',npts,' F=POINT'
        do i=1,npts
        write(iout,100) xf(i,jline+1),yf(i,jline+1),zf(i,jline+1)
        enddo
#endif /* TECPLOT_IO */

        enddo

#ifdef TECPLOT_IO
100     format(3(2x,1pe15.7))
#endif /* TECPLOT_IO */
!----------------------------------------


#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


!       call draw_bounding_boxes(iout)

#ifdef TECPLOT_IO
        if(mype.eq.0) then
        close(unit=iout)
        endif                                ! mype if test
#endif /* TECPLOT_IO */

        return
        end
