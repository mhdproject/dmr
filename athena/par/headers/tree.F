!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------
! tree module



      module tree

      use paramesh_dimensions

      private

      public :: maxblocks_tr
      public :: nchild, nfaces, mchild, mfaces, mdim,  mflags
      public :: mboundaries,nboundaries

! Block limit to be used in manipulating tree when modifying the grid.
#ifndef LIBRARY
      integer, parameter :: maxblocks_tr=4*maxblocks
#else
      integer :: maxblocks_tr
#endif
! Number of children of a node
#ifndef LIBRARY
      integer, parameter :: nchild=2**ndim
#else
      integer :: nchild
#endif

! Number of faces on a grid block
#ifndef LIBRARY
      integer, parameter :: nfaces=2*ndim
#else
      integer :: nfaces
#endif

! Parameters used to define array sizes
      integer, parameter :: mdim=3,mchild=2**mdim,mfaces=2*mdim

! Parameters used to declare the number of block marker flags needed
! --<< USER EDIT >>--
      integer, parameter :: mflags=MFLAGS

! Parameters used to declare the number of boundary regions where boundary
! conditions are to be applied.  Typically: 2*ndim
#ifndef LIBRARY
! --<< USER EDIT >>--
      integer, parameter :: nboundaries=NBOUNDARIES
#else
      integer :: nboundaries
#endif

      integer, parameter :: mboundaries=2*mdim

      public :: neigh,child,which_child,type_of_neigh,type_of_child
      public :: parent,lrefine,lnblocks,new_lnblocks
      public :: nodetype,empty,bflags,newchild,derefine,refine
      public :: stay,work_block,coord,bsize,bnd_box
      public :: grid_xmin,grid_xmax,grid_ymin,grid_ymax
      public :: grid_zmin,grid_zmax
      public :: lrefine_max,lrefine_min
      public :: level_cell_sizes

! Variables for storing tree datastructure
#ifndef LIBRARY
      integer, save :: neigh(2,mfaces,maxblocks_tr)
      integer, save :: child(2,mchild,maxblocks_tr)
      integer, save :: which_child(maxblocks_tr)
      integer, save :: type_of_neigh(mfaces,maxblocks_tr)
      integer, save :: type_of_child(mchild,maxblocks_tr)
      integer, save :: parent(2,maxblocks_tr),lrefine(maxblocks_tr)
#else
      integer, allocatable, save :: neigh(:,:,:)
      integer, allocatable, save :: child(:,:,:)
      integer, allocatable, save :: which_child(:)
      integer, allocatable, save :: type_of_neigh(:,:)
      integer, allocatable, save :: type_of_child(:,:)
      integer, allocatable, save :: parent(:,:),lrefine(:)
#endif
      integer, save :: lnblocks,new_lnblocks
#ifndef LIBRARY
      integer, save :: nodetype(maxblocks_tr)
      integer, save :: empty(maxblocks_tr),bflags(mflags,maxblocks_tr)
      logical, save :: newchild(maxblocks_tr)
      logical, save :: derefine(maxblocks_tr),refine(maxblocks_tr)
      logical, save :: stay(maxblocks_tr)
#else
      integer, allocatable, save :: nodetype(:)
      integer, allocatable, save :: empty(:),bflags(:,:)
      logical, allocatable, save :: newchild(:)
      logical, allocatable, save :: derefine(:),refine(:)
      logical, allocatable, save :: stay(:)
#endif
#ifndef LIBRARY
      real, save :: work_block(maxblocks_tr)
      real, save :: coord(mdim,maxblocks_tr)
      real, save :: bsize(mdim,maxblocks_tr)
      real, save :: bnd_box(2,mdim,maxblocks_tr)
#else
      real, allocatable, save :: work_block(:)
      real, allocatable, save :: coord(:,:)
      real, allocatable, save :: bsize(:,:)
      real, allocatable, save :: bnd_box(:,:,:)
#endif
      real,save :: grid_xmin,grid_xmax
      real,save :: grid_ymin,grid_ymax
      real,save :: grid_zmin,grid_zmax

#ifndef LIBRARY
      real,save :: level_cell_sizes(mdim,maxlevels)
#else
      real, allocatable, save :: level_cell_sizes(:,:)
#endif
      integer, save :: lrefine_max,lrefine_min

! flag to record grid change
      public :: grid_changed, grid_analysed_mpi
      integer, save :: grid_changed, grid_analysed_mpi

! added for surrblks calculation with mpi
      public :: boundary_box,boundary_index
#ifndef LIBRARY
      real, save    :: boundary_box(2,mdim,mboundaries)
      integer, save :: boundary_index(mboundaries)
#else
      real, allocatable,save    :: boundary_box(:,:,:)
      integer, allocatable, save :: boundary_index(:)
#endif

! added for use with mpi block buffering
      public :: strt_buffer,last_buffer
      public :: strt_buffer_tree,last_buffer_tree
      public :: laddress,surr_blks
#ifdef SAVE_MORTS
      public :: surr_morts
#endif
      integer, save :: strt_buffer,last_buffer
      integer, save :: strt_buffer_tree,last_buffer_tree
#ifndef LIBRARY
      integer, save :: surr_blks(3,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
#ifdef SAVE_MORTS
      integer, save :: surr_morts(6,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
#endif
      integer, save :: laddress(1:2,1:maxblocks_alloc)
#else
      integer, allocatable, save :: surr_blks(:,:,:,:,:)
#ifdef SAVE_MORTS
      integer, allocatable, save :: surr_morts(:,:,:,:,:)
#endif
      integer, allocatable, save :: laddress(:,:)
#endif


! arrays to store info about block neighbors which are boundaries
      public :: bc_block_neighs,bc_block_neighs_send
      public :: bc_block_neighs_length
      public :: bc_block_neighs_status
      integer,save,allocatable :: bc_block_neighs(:,:)
      integer,save,allocatable :: bc_block_neighs_send(:,:)
      integer,save             :: bc_block_neighs_length
      integer,save             :: bc_block_neighs_status


! DECLARE variables which are targets
      target refine, derefine, newchild, empty
      target lrefine, nodetype, work_block
      target parent, coord, bsize, neigh
      target child, bnd_box, stay
!--------------------------------------------
!
! A convention is established for numbering the neighbors (or faces
! of a block. The first neighbor is at lower x coordinate, the 
! second at higher x, the third at lower y, fourth at higher y, fifth
! at lower z and the sixth at higher z.
!
! The convention by which the children of a block are numbered is the
! same as the fortran array ordering, so that the first child is
! at lower x, y and z coordinate, the second child is at upper x
! but lower y and z, the third is at lower x, upper y and lower z,
! and so on.
!
! When a block has a refined neighbor we will need to know which children
! of this neighbor are to provide guard cell information. The id's of the
! correct children are stored in kchild using the conventions described 
! above. For example, if we are working on the 3rd neighbor of the
! current block and it is at finer refinement level, then we must access
! the children designated by kchild(:,3), in this case children 1, 2, 5
! and 6.

!--------------------------------------------
!
!
! Description:
! This is the include file for a quad or oct-tree data structure,
! implemented on a parallel computer.
!
! The tree organizes a set of up to maxblocks_tr grids on each processor.
! All the grids are assumed to be cartesian with a uniform size. Each 
! grid has a level of refinement associated with it. The set of level 0
! grids cover the computational domain without overlap. Each grid
! can be the parent of 2**d offspring grids which completely cover 
! the sub-domain of their parents, where d is the physical dimension
! of the simulation. The linear resolution varies by a factor of 2 
! between successive levels of refinement. At no point do we allow the
! resolution to jump by more than one level of refinement.
!
!
! In the following list the index i ranges from 1 to maxblocks. 
!
!       neigh(2,mfaces,i)     local and processor ids of block i's neighbors,
!                               at i's refinement level. If a neighbor does 
!                               not exist both values are set to -1, unless 
!                               that face is at an external domain boundary
!                               where non-periodic boundary conditions are to
!                               be applied, in which case these are set to -20
!                               or less, depending on the boundary conditions
!                               to be applied on the boundary in question.
!       child(2,mchild,i)     local and processor ids of block i's children
!       type_of_child(mchild,i)  stores nodetype of children of block i
!       type_of_neigh(mfaces,i)  stores nodetype of neighbors of block i
!       parent(2,i)           local and processor ids of block i's parent
!       coord(mdim,i)         array storing x,y and z coordinates of the
!                               center of block i.
!       bnd_box(2,mdim,i)     bounding box information for block i. The 
!                               lower edge of block i along the j-th coordinate
!                               axis is at bnd_box(1,j,i) and the upper edge
!                               at bnd_box(2,j,i).
!       bsize(mdim,i)         size of block i in the x, y and z directions.
!       lrefine(i)            refinement level of block i.
!       nodetype(i)           defines the node type, if 1 then the node is a
!                               leaf node, if 2 then the node is a parent but
!                               with at least 1 leaf child, otherwise it is
!                               set to 3 and it does not have any up-to-date
!                               data.
!       empty(i)              used to designate empty blocks, for example
!                               when an obstacle is inserted inside the
!                               computational domain. normal blocks have
!                               empty=0, empty blocks have empty=1.
!
!       bflags(mflags,i)      an array of integer flags which can be used
!                               to control computation on the grid blocks
!                               and which are inherited by children from
!                               their parents.
!       
!       which_child(i)        an integer identifying which part of the
!                               parents volume this child corresponds to.
!       new_child(i)          if true then child has just been produced by
!                               a refinement step, otherwise false.
!       lnblocks              number of blocks on the local processor
!       new_lnblocks          the new number of blocks on the local 
!                               processor after a refinement or derefinement 
!                               step.
!       refine(i)             refinement flag. If set to .true. block i
!                               will be refined during the next call to
!                               REFINE_DEREFINE.
!       derefine(i)           derefinement flag. If set to .true. block i
!                               will be derefined during the next call to
!                               REFINE_DEREFINE, provided this blocks parent
!                               is not marked for refinement.
!       stay(i)               indicates block i refinement level unchanged.
!       work_block()          stores workload weighting associated with block i
!
!
!-----------------------------------------------------------------

      end module tree
