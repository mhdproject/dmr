!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------------
! mpi_morton module



      module mpi_morton

      use paramesh_dimensions

      private

      public :: npts_neigh
      integer, parameter :: npts_neigh = 500000

! 
! variables for storing the morton environment
      public :: pe_remote,pe_source,pe_destination
      public :: mortonbnd,r_mortonbnd
      public :: morton_limits
      public :: mortonenv,ir_buf,is_buf
      public :: no_of_mortonbnds_received
      public :: morton_limits_set
      public :: mpi_tree_set
#ifndef LIBRARY
      integer, save :: 
     .    mortonbnd(6,1:3,1:maxblocks)
#else
      integer, allocatable, save :: 
     .    mortonbnd(:,:,:)
#endif
      integer, save :: no_of_mortonbnds_received
      logical, save :: morton_limits_set
      logical, save :: mpi_tree_set


      integer, save,dimension(:),allocatable :: pe_remote
      integer, save,dimension(:),allocatable :: pe_source
      integer, save,dimension(:),allocatable :: pe_destination
      integer, 
     .  save,dimension(:,:,:,:),allocatable :: r_mortonbnd
      integer, 
     .  save,dimension(:,:,:,:),allocatable :: morton_limits
      integer, save,dimension(:,:),allocatable :: ir_buf
      integer, save,dimension(:,:),allocatable :: is_buf

      public :: commatrix_send, commatrix_recv
      public :: commatrix_guard,commatrix_prol
      public :: commatrix_flux
      public :: commatrix_restrict
      public :: laddress_guard,laddress_prol,laddress_flux
      public :: laddress_restrict
      integer, save,dimension(:),allocatable :: commatrix_send
      integer, save,dimension(:),allocatable :: commatrix_recv
      integer, save,dimension(:,:),allocatable :: commatrix_guard
      integer, save,dimension(:,:),allocatable :: commatrix_prol
      integer, save,dimension(:,:),allocatable :: commatrix_flux
      integer, save,dimension(:,:),allocatable :: 
     .                                         commatrix_restrict
#ifndef LIBRARY
      integer, save,dimension(1:2,1:maxblocks_alloc) 
     .                                          :: laddress_guard
      integer, save,dimension(1:2,1:maxblocks_alloc) 
     .                                          :: laddress_prol
      integer, save,dimension(1:2,1:maxblocks_alloc) 
     .                                          :: laddress_flux
      integer, save,dimension(1:2,1:maxblocks_alloc) 
     .                                          :: laddress_restrict
#else
      integer, allocatable, save,dimension(:,:) 
     .                                          :: laddress_guard
      integer, allocatable, save,dimension(:,:) 
     .                                          :: laddress_prol
      integer, allocatable, save,dimension(:,:) 
     .                                          :: laddress_flux
      integer, allocatable, save,dimension(:,:) 
     .                                          :: laddress_restrict
#endif

! list of block edges which need diagonal info during edge averaging
      public :: edge_mark,no_of_diagonal_edges
      integer, save :: edge_mark(6,4,npts_neigh)
      integer, save :: no_of_diagonal_edges

! a list of blocks to be sent from the local processor
      public :: to_be_sent,to_be_sent_guard,to_be_sent_prol
      public :: to_be_sent_flux
      integer,save,dimension(:,:,:),allocatable :: to_be_sent
      integer,save,dimension(:,:,:),allocatable 
     .                                   :: to_be_sent_guard
      integer,save,dimension(:,:,:),allocatable 
     .                                   :: to_be_sent_prol
      integer,save,dimension(:,:,:),allocatable 
     .                                   :: to_be_sent_flux

      public :: to_be_sent_restrict
      integer,dimension(:,:,:),allocatable 
     .                                   :: to_be_sent_restrict

! a list of blocks to be received by the local processor
      public :: to_be_received
      integer,save,dimension(:,:,:),allocatable :: to_be_received

! used to make searching of laddress more efficient
      public :: ladd_strt,ladd_end
      integer,save,dimension(:),allocatable :: ladd_strt,ladd_end

!new code
      public :: to_be_received_guard
      public :: to_be_received_prol
      public :: to_be_received_flux
      integer,save,dimension(:,:,:),allocatable 
     .                                   :: to_be_received_guard
      integer,save,dimension(:,:,:),allocatable 
     .                                   :: to_be_received_prol
      integer,save,dimension(:,:,:),allocatable 
     .                                   :: to_be_received_flux

      public :: to_be_received_restrict
      integer,dimension(:,:,:),allocatable 
     .                                   :: to_be_received_restrict

      public :: pe_source_guard
      public :: pe_source_prol
      public :: pe_source_flux
      public :: pe_source_restrict
      integer, save,dimension(:),allocatable :: pe_source_guard
      integer, save,dimension(:),allocatable :: pe_source_prol
      integer, save,dimension(:),allocatable :: pe_source_flux
      integer, save,dimension(:),allocatable 
     .                                   :: pe_source_restrict

      public :: message_size_cc
      public :: message_size_fc
      public :: message_size_ec
      public :: message_size_nc
      public :: message_size_wk
      public :: mess_segment_loc
      integer,save,dimension(2*27) :: message_size_cc,message_size_fc
      integer,save,dimension(2*27) :: message_size_ec,message_size_nc
      integer,save,dimension(2*27) :: message_size_wk
      integer,save,dimension(:),allocatable :: mess_segment_loc

      public :: temprecv_buf
!     public :: send_buf, recv_buf
!     public :: temptrecv_buf
      real, save,dimension(:),allocatable :: temprecv_buf
!     real, save,dimension(:),allocatable :: send_buf,recv_buf
!     real, save,dimension(:),allocatable :: temptrecv_buf

      public :: l_datapacked
      logical,save,dimension(5) :: l_datapacked

!new code end

      public :: largest_no_of_blocks,largest_no_of_blocks_guard
      public :: largest_no_of_blocks_prol,largest_no_of_blocks_flux
      public :: largest_no_of_blocks_restrict
      public :: max_no_to_send,max_no_to_send_guard
      public :: max_no_to_send_prol,max_no_to_send_flux
      public :: max_no_to_send_restrict
      public :: strt_guard,strt_prol,strt_flux
      public :: strt_restrict,no_commatrix_guard
      public :: no_commatrix_prol,no_commatrix_flux
      public :: no_commatrix_restrict
      integer,save  :: largest_no_of_blocks
      integer,save  :: largest_no_of_blocks_guard
      integer,save  :: largest_no_of_blocks_prol
      integer,save  :: largest_no_of_blocks_flux
      integer,save  :: largest_no_of_blocks_restrict
      integer,save  :: max_no_to_send
      integer,save  :: max_no_to_send_guard
      integer,save  :: max_no_to_send_prol
      integer,save  :: max_no_to_send_flux
      integer,save  :: max_no_to_send_restrict
      integer,save  :: strt_guard
      integer,save  :: strt_prol
      integer,save  :: strt_flux
      integer,save  :: strt_restrict
      integer,save  :: no_commatrix_guard
      integer,save  :: no_commatrix_prol
      integer,save  :: no_commatrix_flux
      integer,save  :: no_commatrix_restrict

      public :: lperiodicx,lperiodicy,lperiodicz
      logical, save :: lperiodicx
      logical, save :: lperiodicy
      logical, save :: lperiodicz

      public :: treeinfo

      type treeinfo
        real coord(3)
        real bsize(3)
        real bnd_box(2,3)
        integer parent(2)
        integer which_child
        logical newchild
        integer neigh(2,6)
        integer lrefine
        integer nodetype
        integer empty
      end type treeinfo


      end module mpi_morton
!-----------------------------------------------------------------------
