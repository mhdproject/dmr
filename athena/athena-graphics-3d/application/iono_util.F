#define DEBUG1
#define DEBUG3
#define SETUP

c#define DEBUGX


      module hashing
c
      integer emptyloc,primenum
      parameter(emptyloc=-1,primenum=137)
c     "man factor"
c
c Also primenum must be larger than the possible number of nodes on a
c process.
      integer ht(0:primenum-1)
      save

      contains

      subroutine hash_table_init
      implicit none
      ht = emptyloc
      return
      end subroutine

      subroutine hash_table_write
      do i=0,primenum-1
         write(64,'(2i20)') i,ht(i)
      enddo
      return
      end subroutine
      
      subroutine hash_linsrch(key,loc,new_entry)
      implicit none
      integer, intent(in)    :: key
      integer, intent(out)   :: loc
      logical, intent(out)   :: new_entry
      integer i
c
c search the hash table ht(0:nb-1),
c where each bucket has exactly 1 slot,      
c using linear probing.  If HT(loc) = EMPTYLOC, then the "loc"-th bucket
c is emptyloc and KEY can be entered into the table. 
c
c      i = hash_function(key)  where i \in [0 .. (nb-1)]
      i = mod(key,primenum)
      loc = i
      do while ((ht(loc) .ne. key) .and. (ht(loc) .ne. EMPTYLOC))
c     treat the table as circular
         print *,' hashing collision'
         loc = mod(loc+1,primenum)
         if( i .eq. loc ) then 
            print *,' hash table full'
            stop 2118
         endif
      enddo
      if (ht(loc) .eq. emptyloc) then
c        insert a new entry into the hash table
         ht(loc) = key
         write(63,'(2i20)') loc,key
         new_entry = .true.
      else
         new_entry = .false.
      endif
      return
      end subroutine
      end module hashing




      module ionospherefh
      implicit none
      integer mbpoin
      parameter (mbpoin = 300)
      integer mbface
      parameter (mbface = 2*mbpoin)
      real * 8 rearth,rionosphere
      parameter(rearth=6.37,rionosphere=rearth+0.1)

      integer nbface,nbpoin

      integer intmab2d(3,mbface)
c      real * 8 sigmah(mbpoin),sigmap(mbpoin)
c      real * 8 coord2d(3,mbpoin)
c      real * 8 areap(mbpoin)
c      real * 8 geom2d(13,mbface)
c      common /ionosphere/ 
c     $     nbpoin,nbface,
c     $     intmab2d,sigmah,sigmap,
c     $     coord2d,geom2d,areap

      real*8 sigmah(mbpoin),sigmap(mbpoin)
      real*8 coord2d(3,mbpoin)
      real*8 areap(mbpoin)
      real*8 geom2d(13,mbface)
      save
      end module

      module iono_util_data
      use ionospherefh
      implicit none
c#include "tree.fh"
c#include "ionosphere.fh"

      integer ndatas,ndatag
      parameter ( ndatas = 1, ndatag = 1 )

! set data used to define global node ids
      
! maximum no. of elements on a block, on a processor, in final global list
! and in the interim list with redundancy.
      integer max_elem_b,max_elem_l,max_elem_f
c      parameter ( max_elem_b = (nxb*nyb+nxb*nzb+nyb*nzb)*2*2 )
      parameter ( max_elem_b = 3*64*2*2 )
      parameter ( max_elem_l = max_elem_b*5 )

! maximum no. of surface nodes on a block, on a processor, in final global list
! and in the interim list with redundancy.
      integer max_nodes_l
      parameter ( max_nodes_l = max_elem_l*2 )
      
      SAVE
      end module iono_util_data

      
      module iono_data_pass1
      use iono_util_data
      implicit none
      integer nbpoin_l,nbface_l
      integer nodeix_l(max_nodes_l),nodeiy_l(max_nodes_l)
      integer nodeiz_l(max_nodes_l),nodeib_l(max_nodes_l)
      integer ispoin_l(max_nodes_l)
      integer intmat_k(3,max_elem_l)
      integer intmat_l(3,max_elem_l)
      real coord2d_l(3,max_nodes_l)
      real coord2del(3,3,max_elem_l)
      save
      end module iono_data_pass1
      
      module iono_data_pass2
      use iono_util_data
      implicit none
      integer nbpoing
      integer nodeixg(mbpoin),nodeiyg(mbpoin)
      integer nodeizg(mbpoin),nodeibg(mbpoin)
      integer nodeipg(mbpoin)
      integer ispoing(mbpoin)
      integer nbfaceg
      integer intmatg(3,mbface)
      real coord2dg(3,mbpoin)
      save 
      end module iono_data_pass2

      module iono_data_pass3
      use iono_util_data
      implicit none
      integer nbpoinf
      integer nodeixf(mbpoin),nodeiyf(mbpoin)
      integer nodeizf(mbpoin),nodeibf(mbpoin)
      integer nodeipf(mbpoin)
c      integer ispoinf(mbpoin)
      integer ispoinf(3*mbface)

      integer nbfacef
      integer intmatf(3,mbface)
      real coord2df(3,mbpoin)
      save
      end module iono_data_pass3
