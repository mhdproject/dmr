!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_restrict_eblock_marker(rflag)




!------------------------------------------------------------------------
!
! This routine sets a logical flag to mark any leaf blocks and their
! parents where the leaf blocks border on coarser empty leaf blocks. 
! This flag is used inside the restriction routine to selectively
! apply the restriction operator before applying the routine which
! patches up communications in the neighborhood of empty interior
! blocks.
!
! Written :     Peter MacNeice          September 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace


#include "amr_shmem.fh"

      logical rflag(maxblocks)

!------------------------------------

      logical cflag(mchild),lflag
      integer remote_pe,remote_block
      integer cempty
      save cempty,lflag,remote_pe,remote_block

!------------------------------------


      if(lnblocks.gt.0) then



! cycle through the grid blocks on this processor
      do isg = 1,lnblocks


! Is this a parent of a leaf block ?
      if(nodetype(isg).eq.2) then


! Loop over the faces of this block.
       do jf=1,nfaces
       cflag(:) = .false. 

       remote_block = neigh(1,jf,isg)
       remote_pe    = neigh(2,jf,isg)

! If this neighbor exists
       if(remote_block.gt.0) then

       cempty = 0
       call shmem_integer_get ( cempty,
     .       empty(remote_block),1,remote_pe)

! If this neighbor is empty then flag this parent block and record
! which children also need to be flagged.
       if(cempty.eq.1) then
       rflag(isg)=.true.
       if(jf.eq.1) then
       cflag(1)=.true. 
       if(ndim.ge.2) cflag(3)=.true. 
       if(ndim.eq.3) then
       cflag(5)=.true. 
       cflag(7)=.true. 
       endif
       elseif(jf.eq.2) then
       cflag(2)=.true. 
       if(ndim.ge.2) cflag(4)=.true. 
       if(ndim.eq.3) then
       cflag(6)=.true. 
       cflag(8)=.true. 
       endif
       elseif(jf.eq.3) then
       cflag(1)=.true. 
       cflag(2)=.true. 
       if(ndim.eq.3) then
       cflag(5)=.true. 
       cflag(6)=.true. 
       endif
       elseif(jf.eq.4) then
       cflag(3)=.true. 
       cflag(4)=.true. 
       if(ndim.eq.3) then
       cflag(7)=.true. 
       cflag(8)=.true. 
       endif
       elseif(jf.eq.5) then
       cflag(1)=.true. 
       cflag(2)=.true. 
       cflag(3)=.true. 
       cflag(4)=.true. 
       elseif(jf.eq.6) then
       cflag(5)=.true. 
       cflag(6)=.true. 
       cflag(7)=.true. 
       cflag(8)=.true. 
       endif
       endif

       endif

       enddo

! If this parent has been flagged, then loop over its children
! and flag those which were identified in the previous section.
       lflag = .true.
       if(rflag(isg)) then
       do ich = 1,nchild

       if(cflag(ich)) then

       remote_block = child(1,ich,isg)
       remote_pe    = child(2,ich,isg)
       call shmem_logical_put(rflag(remote_block),
     .       lflag,1,remote_pe)

       endif
       enddo

       endif



      endif
      enddo
      endif

      call shmem_barrier_all()

      return
      end subroutine amr_restrict_eblock_marker
