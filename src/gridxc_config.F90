#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!< Contains global (MPI-related) data and exports an
! initialization routine.

module gridxc_config
!
!  Contains global data for internal operation
!  of the library.
!  This makes it non-thread-safe for now
!
! MPI variables
#ifdef HAVE_MPI
  integer,             save :: gridxc_comm
  integer,public,      save :: gridxc_myNode= -1, gridxc_totNodes=-1
#else
  integer,public,      save :: gridxc_myNode= 0, gridxc_totNodes=1
#endif
!
CONTAINS

#ifdef HAVE_MPI
!< Initialization routine. The `comm` argument is not present in the
! serial version.  
SUBROUTINE gridxc_init(comm)

!> MPI communicator for the library
integer, intent(in)           :: comm   ! NOT optional

integer :: mpierr
#else
!< Initialization routine. Basically a no-op for the serial version.
! The MPI version has a `comm` argument.
SUBROUTINE gridxc_init()
#endif
  
#ifdef HAVE_MPI
   gridxc_comm = comm
   call MPI_Comm_Size(comm,gridxc_totNodes,mpierr)
   call MPI_Comm_Rank(comm,gridxc_myNode,mpierr)
   !
#else
   gridxc_myNode = 0
   gridxc_totNodes = 1
#endif

end subroutine gridxc_init

end module gridxc_config

