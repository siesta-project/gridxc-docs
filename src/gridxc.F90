#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!<
!```
! Provides the following main XC routines:
!
!   gridxc_init     ! Initialization
!
!   gridxc_atomXC   ! XC for a spherical charge distribution
!   gridxc_cellXC   ! XC for a periodic unit cell
!   gridxc_setXC    ! Sets XC functional(s) to be used by atomXC and/or cellXC
!
!   (For the latter there are two simpler alternatives:
!   gridxc_setxc_family_authors(family,authors)
!   gridxc_setxc_libxc_ids(nfuncs,libxc_ids)
!   )
!
!   gridxc_getXC    ! Returns the XC functional(s) being used
!
! Real kinds (precision) of arguments to call atomxc and cellxc
!   grid_p ! Precision for grid arrays to call cellxc

! Secondary entry points for testers and lower-level programming
!   gridxc_ldaxc    ! LDA-XC functionals
!   gridxc_ggaxc    ! GGA-XC functionals

! Extra utilities placed here for non-siesta users
! See correspondig modules in the source for usage documentation
!
!   gridxc_nfft          ! Get allowed sizes for FFTs (for VDW functionals)
!   setMeshDistr         ! Set a distribution of mesh points over processors
!   myMeshBox            ! Get my processor mesh box
!
!   setDebugOutputUnit   ! Initialize debug report
!   closeDebugOutputFile ! Print debug report
!```

MODULE gridXC

! Real kinds (precision) of arguments

  USE precision, only: dp      ! Standard real-kind (double) precision
  USE precision, only: grid_p  ! Precision for grid arrays

! Main entry routines of gridXC library
!--------------------------------------------------------------
  use gridxc_config, only: gridxc_init  ! initialization

  ! XC for a spherical charge distribution
  USE m_atomXC, only: gridxc_atomXC => atomXC 

  USE m_cellXC, only: gridxc_cellXC => cellXC  ! XC for a periodic unit cell
  USE xcmod,    only: gridxc_getXC => getXC   ! Returns XC functional(s)
  USE xcmod,    only: gridxc_setXC => setXC   ! Sets XC functional(s)
#ifdef HAVE_LIBXC
  USE xcmod,    only: gridxc_setXC_libxc => setXC_libxc_ids   ! Sets XC functional(s)
#endif

! Secondary entry points for testers and lower-level programming
  USE m_ldaxc,  only: gridxc_ldaxc => ldaxc    ! LDA-XC functionals
  USE m_ggaxc,  only: gridxc_ggaxc => ggaxc    ! GGA-XC functionals

! Extra utilities placed here for convenience
! See correspondig modules for usage documentation
  USE m_fft_gpfa,    only: gridxc_nfft => nfft  ! Get allowed sizes for FFTs
!----------------------------------------------------------------------
!-----------------------------------
#ifdef DEBUG_XC
  USE debugXC,  only: setDebugOutputUnit   ! Set debug report
  USE debugXC,  only: closeDebugOutputFile ! Print debug report
#endif
  USE mesh3d,   only: myMeshBox            ! Get my processor mesh box
  USE mesh3d,   only: setMeshDistr         ! Set a distribution of mesh
                                           ! points over parallel processors
  PUBLIC

END MODULE gridXC
  
