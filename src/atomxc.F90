#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBXC
#include "xc_version.h"
#endif

MODULE gridxc_atom

  implicit none

  PUBLIC:: atomXC  ! Exchange and correlation of a spherical density

  PRIVATE ! Nothing is declared public beyond this point

CONTAINS

!< Finds total exchange-correlation energy and potential for a
! spherical electron density distribution.
!```  
! This version implements the Local (spin) Density Approximation and
! the Generalized-Gradient-Aproximation with the 'explicit mesh 
! functional' approach of White & Bird, PRB 50, 4954 (1994).
! Gradients are 'defined' by numerical derivatives, using 2*nn+1 mesh
!   points, where nn is a parameter defined below
! Ref: L.C.Balbas et al, PRB 64, 165110 (2001)
! Written by J.M.Soler using algorithms developed by 
!   L.C.Balbas, J.L.Martins and J.M.Soler, Dec.1996
! Van der Waals functional added by J.M.Soler, Jul.2008, as explained in
!   G.Roman-Perez and J.M.Soler, PRL 103, 096102 (2009)
! ************************* INPUT ***********************************
! INTEGER irel         : Relativistic exchange? (0=>no, 1=>yes)
! INTEGER nr           : Number of radial mesh points
! INTEGER maxr         : Physical first dimension of Dens and Vxc
! REAL*8  rmesh(nr)    : Radial mesh points. Must be nr.le.maxr
! INTEGER nSpin        : nSpin=1 => unpolarized; nSpin=2 => polarized
! REAL*8  Dens(maxr,nSpin) : Total (nSpin=1) or spin (nSpin=2) electron
!                            density at mesh points
! ************************* OUTPUT **********************************
! REAL*8  Ex              : Total exchange energy
! REAL*8  Ec              : Total correlation energy
! REAL*8  Dx              : IntegralOf( rho * (eps_x - v_x) )
! REAL*8  Dc              : IntegralOf( rho * (eps_c - v_c) )
! REAL*8  Vxc(maxr,nSpin) : (Spin) exch-corr potential
! ************************ UNITS ************************************
! Distances in atomic units (Bohr).
! Densities in atomic units (electrons/Bohr**3)
! Energy unit depending of parameter Eunit below (currently Rydberg)
!```

subroutine atomXC( irel, nr, maxr, rmesh, nSpin, Dens, Ex, Ec, Dx, Dc, Vxc )

! Module routines
  use gridxc_alloc,   only: alloc_default ! Sets (re)allocation defaults
  use gridxc_alloc,   only: de_alloc      ! Deallocates arrays
  use gridxc_mesh1D,  only: derivative    ! Performs numerical derivative
  use gridxc_sys,     only: die           ! Termination routine
  use gridxc_xcmod,   only: getXC         ! Returns the XC functional to be used
  use gridxc_gga,     only: ggaxc         ! General GGA XC routine
  use gridxc_mesh1D,  only: interpolation=>interpolation_local ! Interpolation routine
  use gridxc_lda,     only: ldaxc         ! General LDA XC routine
!  use gridxc_filter,only: kcPhi         ! Finds planewave cutoff of a radial func.
  use gridxc_radfft,  only: radfft        ! Radial fast Fourier transform
  use gridxc_alloc,   only: re_alloc      ! Reallocates arrays
  use gridxc_mesh1D,  only: set_mesh      ! Sets a one-dimensional mesh
  use gridxc_mesh1D,  only: set_interpolation  ! Sets the interpolation method
  use gridxc_vdwxc, only: vdw_decusp    ! Cusp correction to VDW energy
  use gridxc_vdwxc, only: vdw_localxc   ! Local LDA/GGA xc apropriate for vdW flavour
  use gridxc_vdwxc, only: vdw_get_qmesh ! Returns q-mesh for VDW integrals
  use gridxc_vdwxc, only: vdw_phi       ! Returns VDW functional kernel
  use gridxc_vdwxc, only: vdw_set_kcut  ! Fixes k-cutoff in VDW integrals
  use gridxc_vdwxc, only: vdw_theta     ! Returns VDW theta function

! Module types and variables
  use gridxc_precision, only: dp          ! Double precision real type
  use gridxc_alloc,   only: allocDefaults ! Derived type for allocation defaults
  use gridxc_config, only: myNode => gridxc_myNode

#ifdef DEBUG_XC
!  use gridxc_vdwxc, only: qofrho        ! Returns q(rho,grad_rho)

  use gridxc_debugXC, only: udebug        ! Output file unit for debug info
  use gridxc_debugXC, only: setDebugOutputUnit  ! Sets udebug
#endif /* DEBUG_XC */

#ifdef HAVE_LIBXC
#if XC_MAJOR_VERSION >= 4
    use xc_f03_lib_m
#else
    use xc_f90_types_m
    use xc_f90_lib_m
#endif
#endif /* HAVE_LIBXC */

  implicit none

! Argument types and dimensions
  integer, intent(in) :: irel       ! Relativistic exchange? (0=>no, 1=>yes)
  integer, intent(in) :: maxr       ! First dimension of arrays dens and Vxc
  integer, intent(in) :: nr               ! Number of radial mesh points
  integer, intent(in) :: nSpin            ! Number of spin components
  real(dp),intent(in) :: Dens(maxr,nSpin) ! (spin) Density at mesh points
  real(dp),intent(in) :: rmesh(nr)        ! Radial mesh points
  real(dp),intent(out):: Ex               ! Total exchange energy 
  real(dp),intent(out):: Ec               ! Total correlation energy
  real(dp),intent(out):: Dx               ! IntegralOf(rho*(eps_x-v_x))
  real(dp),intent(out):: Dc               ! IntegralOf(rho*(eps_c-v_c))
  real(dp),intent(out):: Vxc(maxr,nSpin)  ! (spin) xc potential

! Subroutine name
  character(len=*),parameter :: myName = 'atomXC '
  character(len=*),parameter :: errHead = myName//'ERROR: '

! Fix energy unit:  Eunit=1.0 => Hartrees,
!                   Eunit=0.5 => Rydbergs,
!                   Eunit=0.03674903 => eV
  real(dp),  parameter   :: Eunit = 0.5_dp

! Fix order of the numerical derivatives: used radial points = 2*nn+1
  integer, parameter :: nn = 5

! Fix number of points for radial FFTs (VDW only)
  integer, parameter :: nk = 512

! Fix the interpolation method to change to the uniform FFT mesh (VDW only)
  character(len=*),parameter:: interp_method = 'lagrange' !('lagrange'|'spline')

! Fix kin. energy leak to find the planewave cutoff of the radial density (VDW)
!  real(dp), parameter :: EtolKc = 0.003_dp ! Ry

! Fix limits to the planewave cutoff of the radial density (VDW only)
  real(dp), parameter :: kcMin  = 20._dp   ! Bohr^-1
  real(dp), parameter :: kcMax  = 50._dp   ! Bohr^-1
  real(dp), parameter :: Dmin   = 1.e-9_dp ! Min density when estimating kc
  real(dp), parameter :: Dcut   = 1.e-9_dp ! Min density for nonzero Vxc

! Fix the maximum number of functionals to be combined
  integer, parameter :: maxFunc = 10

! Max number of spin components
  integer, parameter :: maxSpin = 4

! Local variables and arrays
  logical :: &
    GGA, GGAfunc, VDW, VDWfunc
  integer :: &
    ik, in, in1, in2, iq, ir, is, ix, jn, ndSpin, nf, nq, nXCfunc
  real(dp) :: & 
    dEcdD(maxSpin), dEcdGD(3,maxSpin), dEcuspdD(maxSpin), &
    dExdD(maxSpin), dExdGD(3,maxSpin), dEdDaux(maxSpin),  &
    dEcuspdGD(3,maxSpin), dVcdD(maxSpin,maxSpin), dVxdD(maxSpin,maxSpin)
  real(dp) :: & 
    dGdm(-nn:nn), d2ydx2, dk, dr, Dtot, &
    Eaux, Ecut, epsC, epsCusp, epsNL, epsX, f1, f2, &  
    k(0:nk), kc, kmax, pi, r(0:nk), rmax, x0, xm, xp, y0, ym, yp, &
    XCweightC(maxFunc), XCweightX(maxFunc)
  character(len=20):: &
    XCauth(maxFunc), XCfunc(maxFunc)

  real(dp), pointer :: &
    D(:,:)=>null(), dGDdD(:,:)=>null(), drdm(:)=>null(), dVol(:)=>null(), &
    GD(:,:,:)=>null()
  real(dp), pointer:: &
    dphidk(:,:)=>null(), dtdgd(:,:,:)=>null(), dtdd(:,:)=>null(), &
    phi(:,:)=>null(), tk(:,:)=>null(), tr(:,:)=>null(), &
    uk(:,:)=>null(), ur(:,:)=>null()
  type(allocDefaults):: &
    prevAllocDefaults
#ifdef DEBUG_XC
!  real(dp):: q, dqdrho, dqdgrho(3)
!  real(dp):: epsCtmp, dEcdDtmp(nSpin), dEcdGDtmp(3,nSpin)
!   integer:: fileUnit
!   logical:: fileOpened
!   character(len=32):: fileName
#endif /* DEBUG_XC */

#ifdef HAVE_LIBXC
#if XC_MAJOR_VERSION >= 4
      type(xc_f03_func_t), allocatable :: xc_func(:)
      type(xc_f03_func_info_t), allocatable :: xc_info(:)
#else
      type(xc_f90_pointer_t), allocatable :: xc_func(:), xc_info(:)
#endif
      logical, allocatable                :: is_libxc(:)
      integer :: xc_ispin, xc_id, iostat
#endif
  
#ifdef DEBUG_XC
  call setDebugOutputUnit(myNode)   ! Initialize udebug variable
#endif /* DEBUG_XC */

! Check dimension of arrays Dens and Vxc
  if (maxr<nr) call die(errHead//'maxr<nr')

! Find number of diagonal spin components
  ndSpin = min( nSpin, 2 )

! Get the functional(s) to be used
  call getXC( nXCfunc, XCfunc, XCauth, XCweightX, XCweightC )

  ! Set GGA and VDW switches
  GGA = .false.
  VDW = .false.
  do nf = 1,nXCfunc
    if ( XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
      VDW = .true.
      GGA = .true.
    else if ( XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
      GGA = .true.
    else if ( XCfunc(nf).ne.'LDA' .and. XCfunc(nf).ne.'lda' .and. &
              XCfunc(nf).ne.'LSD' .and. XCfunc(nf).ne.'lsd' ) then
      call die(errHead//'Unknown functional '//XCfunc(nf))
    endif
  enddo ! nf

  ! Set the possible libxc objects to avoid call overhead in grid
#ifdef HAVE_LIBXC
  allocate(is_libxc(nXCfunc))
  is_libxc(:) = .false.
  allocate(xc_func(nXCfunc), xc_info(nXCfunc))
#endif
  
  do nf = 1,nXCfunc
      if ( XCauth(nf)(1:5) == 'LIBXC') then
#ifdef HAVE_LIBXC
      read(XCauth(nf)(7:11),iostat=iostat,fmt=*) xc_id
      if (iostat /= 0) call die("Bad libxc code in " // XCauth(nf))
        IF (nspin == 1) THEN
          xc_ispin = XC_UNPOLARIZED
        ELSE
          xc_ispin = XC_POLARIZED
        ENDIF

#if XC_MAJOR_VERSION >= 4           
        if ((irel == 1) .and. (xc_id == 1)) then
           ! Change LDA_X to LDA_X_REL (532)
           ! This was introduced in libXC v4
           xc_id = 532
        endif
        call xc_f03_func_init(xc_func(nf), xc_id, xc_ispin)
        xc_info(nf) = xc_f03_func_get_info(xc_func(nf))
#else
        call xc_f90_func_init(xc_func(nf), xc_info(nf), xc_id, xc_ispin)
        if ((irel == 1) .and. (xc_id == 1)) then
          ! To set the relativistic flag we need to call this
          ! general function, and include the default 4/3 factor
          ! that gives standard exchange (implicit in the above
          ! initialization call)
           call xc_f90_lda_x_set_par(xc_func(nf), 4.0_dp/3.0_dp,   &
                                    XC_RELATIVISTIC, 0.0_dp)
        endif
#endif
        is_libxc(nf) = .true.
#else
        call die("Libxc not compiled in. Cannot handle "// trim(XCauth(nf)))
#endif /* HAVE_LIBXC */
     endif
  enddo ! nf

! Set routine name for allocations
!  call alloc_default( old=prevAllocDefaults, routine=myName )

! Allocate temporary arrays
  call re_alloc( D,       1,nSpin, 1,nr, myName//'D' )
  call re_alloc( dGDdD,    -nn,nn, 1,nr, myName//'dGDdD' )
  call re_alloc( drdm,             1,nr, myName//'drdm'  )
  call re_alloc( dVol,             1,nr, myName//'dVol'  )
  call re_alloc( GD, 1,3, 1,nSpin, 1,nr, myName//'GD' )

! Find some parameters for the FFT mesh
  pi = 4.0_dp * atan(1.0_dp)
  rmax = rmesh(nr)
  dr = rmax / nk
  kmax = pi / dr
  dk = kmax / nk

! Find density and gradient of density at mesh points
  do ir = 1,nr

    ! Find interval of neighbour points to calculate derivatives
    in1 = max(  1, ir-nn ) - ir
    in2 = min( nr, ir+nn ) - ir

    ! Find weights of numerical derivation, dGdm = d(dF(n)/dn)/dF_m, 
    ! from Lagrange interpolation formula:
    ! F(n) = Sum_m F_m * Prod_{j/=m} (n-j)/(m-j)
    ! dF(n)/dn = Sum_m F_m * Prod_{j/=m,j/=n} (n-j) / Prod_{j/=m} (m-j)
    do in = in1,in2
      if (in==0) then
        dGdm(in) = 0
        do jn = in1,in2
          if (jn/=0) dGdm(in) = dGdm(in) + 1._dp / (0 - jn)
        enddo
      else
        f1 = 1
        f2 = 1
        do jn = in1,in2
          if (jn/=in .and. jn/=0) f1 = f1 * (0  - jn)
          if (jn/=in)             f2 = f2 * (in - jn)
        enddo
        dGdm(in) = f1 / f2
      endif
    enddo

    ! Find dr(m)/dm
    drdm(ir) = 0.0_dp
    do in = in1,in2
      drdm(ir) = drdm(ir) + rmesh(ir+in) * dGdm(in)
    enddo

    ! Find differential of volume. Use trapezoidal integration rule
    dVol(ir) = 4.0_dp * pi * rmesh(ir)**2 * drdm(ir)
    if (ir==1 .or. ir==nr) dVol(ir) = 0.5_dp*dVol(ir)

    ! Find the weights for the derivative d(gradF(n))/d(F(m)), of
    ! the gradient at point n with respect to the value at point m
    ! dGDdD = d((dD(r)/dr)(n)/dD_m = d( (dD(n)/dn) / (dr(n)/dn) ) / dD_m
    !       = d(dD(n)/dn)/dD_m / (dr(n)/dn)
    if (GGA) then
      do in = in1,in2
        dGDdD(in,ir) = dGdm(in) / drdm(ir)
      enddo
    endif

    ! Find density and gradient of density at this point
    do is = 1,nSpin
      D(is,ir) = Dens(ir,is)
    enddo
    if (GGA) then
      do is = 1,nSpin
        GD(1:3,is,ir) = 0.0_dp
        do in = in1,in2
          GD(3,is,ir) = GD(3,is,ir) + dGDdD(in,ir) * Dens(ir+in,is)
        enddo
      enddo
    endif ! GGA

  end do ! ir

! Van der Waals initializations
  if (VDW) then

    ! Allocate VdW arrays
    call vdw_get_qmesh( nq )
    call re_alloc( dtdd,         1,nq, 1,nSpin, myName//'dtdd'  )
    call re_alloc( dtdgd,   1,3, 1,nq, 1,nSpin, myName//'dtdgd' )
    call re_alloc( dphidk, 1,nq, 1,nq,          myName//'dphidk')
    call re_alloc( phi,    1,nq, 1,nq,          myName//'phi'   )
    call re_alloc( tk,     0,nk, 1,nq,          myName//'tk'    )
    call re_alloc( tr,     1,nr, 1,nq,          myName//'tr'    )
    call re_alloc( uk,     0,nk, 1,nq,          myName//'uk'    )
    call re_alloc( ur,     1,nr, 1,nq,          myName//'ur'    )

    ! Find planewave cutoff of density
!    kc = kcPhi( 0, nr, rmesh(:), sum(dens(:,1:ndSpin),2), EtolKc )
!    kc = min( kc, kcmax )
#ifdef DEBUG_XC
!   write(udebug,'(a,f12.6)') myName//'kc =', kc
#endif /* DEBUG_XC */

    ! Estimate the planewave cutoff of the density
    Ecut = 0
    do ir = 2,nr-1
      Dtot = sum(dens(ir,1:ndSpin))
      if (Dtot < Dmin) cycle
      xm = rmesh(ir-1)
      x0 = rmesh(ir)
      xp = rmesh(ir+1)
      ym = rmesh(ir-1)*sum(dens(ir-1,1:ndSpin))
      y0 = rmesh(ir)  *sum(dens(ir  ,1:ndSpin))
      yp = rmesh(ir+1)*sum(dens(ir+1,1:ndSpin))
      d2ydx2 = ( (yp-y0)/(xp-x0) - (y0-ym)/(x0-xm) ) * 2 / (xp-xm)
      Ecut = max( Ecut, abs(d2ydx2/y0) )
    end do ! ir
    kc = sqrt(Ecut)
    kc = max( kc, kcMin )
    kc = min( kc, kcMax )
#ifdef DEBUG_XC
!   write(udebug,'(a,f12.6)') myName//'kc =', kc
#endif /* DEBUG_XC */

    ! Set mesh cutoff to filter VdW kernel
    call vdw_set_kcut( kc )

    ! Find expansion of theta(q(ir))
    do ir = 1,nr
      call vdw_theta( nSpin, D(:,ir), GD(:,:,ir), tr(ir,1:nq), dtdd, dtdgd )
    end do ! ir

    ! Find uniform meshes for FFTs
    forall(ir=0:nk) r(ir) = ir*dr
    forall(ik=0:nk) k(ik) = ik*dk

    ! Find theta_iq(r) in a uniform radial mesh
    call set_interpolation( interp_method )
    call set_mesh( nr, rmesh )  ! 'From' mesh of tr array
    do iq = 1,nq
      tk(0:nk,iq) = interpolation( nk+1, r(0:nk), nr, tr(1:nr,iq) )
    end do

#ifdef DEBUG_XC
!    ! Write density
!    open(1,file='d.out')
!    do ir = 1,nr
!      write(1,'(3f15.9)') rmesh(ir), D(:,ir)
!    end do
!    close(1)
!    ! Write q(rho,grad_rho)
!    open(1,file='q.out')
!    do ir = 1,nr
!      call qofrho( sum(d(:,ir)), sum(gd(:,:,ir),2), q, dqdrho, dqdgrho )
!      write(1,'(3f15.9)') rmesh(ir), q
!    end do
!    close(1)
!    ! Write theta
!    open(1,file='trmesh.out')
!    do ir = 1,nr
!      write(1,'(30f15.9)') rmesh(ir), tr(ir,1:nq)
!    end do
!    close(1)
!    open(1,file='tr.out')
!    do ir = 0,nk
!      write(1,'(30f15.9)') r(ir), tk(ir,1:nq)
!    end do
!    close(1)
#endif /* DEBUG_XC */

    ! Fourier-transform theta_iq(r) for each iq
    do iq = 1,nq
      call radfft( 0, nk, rmax, tk(0:nk,iq), tk(0:nk,iq) )
    end do ! iq

    ! Find u_iq(r) = Sum_iq' Int_dr' phi_iq_iq'(r,r')*theta_iq'(r')
    ! by convolution in reciprocal space
    do ik = 0,nk
      ! Find Fourier transform of VdW kernel phi_iq_iq'(r,r')
      call vdw_phi( k(ik), phi, dphidk )
      ! Find Fourier transform of u_iq(r)
      uk(ik,1:nq) = matmul( tk(ik,1:nq), phi(1:nq,1:nq) )
    end do ! ik

    ! Inverse Fourier transform of u_iq(r) for each iq
    do iq = 1,nq
      call radfft( 0, nk, kmax, uk(0:nk,iq), uk(0:nk,iq) )
    end do ! iq

    ! Find u_iq(r) in the original radial mesh
    call set_mesh( nk+1, xmin=0._dp, xmax=rmax )  ! 'From' mesh of uk array
    do iq = 1,nq
      ur(1:nr,iq) = interpolation( nr, rmesh(1:nr), nk+1, uk(0:nk,iq) )
    end do ! iq

#ifdef DEBUG_XC
!    ! Write u
!    open(1,file='ur.out')
!    do ir = 0,nk
!      write(1,'(30f15.9)') r(ir), uk(ir,1:nq)
!    end do
!    close(1)
!    open(1,file='urmesh.out')
!    do ir = 1,nr
!      write(1,'(30f15.9)') rmesh(ir), ur(ir,1:nq)
!    end do
!    close(1)
#endif /* DEBUG_XC */

  end if ! (VDW)

! Initialize output
  Ex = 0.0_dp
  Ec = 0.0_dp
  Dx = 0.0_dp
  Dc = 0.0_dp
  Vxc(1:nr,1:nSpin) = 0.0_dp

! Loop over mesh points
  do ir = 1,nr

    ! Find interval of neighbour points to calculate derivatives
    in1 = max(  1, ir-nn ) - ir
    in2 = min( nr, ir+nn ) - ir

    ! Loop over exchange-correlation functions
    do nf = 1,nXCfunc
      ! Is this a GGA or VDW?
      if (XCfunc(nf).eq.'VDW' .or. XCfunc(nf).eq.'vdw') then
        VDWfunc = .true.
        GGAfunc = .true.
      else if (XCfunc(nf).eq.'GGA' .or. XCfunc(nf).eq.'gga') then
        VDWfunc = .false.
        GGAfunc = .true.
      else
        VDWfunc = .false.
        GGAfunc = .false.
      endif

      ! Find exchange and correlation energy densities and their 
      ! derivatives with respect to density and density gradient
      if (VDWfunc) then

        ! Local exchange-corr. part from the apropriate LDA/GGA functional
        call vdw_localxc( irel, nSpin, D(:,ir), GD(:,:,ir), epsX, epsC, &
                          dExdD, dEcdD, dExdGD, dEcdGD )

#ifdef DEBUG_XC
!        ! Select only non local correlation energy and potential
!        epsX = 0
!        epsC = 0
!        dExdD = 0
!        dEcdD = 0
!        dExdGD = 0
!        dEcdGD = 0
#endif /* DEBUG_XC */

        ! Local cusp correction to nonlocal VdW energy integral
        call vdw_decusp( nSpin, D(:,ir), GD(:,:,ir), &
                         epsCusp, dEcuspdD, dEcuspdGD )

#ifdef DEBUG_XC
!        ! Select only non local correlation energy and potential
!        epsCusp = 0
!        dEcuspdD = 0
!        dEcuspdGD = 0
#endif /* DEBUG_XC */

        ! Find expansion of theta(q(r)) for VdW
        call vdw_theta( nSpin, D(:,ir), GD(:,:,ir), tr(ir,1:nq), dtdd, dtdgd )

        ! Add nonlocal VdW energy contribution and its derivatives
        Dtot = sum(D(1:ndSpin,ir))
        epsNL = epsCusp + 0.5_dp*sum(ur(ir,1:nq)*tr(ir,1:nq))/(Dtot+tiny(Dtot))
        epsC = epsC + epsNL
        do is = 1,nSpin
          dEcdD(is) = dEcdD(is) + dEcuspdD(is) &
                    + sum(ur(ir,1:nq)*dtdd(1:nq,is))
          do ix = 1,3
            dEcdGD(ix,is) = dEcdGD(ix,is) + dEcuspdGD(ix,is) &
                          + sum(ur(ir,1:nq)*dtdgd(ix,1:nq,is))
          end do ! ix
        end do ! is

      else if (GGAfunc) then
        call ggaxc( XCauth(nf), irel, nSpin, D(:,ir), GD(:,:,ir),  &
                    epsX, epsC, dExdD, dEcdD, dExdGD, dEcdGD       &
#ifdef HAVE_LIBXC
                    , is_libxc(nf), xc_func(nf), xc_info(nf) )
#else
                    )
#endif
                    
#ifdef DEBUG_XC
!        call ldaxc( 'PW92', irel, nSpin, D(:,ir), Eaux, epsC,  &
!                    dEdDaux, dEcdD, dVxdD, dVcdD )
!        call ggaxc( XCauth(nf), irel, nSpin, D(:,ir), GD(:,:,ir),  &
!                    epsX, epsCtmp, dExdD, dEcdDtmp, dExdGD, dEcdGDtmp )
#endif /* DEBUG_XC */
      else
        call ldaxc( XCauth(nf), irel, nSpin, D(:,ir), epsX, epsC, &
                    dExdD, dEcdD, dVxdD, dVcdD                    &
#ifdef HAVE_LIBXC
                    , is_libxc(nf), xc_func(nf), xc_info(nf) )
#else
                    )
#endif
      endif

      ! Scale terms by weights
      epsX = epsX*XCweightX(nf)
      epsC = epsC*XCweightC(nf)
      dExdD(1:nSpin) = dExdD(1:nSpin)*XCweightX(nf)
      dEcdD(1:nSpin) = dEcdD(1:nSpin)*XCweightC(nf)
      if (GGAfunc) then
        dExdGD(1:3,1:nSpin) = dExdGD(1:3,1:nSpin)*XCweightX(nf)
        dEcdGD(1:3,1:nSpin) = dEcdGD(1:3,1:nSpin)*XCweightC(nf)
      endif

      ! Add contributions to exchange-correlation energy and its
      ! derivatives with respect to density at all points
      do is = 1,nSpin
        Ex = Ex + dVol(ir)*D(is,ir)*epsX
        Ec = Ec + dVol(ir)*D(is,ir)*epsC
        Dx = Dx + dVol(ir)*D(is,ir)*(epsX - dExdD(is))
        Dc = Dc + dVol(ir)*D(is,ir)*(epsC - dEcdD(is))
        Vxc(ir,is) = Vxc(ir,is) + dVol(ir)*(dExdD(is) + dEcdD(is))
        if (GGAfunc) then
          do in = in1,in2
            Dx= Dx - dVol(ir)*D(is,ir+in)*dExdGD(3,is)*dGDdD(in,ir)
            Dc= Dc - dVol(ir)*D(is,ir+in)*dEcdGD(3,is)*dGDdD(in,ir)
            Vxc(ir+in,is) = Vxc(ir+in,is) + &
              dVol(ir)*(dExdGD(3,is) + dEcdGD(3,is))*dGDdD(in,ir)
          enddo ! in
        endif ! (GGAfunc)
      enddo ! is

    enddo ! nf

  enddo ! ir

#ifdef HAVE_LIBXC
  do nf = 1,nXCfunc
    if (is_libxc(nf)) then
#if XC_MAJOR_VERSION >= 4
      call xc_f03_func_end(xc_func(nf))
#else
      call xc_f90_func_end(xc_func(nf))
#endif
    endif
  enddo
  deallocate(xc_func,xc_info,is_libxc)
#endif

! Divide by volume element to obtain the potential (per electron)
  do is = 1,nSpin
    ! Avoid ir=1 => r=0 => dVol=0
    Vxc(2:nr,is) = Vxc(2:nr,is) / dVol(2:nr)
    ! Extrapolate to the origin from first two points, requiring dVxc/di=0
    Vxc(1,is) = (4*Vxc(2,is) - Vxc(3,is)) / 3
  enddo ! is

! Make Vxc=0 if VDWfunctl and Dens<Dcut, to avoid singularities
  if (VDW) then
    do ir = 1,nr
      Dtot = sum(Dens(ir,1:ndSpin))
      if (Dtot<Dcut) Vxc(ir,:) = 0
    end do
  end if ! (VDWfunctl)

! Divide by energy unit
  Ex = Ex / Eunit
  Ec = Ec / Eunit
  Dx = Dx / Eunit
  Dc = Dc / Eunit
  Vxc(1:nr,1:nSpin) = Vxc(1:nr,1:nSpin) / Eunit

#ifdef DEBUG_XC
!    ! Write potential
!    open(1,file='v.out')
!    do ir = 1,nr
!      write(1,'(3f15.9)') rmesh(ir), Vxc(ir,1:nSpin)
!    end do
!    close(1)
#endif /* DEBUG_XC */

! Deallocate VDW-related arrays
  if (VDW) then
    call de_alloc( ur,       myName//'ur' )
    call de_alloc( uk,       myName//'uk' )
    call de_alloc( tr,       myName//'tr' )
    call de_alloc( tk,       myName//'tk' )
    call de_alloc( phi,      myName//'phi' )
    call de_alloc( dphidk,   myName//'dphidk' )
    call de_alloc( dtdgd,    myName//'dtdgd' )
    call de_alloc( dtdd,     myName//'dtdd' )
  endif ! (VDW)

! Deallocate temporary arrays
  call de_alloc( GD,    myName//'GD' )
  call de_alloc( dVol,  myName//'dVol' )
  call de_alloc( drdm,  myName//'drdm' )
  call de_alloc( dGDdD, myName//'dGDdD' )
  call de_alloc( D,     myName//'D' )

! Restore previous allocation defaults
!  call alloc_default( restore=prevAllocDefaults )

#ifdef DEBUG_XC
!  fileUnit = 57
!  fileName = 'atomxc.vxc'
!  inquire(file=fileName,opened=fileOpened)
!  if (.not.fileOpened) open(unit=fileUnit,file=fileName)
!  write(fileUnit,*) nr, nSpin
!  do ir = 1,nr
!    write(fileUnit,'(5e15.6)') &
!      rmesh(ir), Dens(ir,1:nSpin), Vxc(ir,1:nSpin)
!  end do
#endif /* DEBUG_XC */

end subroutine atomXC

END MODULE gridxc_atom


