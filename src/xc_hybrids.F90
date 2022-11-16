#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!@LICENSE
!
!     MODULE xc_hybrids
!
! Copyright Honghui Shang, 2010
!           Xinming Qin, 2013
!           Yann Pouillon, 2017, 2018, 2019
!
! This module implements the semi-local part of the HSE06 and PBE0
! exchange-correlation functionals.
!
! Standalone routines written by Honghui Shang, April 01, 2010
! Edited by Xinming Qin, December 12, 2013
! Ported to SIESTA 4.x by Yann Pouillon, 2017-2018
! Refactored into a Fortran module by Yann Pouillon, September, 2019
!-----------------------------------------------------------------------
MODULE gridxc_hybrids

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hsexc  ! HSE06 XC functional
  PUBLIC :: pbe0xc ! PBE0 XC functional

CONTAINS

  ! ********************************************************************
  ! Implements HSE Generalized-Gradient-Approximation.
  ! Ref: JOURNAL OF CHEMICAL PHYSICS VOLUME 120, NUMBER 16 (2004)
  ! edited from  L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
  ! Coded by shanghui, 2010.04.01
  ! Modified by xmqin, 2013,09.12
  ! Set omega as an input parameter read from .fdf
  ! ******** INPUT *****************************************************
  ! INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
  ! INTEGER nspin          : Number of spin polarizations (1 or 2)
  ! REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
  !                           spin electron density (if nspin=2)
  ! REAL*8  GDens(3,nspin) : Total or spin density gradient
  ! ******** OUTPUT ****************************************************
  ! REAL*8  EX             : Exchange energy density
  ! REAL*8  EC             : Correlation energy density
  ! REAL*8  DEXDD(nspin)   : Partial derivative
  !                           d(DensTot*Ex)/dDens(ispin),
  !                           where DensTot = Sum_ispin( Dens(ispin) )
  !                          For a constant density, this is the
  !                          exchange potential
  ! REAL*8  DECDD(nspin)   : Partial derivative
  !                           d(DensTot*Ec)/dDens(ispin),
  !                           where DensTot = Sum_ispin( Dens(ispin) )
  !                          For a constant density, this is the
  !                          correlation potential
  ! REAL*8  DEXDGD(3,nspin): Partial derivative
  !                           d(DensTot*Ex)/d(GradDens(i,ispin))
  ! REAL*8  DECDGD(3,nspin): Partial derivative
  !                           d(DensTot*Ec)/d(GradDens(i,ispin))
  ! ********* UNITS ****************************************************
  ! Lengths in Bohr
  ! Densities in electrons per Bohr**3
  ! Energies in Hartrees
  ! Gradient vectors in cartesian coordinates
  ! ********* ROUTINES CALLED ******************************************
  ! EXCHNG, PW92C
  ! ********************************************************************
  SUBROUTINE HSEXC(Irel, Nspin, Dens, Gdens, Ex, Ec, Dexdd, Decdd, &
               Dexdgd, Decdgd)
 
    USE gridxc_precision, ONLY: dp
    USE gridxc_lda, ONLY: EXCHNG, PW92C
    USE gridxc_xwpbe, ONLY: XWPBE

    IMPLICIT NONE
 
    INTEGER, INTENT(IN) :: Irel, Nspin
    REAL(dp), INTENT(IN) :: Dens(Nspin), Gdens(3,Nspin)
    REAL(dp), INTENT(INOUT) :: Decdd(Nspin), Decdgd(3,Nspin), &
      Dexdd(Nspin), Dexdgd(3,Nspin), Ec, Ex
 
    ! Internal variables
    INTEGER :: is, ix
 
    REAL(dp) :: a, beta, d(2), dadd, decudd, df1dd, &
      df2dd, df3dd, df4dd, df1dgd, df3dgd, df4dgd, &
      dfcdd(2), dfcdgd(3,2), dfdd, dfdgd, dfxdd(2), &
      dfxdgd(3,2), dhdd, dhdgd, dkfdd, dksdd, dpdd, &
      dpdz, drsdd, ds(2), dsdd, dsdgd, dt, dtdd, dtdgd, &
      dzdd(2), ecunif, exunif, f, f1, f2, f3, &
      f4, fc, fx, gamma, gd(3,2), gdm(2), &
      gdms, gdmt, gds, gdt(3), h, kappa, &
      kf, kfs, ks, mu, phi, pi, rs, s, t, &
      vcunif(2), vxunif(2), zeta
    REAL(dp) :: d_left, kfs_left, s_left, fsr_left, d_right, &
      kfs_right, s_right, fsr_right
    REAL(dp) :: fx_pbe, fx_cp2k, dfxdd_cp2k(2), dfxdgd_cp2k(3,2)
 
    ! Lower bounds of density and its gradient to avoid divisions by zero
    REAL(dp), PARAMETER :: DENMIN=1.D-12
    REAL(dp), PARAMETER :: GDMIN=1.D-12
 
    ! Fix some numerical parameters
    REAL(dp), PARAMETER :: FOUTHD=4.D0/3.D0
    REAL(dp), PARAMETER :: HALF=0.5D0
    REAL(dp), PARAMETER :: THD=1.D0/3.D0
    REAL(dp), PARAMETER :: THRHLF=1.5D0
    REAL(dp), PARAMETER :: TWO=2.D0
    REAL(dp), PARAMETER :: TWOTHD=2.D0/3.D0
 
    ! Fix some more numerical constants
    pi = 4*ATAN(1.D0)
    beta = 0.066725D0
    gamma = (1-LOG(TWO))/pi**2
    mu = beta*pi**2/3
    kappa = 0.804D0
 
    ! Translate density and its gradient to new variables
    IF ( Nspin==1 ) THEN
       d(1) = HALF*Dens(1)
       d(2) = d(1)
       dt = MAX(DENMIN,Dens(1))
       DO ix = 1, 3
          gd(ix,1) = HALF*Gdens(ix,1)
          gd(ix,2) = gd(ix,1)
          gdt(ix) = Gdens(ix,1)
       ENDDO
    ELSE
       d(1) = Dens(1)
       d(2) = Dens(2)
       dt = MAX(DENMIN,Dens(1)+Dens(2))
       DO ix = 1, 3
          gd(ix,1) = Gdens(ix,1)
          gd(ix,2) = Gdens(ix,2)
          gdt(ix) = Gdens(ix,1) + Gdens(ix,2)
       ENDDO
    ENDIF
    gdm(1) = SQRT(gd(1,1)**2+gd(2,1)**2+gd(3,1)**2)
    gdm(2) = SQRT(gd(1,2)**2+gd(2,2)**2+gd(3,2)**2)
    gdmt = SQRT(gdt(1)**2+gdt(2)**2+gdt(3)**2)
    gdmt = MAX(GDMIN,gdmt)
 
    ! Find local correlation energy and potential
    CALL PW92C(2, d, ecunif, vcunif)
 
    ! Find total correlation energy
    rs = (3/(4*pi*dt))**THD
    kf = (3*pi**2*dt)**THD
    ks = SQRT(4*kf/pi)
    zeta = (d(1)-d(2))/dt
    zeta = MAX(-1.D0+DENMIN,zeta)
    zeta = MIN(1.D0-DENMIN,zeta)
    phi = HALF*((1+zeta)**TWOTHD+(1-zeta)**TWOTHD)
    t = gdmt/(2*phi*ks*dt)
    f1 = ecunif/gamma/phi**3
    f2 = EXP(-f1)
    a = beta/gamma/(f2-1)
    f3 = t**2 + a*t**4
    f4 = beta/gamma*f3/(1+a*f3)
    h = gamma*phi**3*LOG(1+f4)
    fc = ecunif + h
 
    ! Find correlation energy derivatives
    drsdd = -(THD*rs/dt)
    dkfdd = THD*kf/dt
    dksdd = HALF*ks*dkfdd/kf
    dzdd(1) = 1/dt - zeta/dt
    dzdd(2) = -(1/dt) - zeta/dt
    dpdz = HALF*TWOTHD*(1/(1+zeta)**THD-1/(1-zeta)**THD)
    DO is = 1, 2
       decudd = (vcunif(is)-ecunif)/dt
       dpdd = dpdz*dzdd(is)
       dtdd = (-t)*(dpdd/phi+dksdd/ks+1/dt)
       df1dd = f1*(decudd/ecunif-3*dpdd/phi)
       df2dd = (-f2)*df1dd
       dadd = (-a)*df2dd/(f2-1)
       df3dd = (2*t+4*a*t**3)*dtdd + dadd*t**4
       df4dd = f4*(df3dd/f3-(dadd*f3+a*df3dd)/(1+a*f3))
       dhdd = 3*h*dpdd/phi
       dhdd = dhdd + gamma*phi**3*df4dd/(1+f4)
       dfcdd(is) = vcunif(is) + h + dt*dhdd
 
       DO ix = 1, 3
          dtdgd = (t/gdmt)*gdt(ix)/gdmt
          df3dgd = dtdgd*(2*t+4*a*t**3)
          df4dgd = f4*df3dgd*(1/f3-a/(1+a*f3))
          dhdgd = gamma*phi**3*df4dgd/(1+f4)
          dfcdgd(ix,is) = dt*dhdgd
       ENDDO
    ENDDO
 
    ! Find exchange energy and potential
    CALL XWPBE(Nspin, Dens, Gdens, fx_cp2k, dfxdd_cp2k, dfxdgd_cp2k)
    fx_cp2k = fx_cp2k/dt
 
    ! Find exchange energy and potential
    fx = 0
    DO is = 1, 2
       ds(is) = MAX(DENMIN,2*d(is))
       gdms = MAX(GDMIN,2*gdm(is))
       kfs = (3*pi**2*ds(is))**THD
       s = gdms/(2*kfs*ds(is))
       f1 = 1 + mu*s**2/kappa
       f = 1 + kappa - kappa/f1
       !
       ! Note nspin=1 in call to exchng...
       !
       CALL EXCHNG(Irel, 1, ds(is), exunif, vxunif(is))
       fx = fx + ds(is)*exunif*f
 
       dkfdd = THD*kfs/ds(is)
       dsdd = s*(-(dkfdd/kfs)-1/ds(is))
       df1dd = 2*(f1-1)*dsdd/s
       dfdd = kappa*df1dd/f1**2
       dfxdd(is) = vxunif(is)*f + ds(is)*exunif*dfdd
 
       DO ix = 1, 3
          gds = 2*gd(ix,is)
          dsdgd = (s/gdms)*gds/gdms
          df1dgd = 2*mu*s*dsdgd/kappa
          dfdgd = kappa*df1dgd/f1**2
          dfxdgd(ix,is) = ds(is)*exunif*dfdgd
       ENDDO
    ENDDO
    fx = HALF*fx/dt
 
    ! Set output arguments
    Ex = fx_cp2k
    Ec = fc
    DO is = 1, Nspin
       Dexdd(is) = dfxdd_cp2k(is)
       Decdd(is) = dfcdd(is)
       DO ix = 1, 3
          Dexdgd(ix,is) = dfxdgd_cp2k(ix,is)
          Decdgd(ix,is) = dfcdgd(ix,is)
       ENDDO
    ENDDO
 
  END SUBROUTINE HSEXC

  ! ********************************************************************
  ! Implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
  ! Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
  ! Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
  ! Add an fator=0.75 for Exchange term, Correlaton term is unchanged.
  ! ******** INPUT *****************************************************
  ! INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
  ! INTEGER nspin          : Number of spin polarizations (1 or 2)
  ! REAL*8  Dens(nspin)    : Total electron density (if nspin=1) or
  !                           spin electron density (if nspin=2)
  ! REAL*8  GDens(3,nspin) : Total or spin density gradient
  ! ******** OUTPUT ****************************************************
  ! REAL*8  EX             : Exchange energy density
  ! REAL*8  EC             : Correlation energy density
  ! REAL*8  DEXDD(nspin)   : Partial derivative
  !                           d(DensTot*Ex)/dDens(ispin),
  !                           where DensTot = Sum_ispin( Dens(ispin) )
  !                          For a constant density, this is the
  !                          exchange potential
  ! REAL*8  DECDD(nspin)   : Partial derivative
  !                           d(DensTot*Ec)/dDens(ispin),
  !                           where DensTot = Sum_ispin( Dens(ispin) )
  !                          For a constant density, this is the
  !                          correlation potential
  ! REAL*8  DEXDGD(3,nspin): Partial derivative
  !                           d(DensTot*Ex)/d(GradDens(i,ispin))
  ! REAL*8  DECDGD(3,nspin): Partial derivative
  !                           d(DensTot*Ec)/d(GradDens(i,ispin))
  ! ********* UNITS ****************************************************
  ! Lengths in Bohr
  ! Densities in electrons per Bohr**3
  ! Energies in Hartrees
  ! Gradient vectors in cartesian coordinates
  ! ********* ROUTINES CALLED ******************************************
  ! EXCHNG, PW92C
  ! ********************************************************************
  SUBROUTINE PBE0XC(Irel, Nspin, Dens, Gdens, Ex, Ec, Dexdd, Decdd, &
               Dexdgd, Decdgd)
 
    USE gridxc_precision, ONLY: dp
    USE gridxc_lda, ONLY: EXCHNG, PW92C

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Irel, Nspin
    REAL(dp), INTENT(IN) :: Dens(Nspin), Gdens(3,Nspin)
    REAL(dp), INTENT(INOUT) :: Decdd(Nspin), Decdgd(3,Nspin), &
      Dexdd(Nspin), Dexdgd(3,Nspin), Ec, Ex
 
    ! Internal variables
    INTEGER :: is, ix
 
    REAL(dp) :: a, beta, d(2), dadd, decudd, df1dd,    &
      df2dd, df3dd, df4dd, df1dgd, df3dgd, df4dgd,       &
      dfcdd(2), dfcdgd(3,2), dfdd, dfdgd, dfxdd(2),       &
      dfxdgd(3,2), dhdd, dhdgd, dkfdd, dksdd, dpdd,      &
      dpdz, drsdd, ds(2), dsdd, dsdgd, dt, dtdd, dtdgd,&
      dzdd(2), ecunif, exunif, f, f1, f2, f3, &
      f4, fc, fx, gamma, gd(3,2), gdm(2),       &
      gdms, gdmt, gds, gdt(3), h, kappa,  &
      kf, kfs, ks, mu, phi, pi, rs, s, t,       &
      vcunif(2), vxunif(2), zeta
 
    ! Lower bounds of density and its gradient to avoid divisions by zero
    REAL(dp), PARAMETER :: DENMIN=1.D-12
    REAL(dp), PARAMETER :: GDMIN=1.D-12
 
    ! Fix some numerical parameters
    REAL(dp), PARAMETER :: FOUTHD=4.D0/3.D0
    REAL(dp), PARAMETER :: HALF=0.5D0
    REAL(dp), PARAMETER :: THD=1.D0/3.D0
    REAL(dp), PARAMETER :: THRHLF=1.5D0
    REAL(dp), PARAMETER :: TWO=2.D0
    REAL(dp), PARAMETER :: TWOTHD=2.D0/3.D0
 
    ! Exchange factor 3/4---XC_PBE0=3/4*X_PBE+C_PBE+1/4X_HFX
    REAL(dp), PARAMETER :: FACTOR=0.75D0
 
    ! Fix some more numerical constants
    pi = 4*ATAN(1.D0)
    beta = 0.066725D0
    gamma = (1-LOG(TWO))/pi**2
    mu = beta*pi**2/3
    kappa = 0.804D0
 
    ! Translate density and its gradient to new variables
    IF ( Nspin==1 ) THEN
       d(1) = HALF*Dens(1)
       d(2) = d(1)
       dt = MAX(DENMIN,Dens(1))
       DO ix = 1, 3
          gd(ix,1) = HALF*Gdens(ix,1)
          gd(ix,2) = gd(ix,1)
          gdt(ix) = Gdens(ix,1)
       ENDDO
    ELSE
       d(1) = Dens(1)
       d(2) = Dens(2)
       dt = MAX(DENMIN,Dens(1)+Dens(2))
       DO ix = 1, 3
          gd(ix,1) = Gdens(ix,1)
          gd(ix,2) = Gdens(ix,2)
          gdt(ix) = Gdens(ix,1) + Gdens(ix,2)
       ENDDO
    ENDIF
    gdm(1) = SQRT(gd(1,1)**2+gd(2,1)**2+gd(3,1)**2)
    gdm(2) = SQRT(gd(1,2)**2+gd(2,2)**2+gd(3,2)**2)
    gdmt = SQRT(gdt(1)**2+gdt(2)**2+gdt(3)**2)
    gdmt = MAX(GDMIN,gdmt)
 
    ! Find local correlation energy and potential
    CALL PW92C(2,d,ecunif,vcunif)
 
    ! Find total correlation energy
    rs = (3/(4*pi*dt))**THD
    kf = (3*pi**2*dt)**THD
    ks = SQRT(4*kf/pi)
    zeta = (d(1)-d(2))/dt
    zeta = MAX(-1.D0+DENMIN,zeta)
    zeta = MIN(1.D0-DENMIN,zeta)
    phi = HALF*((1+zeta)**TWOTHD+(1-zeta)**TWOTHD)
    t = gdmt/(2*phi*ks*dt)
    f1 = ecunif/gamma/phi**3
    f2 = EXP(-f1)
    a = beta/gamma/(f2-1)
    f3 = t**2 + a*t**4
    f4 = beta/gamma*f3/(1+a*f3)
    h = gamma*phi**3*LOG(1+f4)
    fc = ecunif + h
 
    ! Find correlation energy derivatives
    drsdd = -(THD*rs/dt)
    dkfdd = THD*kf/dt
    dksdd = HALF*ks*dkfdd/kf
    dzdd(1) = 1/dt - zeta/dt
    dzdd(2) = -(1/dt) - zeta/dt
    dpdz = HALF*TWOTHD*(1/(1+zeta)**THD-1/(1-zeta)**THD)
    DO is = 1, 2
       decudd = (vcunif(is)-ecunif)/dt
       dpdd = dpdz*dzdd(is)
       dtdd = (-t)*(dpdd/phi+dksdd/ks+1/dt)
       df1dd = f1*(decudd/ecunif-3*dpdd/phi)
       df2dd = (-f2)*df1dd
       dadd = (-a)*df2dd/(f2-1)
       df3dd = (2*t+4*a*t**3)*dtdd + dadd*t**4
       df4dd = f4*(df3dd/f3-(dadd*f3+a*df3dd)/(1+a*f3))
       dhdd = 3*h*dpdd/phi
       dhdd = dhdd + gamma*phi**3*df4dd/(1+f4)
       dfcdd(is) = vcunif(is) + h + dt*dhdd
 
       DO ix = 1, 3
          dtdgd = (t/gdmt)*gdt(ix)/gdmt
          df3dgd = dtdgd*(2*t+4*a*t**3)
          df4dgd = f4*df3dgd*(1/f3-a/(1+a*f3))
          dhdgd = gamma*phi**3*df4dgd/(1+f4)
          dfcdgd(ix,is) = dt*dhdgd
       ENDDO
    ENDDO
 
    ! Find exchange energy and potential
    fx = 0
    DO is = 1, 2
       ds(is) = MAX(DENMIN,2*d(is))
       gdms = MAX(GDMIN,2*gdm(is))
       kfs = (3*pi**2*ds(is))**THD
       s = gdms/(2*kfs*ds(is))
       f1 = 1 + mu*s**2/kappa
       f = 1 + kappa - kappa/f1
!
!       Note nspin=1 in call to exchng...
!
       CALL EXCHNG(Irel,1,ds(is),exunif,vxunif(is))
       fx = fx + ds(is)*exunif*f
 
       dkfdd = THD*kfs/ds(is)
       dsdd = s*(-(dkfdd/kfs)-1/ds(is))
       df1dd = 2*(f1-1)*dsdd/s
       dfdd = kappa*df1dd/f1**2
       dfxdd(is) = vxunif(is)*f + ds(is)*exunif*dfdd
 
       DO ix = 1, 3
          gds = 2*gd(ix,is)
          dsdgd = (s/gdms)*gds/gdms
          df1dgd = 2*mu*s*dsdgd/kappa
          dfdgd = kappa*df1dgd/f1**2
          dfxdgd(ix,is) = ds(is)*exunif*dfdgd
       ENDDO
    ENDDO
    fx = HALF*fx/dt
 
    ! Set output arguments
    Ex = fx*FACTOR
    Ec = fc
    DO is = 1, Nspin
       Dexdd(is) = dfxdd(is)*FACTOR
       Decdd(is) = dfcdd(is)
       DO ix = 1, 3
          Dexdgd(ix,is) = dfxdgd(ix,is)*FACTOR
          Decdgd(ix,is) = dfcdgd(ix,is)
       ENDDO
    ENDDO
 
  END SUBROUTINE PBE0XC

END MODULE gridxc_hybrids
