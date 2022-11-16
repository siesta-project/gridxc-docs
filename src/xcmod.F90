#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_LIBXC
#include "xc_version.h"
#endif

!< Stores the information about the XC functional to use,
! and provides routines to set and to get it.
!
! Note: data is currently global. In future, it will be put into a derived type and
! initialized and passed around in a single handle.
!
module gridxc_xcmod

  use gridxc_precision, only: dp              ! Double precision real kind
  use gridxc_sys,       only: die             ! Termination routine
  use gridxc_vdwxc,   only: vdw_set_author  ! Sets vdW functional flavour

  implicit none

public:: &
  setXC_family_authors, &! Sets single XC functional in family/author style
  setXC, &! Sets XC functional(s) to be used (comprehensive)
  getXC   ! Returns the XC functional(s) being used

#ifdef HAVE_LIBXC
public :: setXC_libxc_ids! Sets XC functionals using libxc ids
#endif

private ! Nothing is declared public beyond this point

! These data should be put into a derived type and
! initialized and passed around in a single handle, instead
! of being global

  integer, parameter :: maxFunc = 20
  integer,           save :: nXCfunc=0
  character(len=50), save :: XCauth(MaxFunc), XCfunc(MaxFunc)
  real(dp),          save :: XCweightX(MaxFunc), XCweightC(MaxFunc)

contains

  
!< Sets the xc functional(s) to be used by atomxc and/or cellxc. 
! The allowed functional/author values can be seen [here](|url|/page/02_user_guide/02_api/functionals.html).
  
!< Usage example:
!  
!     call setXC( 1, (/'GGA'/), (/'PBE'/), (/1._dp/), (/1._dp/) )
!  
!    Stops with an error message if n is larger than the internal module parameter maxFunc


  subroutine setXC( n, func, auth, wx, wc )
    implicit none

    !> Number of functionals to use
    integer,         intent(in):: n
    !> Functional family labels (LDA, GGA, etc)
    character(len=*),intent(in):: func(n) 
    !> Functional author labels
    character(len=*),intent(in):: auth(n) 
    !> Functional weights for exchange  (sum(wx) must be 1)
    real(dp),        intent(in):: wx(n)   
    !> Functional weights for correlation (sum(wc) must be 1)
    real(dp),        intent(in):: wc(n)

    integer:: i, j

    if (n>maxFunc) call die('setXC: ERROR: parameter maxFunc too small')
    nXCfunc = n
    XCfunc(1:n) = func(1:n)
    XCauth(1:n) = auth(1:n)
    XCweightX(1:n) = wx(1:n)
    XCweightC(1:n) = wc(1:n)
    do i = 1,n
      if (XCfunc(i)=='VDW' .or. XCfunc(i)=='vdw' .or. XCfunc(i)=='vdW') then
        XCfunc(i) = 'VDW'
        do j = 1,i-1
          if (XCfunc(j)=='VDW' .and. XCauth(j)/=XCauth(i)) &
            call die('setXC ERROR: mixing different VDW authors not allowed')
        end do ! j
        call vdw_set_author( XCauth(i) )
      end if ! (XCfunc(i)=='VDW')
     
      if (XCauth(i)(1:6) == "LIBXC-") then
         call process_libxc_spec(XCfunc(i),XCauth(i))
      endif
      
    end do ! i
   
  end subroutine setXC

  !< Sets a single XC functional in family/author style.
  !  It can only be used for built-in functionals and
  !  'XC' libXC functionals
  subroutine setXC_family_authors( family, auth )

    implicit none
    character(len=*),intent(in):: family
    character(len=*),intent(in):: auth

    call setXC(1, [family], [auth], [1.0_dp], [1.0_dp])
  end subroutine setXC_family_authors

#ifdef HAVE_LIBXC
!< Sets the XC info using libxc numerical codes.
!  It can only be used if LibXC support is compiled in
  
  subroutine setXC_libxc_ids( nfuncs, libxc_ids)
#if XC_MAJOR_VERSION >= 4
    use xc_f03_lib_m
#else 
    use xc_f90_types_m
    use xc_f90_lib_m
#endif

    implicit none
    !> number of functionals
    integer, intent(in) :: nfuncs
    !> numerical libxc codes
    integer, intent(in) :: libxc_ids(nfuncs)

    ! automatic arrays
    character(len=11)   :: family(nfuncs), auth(nfuncs)
    real(dp)            :: weight_x(nfuncs), weight_c(nfuncs)
    
#if XC_MAJOR_VERSION >= 4
    type(xc_f03_func_t) :: xc_func
    type(xc_f03_func_info_t) :: xc_info
#else
    type(xc_f90_pointer_t) :: xc_func, xc_info
#endif
    integer :: xc_ispin
    integer :: kind, ifamily
    integer :: i

    do i = 1, nfuncs
       !
       ! Determine the kind of functional to assign weights correctly
       !
       xc_ispin = XC_UNPOLARIZED 
       ! 'unpolarized' is the least stringent option: for non-polarized
       ! functionals, the other option might result in an error
#if XC_MAJOR_VERSION >= 4
       call xc_f03_func_init(xc_func, libxc_ids(i), xc_ispin)
       xc_info = xc_f03_func_get_info(xc_func)
       kind = xc_f03_func_info_get_kind(xc_info)
       ifamily = xc_f03_func_info_get_family(xc_info)
#else
       call xc_f90_func_init(xc_func, xc_info, libxc_ids(i), xc_ispin)
       kind = xc_f90_info_kind(xc_info)
       ifamily = xc_f90_family_from_id(libxc_ids(i))
#endif
       
       select case (kind)
       case (XC_CORRELATION)
          weight_x(i) = 0.0_dp
          weight_c(i) = 1.0_dp
       case (XC_EXCHANGE)
          weight_x(i) = 1.0_dp
          weight_c(i) = 0.0_dp
       case (XC_EXCHANGE_CORRELATION)
          weight_x(i) = 1.0_dp
          weight_c(i) = 1.0_dp
       case default
          call die("Functional kind not supported")
       end select
#if XC_MAJOR_VERSION >= 4
       call xc_f03_func_end(xc_func)
#else
       call xc_f90_func_end(xc_func)
#endif

       select case ( ifamily )
       case (XC_FAMILY_LDA)
          family(i) = "LDA"
       case (XC_FAMILY_GGA)
          family(i) = "GGA"
       case (XC_FAMILY_HYB_GGA)
          family(i) = "GGA"
       ! There is probably a (negative) case for bad id ...
       case (XC_FAMILY_UNKNOWN)
          call die("Bad libxc functional code")
       case default
          family(i) = "other"
       end select
       write(auth(i),"(a,i5.5)") "LIBXC-", libxc_ids(i)
    end do

    if (sum(weight_x(1:nfuncs)) /= 1.0_dp) then
       call die("Wrong exchange weights")
    endif
    if (sum(weight_c(1:nfuncs)) /= 1.0_dp) then
       call die("Wrong correlation weights")
    endif
    
    call setXC(nfuncs, family, auth, weight_x, weight_c)
    
  end subroutine setXC_libxc_ids
# endif /* HAVE_LIBXC */

!< Returns the xc functional(s) information previously set and stored
!  in the module variables.
  
  subroutine getXC( n, func, auth, wx, wc )
    implicit none
    !> Number of functionals
    integer,         optional,intent(out):: n       
    !> Functional name labels
    character(len=*),optional,intent(out):: func(:) 
    !> Functional author labels
    character(len=*),optional,intent(out):: auth(:) 
    !> Functl weights for exchange
    real(dp),        optional,intent(out):: wx(:)   
    !> Functl weights for correlation
    real(dp),        optional,intent(out):: wc(:)   

    integer:: nf
    nf = nXCfunc
    if (present(n)) n = nf
    if (present(func)) then
       if (size(func)>=nf) func(1:nf) = XCfunc(1:nf)
    end if
    if (present(auth)) then
       if (size(auth)>=nf) auth(1:nf) = XCauth(1:nf)
    end if
    if (present(wx)) then
       if (size(wx)  >=nf) wx(1:nf)   = XCweightX(1:nf)
    end if
    if (present(wc)) then
       if (size(wc)  >=nf) wc(1:nf)   = XCweightC(1:nf)
    end if
  end subroutine getXC

#ifndef HAVE_LIBXC               
  subroutine process_libxc_spec(func,auth)
    character(len=*), intent(in)    :: func
    character(len=*), intent(inout) ::  auth

    call die("Libxc not compiled in. Cannot handle " //  &
              trim(func) // " " // trim(auth))
  end subroutine process_libxc_spec
  
#else
  
  subroutine process_libxc_spec(func,auth)

#if XC_MAJOR_VERSION >= 4
    use xc_f03_lib_m
#else
    use xc_f90_types_m
    use xc_f90_lib_m
#endif

    character(len=*), intent(in)    :: func
    character(len=*), intent(inout) ::  auth

    integer :: iostat, xc_id, idx, xc_id_from_symbol
    integer :: family
    character(len=50) :: symbolic_name

    ! Fields are of the form LIBXC-XXXX-SYMBOL
    ! where -SYMBOL is optional if XXXX is a meaningful code
    idx = index(auth(7:),"-")
    if (idx /= 0) then
      ! We have code and symbol fields
      read(auth(7:7+idx-2),iostat=iostat,fmt=*) xc_id
      symbolic_name = auth(7+idx:)
#if XC_MAJOR_VERSION >= 4
      xc_id_from_symbol = xc_f03_functional_get_number(symbolic_name)
#else
      xc_id_from_symbol = xc_f90_functional_get_number(symbolic_name)
#endif
      if (xc_id == 0) then
        ! A zero in the code field signals that we want
        ! to fall back on the symbolic name field 
        if (xc_id_from_symbol < 0) then
          call die("Cannot get xc_id from " // &
              trim(symbolic_name))
        else
          xc_id = xc_id_from_symbol
        endif
      else
        ! Check consistency
        if (xc_id /= xc_id_from_symbol) then
          call die("Conflicting code field for " // &
              trim(symbolic_name))
        endif
      endif
      ! Normalize the internal representation 
      write(auth,"(a,i5.5,'-',a)") &
          "LIBXC-",xc_id, trim(symbolic_name)
    else
      ! Just a code field
      read(auth(7:),iostat=iostat,fmt=*) xc_id
      if (iostat /= 0) call die("Bad libxc code in " &
          // trim(auth))
      ! Normalize the internal representation 
      write(auth,"(a,i5.5)") "LIBXC-",xc_id
    endif

#if XC_MAJOR_VERSION >= 4
    family = xc_f03_family_from_id (xc_id)
#else
    family = xc_f90_family_from_id (xc_id)
#endif

    select case ( family )
    case (XC_FAMILY_LDA)
      if (func /= "LDA") call die("Family mismatch in " // &
          trim(func) // " " // trim(auth))
    case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      if (func /= "GGA") call die("Family mismatch in " // &
          trim(func) // " " // trim(auth))
    case default
      call die("Unsupported Libxc family or functional")
    end select

  end subroutine process_libxc_spec

#endif
  
end module gridxc_xcmod
