module qg_arrays                !-*-f90-*-

  !************************************************************************
  ! This module contains all of the dynamically allocated field variables,
  ! and the routine to allocate memory for them.  Also sets the preset
  ! K^2 kx and ky grids.
  !
  ! Routines: Setup_fields
  !
  ! Dependencies: qg_params, io_tools
  !************************************************************************

  ! *******************
  ! 26/1/2022
  ! added ktdrag_=ksqd**(0.5*(1.-n_hypo)) for hypofriction
  ! *******************

  implicit none

  complex,dimension(:,:),allocatable   :: q,psi,rhs,thetag,div
  complex,dimension(:,:,:),allocatable   :: u_a,u_phi,uvd_total
  complex,dimension(:,:),  allocatable   :: force_o
  real,   dimension(:,:),  allocatable   :: filter
  real,   dimension(:,:),  allocatable   :: ksqd_,kx_,ky_,ktdrag_
  integer,dimension(:),    allocatable   :: kxv, kyv, lin2kx, lin2ky,psiq
! *************************************************************
! 28/3/2022
  real,   dimension(:),  allocatable   :: xtra,ytra
  real,   dimension(:),  allocatable   :: xtra_1,ytra_1
  real,   dimension(:),  allocatable   :: xtra_0,ytra_0
  real,   dimension(:,:),  allocatable   :: u,v
  real,   dimension(:,:),  allocatable   :: u_1,v_1
  real,   dimension(:,:),  allocatable   :: u_0,v_0
! *************************************************************


  save


contains

  subroutine Setup_fields

    use qg_params, only: kmax,use_forcing,&
                         nmask,nx,ny,n_hypo,&
                         ntra,lagrangian
    integer :: kx, ky

    allocate(q(-kmax:kmax,0:kmax))
    allocate(rhs(-kmax:kmax,0:kmax),uvd_total(nx,ny,3))
    allocate(u_a(-kmax:kmax,0:kmax,2),u_phi(-kmax:kmax,0:kmax,2))
    allocate(filter(-kmax:kmax,0:kmax))
    allocate(ksqd_(-kmax:kmax,0:kmax),kx_(-kmax:kmax,0:kmax))
    allocate(ky_(-kmax:kmax,0:kmax))
    allocate(ktdrag_(-kmax:kmax,0:kmax))
    allocate(kxv(2*kmax+1),kyv(kmax+1))
    q=0.; rhs=0. 

! *************************************************************
! 26/3/2019-3/4/2019
! arrays for Lagrangian trajectories
    if(lagrangian) then
        allocate(xtra(ntra),ytra(ntra))
        allocate(xtra_1(ntra),ytra_1(ntra))
        allocate(xtra_0(ntra),ytra_0(ntra))
        allocate(u(nx,ny),v(nx,ny))
        allocate(u_1(nx,ny),v_1(nx,ny))
        allocate(u_0(nx,ny),v_0(nx,ny))
        u=0.;v=0.
        u_1=0.;v_1=0.
        u_0=0.;v_0=0.
    endif
! *************************************************************

    allocate(psi(-kmax:kmax,0:kmax))

    allocate(thetag(nx,ny))
    allocate(div(-kmax:kmax,0:kmax))
    psi = 0.

    if (use_forcing) then
       allocate(force_o(-kmax:kmax,0:kmax))
       force_o = 0.
    endif

    ! Store values of kx, ky and k^2 in 2d arrays.

    kxv = (/ (kx,kx=-kmax,kmax) /)
    kyv = (/ (ky,ky=0,kmax) /)
    kx_ = float(spread(kxv,2,kmax+1))
    ky_ = float(spread(kyv,1,2*kmax+1))
    ksqd_ = kx_**2 + ky_**2
    ksqd_(0,0) = .1    ! 0,0 never used - this way can divide by K^2 w/o worry

    ! **********************
    ! 26/1/2022
    ktdrag_=ksqd_**(0.5*(1.-n_hypo))
    ! **********************

    ! Get # of true elements in de-aliasing mask and set up index maps so
    ! that we can choose only to invert dynamic part of psi

    filter = 1.                ! Set up de-aliasing part of filter first
    where (ksqd_>= (8./9.)*(kmax+1)**2) filter = 0.
    filter(-kmax:0,0) = 0.

    nmask = sum(int(filter))   ! Get non-zero total of filter

    allocate(lin2kx(nmask),lin2ky(nmask))
    lin2kx = pack(kx_,MASK=(filter>0))
    lin2ky = pack(ky_,MASK=(filter>0))

  end subroutine Setup_fields

end module qg_arrays
