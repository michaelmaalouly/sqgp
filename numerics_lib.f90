module numerics_lib

  !*********************************************************************
  ! This module contains an assortment of low level numerical routines
  ! The idea is to make a close correspondance with what one can do naturally
  ! in Matlab.
  !
  ! Routines: Ring_integral, Spline, Splint, Sortvec, Findzeros,
  !           sub2ind, Ran, Tridiag, Tridiag_cyc, Tri2mat, Diag
  !
  ! Dependencies: io_tools (but may use op_rules in future)
  ! 
  ! See headers on routines for descriptions
  !
  !*********************************************************************

  implicit none
  private

  public :: Ring_integral, Spline, Splint, Sortvec, Findzeros, Sub2ind, &
       Ran, ranp, Diag, Tridiag, Tridiag_cyc, Tri2mat

  interface Ring_integral
     module procedure Ring_integral2, Ring_integral3
  end interface

  interface Diag
     module procedure Mat2diag, Diag2mat
  end interface

contains

  !*********************************************************************

  function Ring_integral2(f,xv,yv,nf)  result(int)

    !*********************************************************************
    ! Calculate integral of 2d field f over angle, where 
    ! f = f(x:[xv],y:[yv]), xv and yv are rank 1 integer arrays
    !*********************************************************************

    real,dimension(:,:),intent(in)         :: f
    integer,dimension(:),intent(in)        :: xv, yv
    integer,intent(in)                     :: nf
    real,dimension(nf)                     :: int
    ! local
    integer                                :: nx, ny, rad
    integer,dimension(size(f,1),size(f,2)) :: radius_arr

    nx = size(xv); ny = size(yv)
    radius_arr = floor(sqrt(float(spread(xv,2,ny)**2+spread(yv,1,nx)**2))+.5)
    int = (/ (sum(f,MASK=(radius_arr-rad==0)), rad=1,nf) /)
      
  end function Ring_integral2

  !*********************************************************************

  function Ring_integral3(f,xv,yv,nf)  result(int)

    !*********************************************************************
    ! Calculate integral of 2d field f over angle, where 
    ! f = f(x:[xv],y:[yv]), xv and yv are rank 1 integer arrays
    !*********************************************************************

    real,dimension(:,:,:),intent(in)       :: f
    integer,dimension(:),intent(in)        :: xv, yv
    integer,intent(in)                     :: nf
    real,dimension(1:nf,size(f,3))         :: int
    ! local
    integer                                :: n, nx, ny, rad
    integer,dimension(size(f,1),size(f,2)) :: radius_arr

    nx = size(xv); ny = size(yv)

    radius_arr = floor(sqrt(float(spread(xv,2,ny)**2+spread(yv,1,nx)**2))+.5)

    do n = 1,size(f,3)
       int(:,n) = (/ (sum(f(:,:,n),MASK=(radius_arr-rad==0)), rad=1,nf) /)
    enddo

  end function Ring_integral3

  !*********************************************************************

  function Spline(x,y) result(y2)

    !*********************************************************************
    ! Implementation of Numerical Recipies Spline calculation.
    ! y2 = projection of x onto y
    !*********************************************************************
    
    real,dimension(:),intent(in) :: x,y
    real,dimension(size(x,1))    :: y2,u
    integer                      :: nx, n
    real                         :: sig, p

    nx = size(x)

    y2 = 0.      ! Use natural spline BC
    u  = 0. 

    do n = 2,nx-1
       sig = (x(n)-x(n-1))/(x(n+1)-x(n-1))
       p   = sig*y2(n-1) + 2
       y2(n) = (sig-1.)/p
       u(n) = (6. * &
           ( (y(n+1)-y(n))/(x(n+1)-x(n)) - (y(n)-y(n-1))/(x(n)-x(n-1)) ) &
           /(x(n+1)-x(n-1)) - sig*u(n-1))/p
    enddo

    do n = nx-1,1,-1
       y2(n) = y2(n)*y2(n+1)+u(n)
    enddo

  end function Spline

  !*********************************************************************

  function Splint(xv,yv,y2v,x) result(y)

    !*********************************************************************
    ! Implementation of Numerical Recipies Splint calculation for Spline
    ! Computes array of second derivatives, d2f/dy2, of yvec = f(xvec).
    !*********************************************************************
    
    use io_tools, only: Message

    real,dimension(:),intent(in) :: xv,yv,y2v
    real,intent(in)              :: x
    real                         :: y
    integer                      :: nx, k, klo, khi
    real                         :: h, a, b

    nx = size(xv)

    klo = 1
    khi = nx
    do while (khi-klo > 1) 
       k = int((khi+klo)/2)
       if (xv(k) > x) then
          khi = k
       else
          klo = k
       endif
    enddo

    h = xv(khi) - xv(klo)
    if (h==0) call Message('Splint: bad xv',fatal='y')
    a = (xv(khi)-x)/h
    b = (x-xv(klo))/h
    y = a*yv(klo) + b*yv(khi) + ((a**3-a)*y2v(klo)+(b**3-b)*y2v(khi))*(h**2)/6.

  end function Splint

  !*************************************************************************

  function Findzeros(f,x)

    !*********************************************************************
    ! Returns a vector containing the best approximation 
    ! of the zeros of the 1-d vector f.  Optional argument x is the  f
    ! same size asis assumed to contain the values of the independent 
    ! variable of whichf is a function.  Uses simple linear interpolation, 
    ! so this is intended for high resolution functions.
    !*********************************************************************

    use io_tools

    real,dimension(:),intent(IN)        :: f,x

    ! Local
    real,dimension(size(f))             :: Findzeros
    integer,dimension(:),allocatable    :: s
    real,dimension(:),allocatable       :: v
    integer                             :: numzs, n, nv, j
    real                                :: slope

    if (size(f)/=size(x)) call Message('Findzeros: Bad vector sizes',fatal='y')
    nv = size(f)
    allocate(s(nv),v(nv))

    s=0; v=0; numzs=0

    do n=2,nv
       if (f(n)/=sign(f(n),f(n-1))) then ! Sign change occured
          numzs=numzs+1
          s(numzs)=n                     ! Collect indeces just after sign chg
       endif
    enddo

    do j=1,numzs
       n = s(j)                             ! Loop through indeces of sign chgs
       slope = (f(n)-f(n-1))/(x(n)-x(n-1))  ! Get local slope
       v(j) = -f(n-1)/slope + x(n-1)        ! Linear interp of f wrt x to 0
    enddo

    Findzeros = v

    deallocate(s,v)

  end function Findzeros

  !*********************************************************************

  subroutine Sortvec(vin,ind_map)

    !*********************************************************************
    ! Sort rank 1 vector vin in ascending order and store in vout.
    ! Store map of initial indeces to final indeces in ind_map,
    ! such that vout(ind_map(i)) = vin(i)
    !*********************************************************************

    real,dimension(:),intent(inout)   :: vin
    integer,dimension(:),intent(out)  :: ind_map

    ! Local
    real,dimension(size(vin))         :: vout
    real                              :: temp
    integer                           :: i,j,nv

    nv = size(vin)
    vout = vin

    do j = 2,nv                      ! Vintage bubble sort
       do i = j,2,-1
          if (vout(i) < vout(i-1)) then
             temp = vout(i)
             vout(i) = vout(i-1)
             vout(i-1) = temp
          endif
       enddo
    enddo

    do j = 1,nv        
       do i = 1,nv
          if (vout(i) == vin(j)) ind_map(i) = j
       enddo
    enddo

    vin = vout

  end subroutine Sortvec

  !*********************************************************************

  function Sub2ind(subarr,nmax)  result(ind)

    !*********************************************************************
    ! Convert an N-dimensional array subscript, passed in array 'subarr', 
    ! each of whose maximum extent is nmax, to a single index.  
    ! size(subarr) = N
    !*********************************************************************

    integer,dimension(:),intent(in) :: subarr
    integer,intent(in)              :: nmax
    integer                         :: ind
    integer                         :: j

    ind = 1 + sum((/ (nmax**j, j=0,size(subarr)-1) /)*(subarr-1))

  end function Sub2ind

  !*********************************************************************

  function Ran(idum,nx,ny) result(ranf)

    !*********************************************************************
    ! Random number generator -- use this one so that one can repeat runs 
    ! exactly, independant of machine
    !*********************************************************************

    integer,intent(inout)           :: idum
    integer,intent(in)              :: nx, ny
    real,dimension(nx,ny)           :: ranf

    ! Local
    integer                         :: j,x,y
    integer,parameter               :: m=714025,ia=1366,ic=150889,irsize=97
    real,parameter                  :: rm=1.4005112e-6
    integer,dimension(irsize),save  :: ir
    integer,save                    :: iff = 0, iy
    
    if ((idum<0).or.(iff==0)) then
       iff=1
       idum=mod(ic-idum,m)
       do j=1,irsize
          idum=mod(ia*idum+ic,m)
          ir(j)=idum
       enddo
       idum=mod(ia*idum+ic,m)
       iy=idum
    endif

    do y = 1,ny
       do x = 1,nx    
          j=1+(irsize*iy)/m
          if ((j>irsize).or.(j<1)) stop 'j out of bounds in ran2'
          iy=ir(j)
          ranf(x,y)=iy*rm
          idum=mod(ia*idum+ic,m)
          ir(j)=idum
       enddo
    enddo

  end function Ran

! **********************************************************************
! 20/4/2022
! Uniform random number generator
! **********************************************************************
! Returns random number uniform in [0,1]
! period > 2e18
! call with idump a negative integer to initialize
! do not alter idump in successive calls

function ranp(idump)
!implicit none

integer,intent(inout) :: idump
integer,parameter :: im1=2147483563,im2=2147483399,imm1=im1-1,ia1= 40014,ia2=40692
integer,parameter :: iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32,ndiv=1+imm1/ntab
real*8, parameter :: am=1.d0/im1,eps=1.d-7,rnmx=1.-eps
integer :: idum2,j,k,iv(ntab),iy
real*8 ::ranp
save iv,iy,idum2
data idum2/123456789/,iv/ntab*0/,iy/0/

if (idump<=0) then
 idump=max(-idump,1)
 idum2=idump
 do j=ntab+8,1,-1
  k=idump/iq1
  idump=ia1*(idump-k*iq1)-k*ir1
  if (idump<0) idump=idump+im1
  if (j<=ntab) iv(j)=idump
 end do
 iy=iv(1)
end if

! start here when not initializing
k=idump/iq1
idump=ia1*(idump-k*iq1)-k*ir1
if (idump<0) idump=idump+im1
k=idum2/iq2
idum2=ia2*(idum2-k*iq2)-k*ir2
if (idum2<0) idum2=idum2+im2
j=1+iy/ndiv
iy=iv(j)-idum2
iv(j)=idump
if (iy<1) iy=iy+imm1
ranp=min(am*iy,rnmx)

return
end function ranp



  !*********************************************************************

  function mat2diag(mat,noff)   result(d)

    !*********************************************************************
    ! Make array with noff'th diagonal = values in 1d vector d (all
    ! other elements 0)
    !*********************************************************************

    use io_tools

    real,dimension(:,:),intent(in)         :: mat
    integer,intent(in)                     :: noff
    real,dimension(size(mat,1)-abs(noff))  :: d

    integer                                :: nmax, n

    if (size(mat,1)/=size(mat,2)) &
         call Message('Err:mat2diag: non-square array',fatal='y')

    nmax = size(mat,1)-abs(noff)

    if (noff>=0) then
       d = (/ (mat(n,n+noff), n=1,nmax) /)
    else
       d = (/ (mat(n-noff,n), n=1,nmax) /)
    endif

  end function mat2diag

  !*********************************************************************

  function diag2mat(d,noff)   result(mat)

    !*********************************************************************
    ! Make array with noff'th diagonal = values in 1d vector d (all
    ! other elements 0)
    !*********************************************************************

    real,dimension(:),intent(in)  :: d
    integer,intent(in)            :: noff
    real,dimension(size(d)+abs(noff),size(d)+abs(noff)) :: mat

    integer                       :: nmax, n

    nmax = size(d)+abs(noff)
    mat = 0. 

    if (noff>=0) then
       do n = 1,nmax
          mat(n,n+noff) = d(n)
       enddo
    else
       do n = 1,nmax
          mat(n-noff,n) = d(n)
       enddo
    endif

  end function diag2mat

  !*********************************************************************

  function tri2mat(tri) result(mat)

    !*********************************************************************
    ! Make a full matrix from a NZ x (-1:1) set of vectors, assuming these
    ! make are the coefficients of a tridiagonal matrix.  
    ! Assumes tri(1,-1) = mat(1,nz) and tri(nz,1) = mat(nz,1)
    !*********************************************************************

    use io_tools

    real,dimension(:,-1:),intent(in)        :: tri
    real,dimension(size(tri,1),size(tri,1)) :: mat
    integer                                 :: nz, n

    nz = size(tri,1)
    if (size(tri,2)/=3) &
       call Message('Error:numerics_lib:size(tri,2)/=3',fatal='y')

!!$    mat1 = 0.
!!$    mat1 = diag2mat(tri(2:nz,-1),-1) 
!!$    call write_field(mat1,'mat1')
!!$    mat1 = mat1 + diag2mat(tri(:,0),0) 
!!$    mat1 = mat1+ diag2mat(tri(1:nz-1,1),1)
!!$    call write_field(mat1,'mat1')

    mat = 0.
    do n = 1,nz-1
       mat(n,n) = tri(n,0)
       mat(n,n+1) = tri(n,1)
       mat(n+1,n) = tri(n+1,-1)
    enddo
    mat(nz,nz) = tri(nz,0)

    if (tri(1,-1)/=0) mat(1,nz) = tri(1,-1)
    if (tri(nz,1)/=0) mat(nz,1) = tri(nz,1)


  end function tri2mat

  !*********************************************************************

  function tridiag(f_in,op) result(f_out)

    !*********************************************************************
    ! Implementation of Numerical Recipies tridiagonal matrix inverter.
    ! Input 'op' is dimesioned op(1:nv,-1:1) and holds the three diagonals
    ! of the operator matrix (op(:,-1) is the subdiagonal, etc...)
    !*********************************************************************
    
    use io_tools

    complex,dimension(:),intent(in) :: f_in
    real,dimension(:,-1:),intent(in):: op
    complex,dimension(size(f_in,1)) :: f_out
    ! Local 
    complex,dimension(size(f_in,1)) :: gamma
    complex                         :: bet
    integer                         :: nv, n

    nv = size(f_in,1)

    if (size(op,1)/=nv) call Message('Error:Tridiag: Operator size not&
                                     & congruent with input field',fatal='y')
    if (size(op,2)/=3) call Message('Error:Tridiag: Dimension 2 of&
                                     & operator must = 3',fatal='y')
    bet = op(1,0)
    f_out(1) = f_in(1)/bet
    do n = 2,nv
       gamma(n) = op(n-1,1)/bet
       bet = op(n,0) - op(n,-1)*gamma(n)
       if (bet==0) stop 'Singular inversion in tridiag'
       f_out(n) = (f_in(n) - op(n,-1)*f_out(n-1))/bet
    enddo
    do n = nv-1,1,-1
       f_out(n) = f_out(n) - gamma(n+1)*f_out(n+1)
    enddo
    
  end function tridiag

  !*********************************************************************

  function tridiag_cyc(f_in,op) result(f_out)

    !*********************************************************************
    ! Implementation of Numerical Recipies CYCLIC tridiagonal matrix inverter.
    ! Input 'op' is dimesioned op(1:nv,-1:1) and holds the three diagonals
    ! of the operator matrix (op(:,-1) is the subdiagonal, etc...)
    ! Additionally, this assumes op(1,-1)/=0, and similarly for op(nv,1).
    ! and assumes that these values are the wrap around coefficients, so that:
    ! f_in(1) = op(1,-1)*f_out(nv) + op(1,0)*f_out(1) + op(1,1)*f_out(2) and
    ! f_in(nv) = op(nv,-1)*f_out(nv-1) + op(nv,0)*f_out(nv) + op(nv,1)*f_out(1)
    !*********************************************************************
    
    use io_tools

    complex,dimension(:),intent(in) :: f_in
    real,dimension(:,-1:),intent(in):: op
    complex,dimension(size(f_in,1)) :: f_out
    ! Local 
    complex,dimension(size(f_in,1))       :: u, z
    real,dimension(size(op,1),-1:size(op,2)-2) :: opnew
    complex                               :: factr
    real                                  :: gamma, alpha, beta
    integer                               :: nv

    nv = size(f_in,1)

    if (size(op,1)/=nv) call Message('Error:Tridiag_cyc: Operator size not&
                                     & congruent with input field',fatal='y')
    if (size(op,2)/=3) call Message('Error:Tridiag_cyc: Dimension 2 of&
                                     & operator must = 3',fatal='y')

    alpha = op(nv,1)
    beta  = op(1,-1)
    gamma = -op(1,0)
    opnew(:,-1) = op(:,-1)
    opnew(:,1)  = op(:,1)
    opnew(2:nv-1,0) = op(2:nv-1,0)
    opnew(1,0)  = op(1,0) - gamma
    opnew(nv,0) = op(nv,0) - alpha*beta/gamma

    f_out = tridiag(f_in,opnew)

    u = 0.
    u(1) = gamma
    u(nv) = alpha

    z = tridiag(u,opnew)
    factr = (f_out(1) + beta*f_out(nv)/gamma)/(1.+z(1)+beta*z(nv)/gamma)

    f_out = f_out - factr*z

  end function tridiag_cyc

  !*********************************************************************

end module numerics_lib
