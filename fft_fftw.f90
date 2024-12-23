module fft_pak                    
! Arrow of time adjusted......8 Nov 04
!
! This one is for the FFTW 3.0.1 fma package, uniprocessor, double precision
! This one tries to use the same plan_f, plan_b all the time. 
!
! Routine: Getfft
!

!  Oliver's notes: Here it is.  There is one explicit kind selection
! for the Fortran integer plan that is used as a C pointer by fftw; this
! must be 32 bit on my G4.  Otherwise this is almost exactly as fftw2,
! so I don't know why they made such a fuss about not making it
! backwards compatible.


! This is from fftw3.f:

  include 'fftw3.f'
  include 'fftw3_mkl.f'



integer*8                :: plan_f,plan_b
!  FFTW needs plan to be pointer on machine.  32bits for G4? Yes.
!integer (kind=3)                :: plan_f,plan_b
integer                  :: hres
integer,dimension(2)     :: nsize

private
save

public :: Init_fft, fft

  interface fft
     module procedure fft2, fft3
  end interface

contains

  !********************************************************************

  subroutine Init_fft(kmax)
  use io_tools,only: Message

    ! Initialize fft routine
    ! Does little for version 3
    integer,intent(in)                     :: kmax
    complex,dimension(:,:), allocatable    :: f,fr
    integer iret,n_cpu

    hres = 2*kmax+2
    nsize(1) = hres
    nsize(2) = hres

    n_cpu=MKL_GET_MAX_THREADS()
    call dfftw_init_threads(iret)
    call Message('iret = ',tag=iret)
    call Message('n_cpu = ',tag=n_cpu)
    call dfftw_plan_with_nthreads(n_cpu)

    allocate(f(nsize(1),nsize(2)))
    allocate(fr(nsize(1),nsize(2)))

    ! FFTW_PATIENT -> tres long
    ! FFT_ESTIMATE
    ! FFT_MEASURE
    f=0.
    fr=0.
    call dfftw_plan_dft_2d(plan_b,nsize(1),nsize(2),f,fr,&
          FFTW_BACKWARD,FFTW_MEASURE) 
    call dfftw_plan_dft_2d(plan_f,nsize(1),nsize(2),f,fr,&
          FFTW_FORWARD,FFTW_MEASURE)
    deallocate(f)
    deallocate(fr)


  end subroutine Init_fft

  !********************************************************************

  function fft2(f,dirn) result(fr)

    ! Calculate 2d complex to complex fft.  dirn = +1 or -1.
    ! these values correspond to sign of exponent in spectral
    ! sum - arrow of time?

  complex,dimension(:,:)                  :: f
  complex,dimension(size(f,1),size(f,2))  :: fr
  integer,intent(in)                      :: dirn
  real                                    :: scale=1.0


  fr=0.

  if (dirn==-1) then
     scale=1.0
     call dfftw_execute_dft(plan_b,f,fr)
  elseif (dirn==1)  then
     scale=1.0/float(hres)**2
     call dfftw_execute_dft(plan_f,f,fr)
  endif
 
  fr = scale*fr




  end function fft2




  function fft3(f,dirn) result(fr)

    ! Calculate 2d complex to complex fft.  dirn = +1 or -1.
    ! these values correspond to sign of exponent in spectral
    ! sum (see man page on zzfft2d).

    complex,dimension(:,:,:)                          :: f
    complex,dimension(size(f,1),size(f,2),size(f,3))  :: fr
    complex,dimension(size(f,1),size(f,2))  :: g,gr
    integer,intent(in)                                :: dirn
    real                                      :: scale


!    print*,'titi'
    fr=0.
    if (dirn==1) scale=1.0/float(nx*ny)
    if (dirn==-1) scale=1.0

    do n = 1,size(f,3)
!       gr=0.
!       g=0.
!       g(1:nx,1:ny)=f(1:nx,1:ny,n)
!       call ccfft2d(dirn,nx,ny,scale,g,nx1,gr,nx1,table,work,0)

       g=f(:,:,n)
       gr=fft2(g,dirn)
       fr(:,:,n)=gr(:,:)

    enddo


  end function fft3


  !********************************************************************

  function fft3bis(f,dirn) result(fr)

    ! Calculate 3d complex to complex fft.  dirn = +1 or -1.
    ! these values correspond to sign of exponent in spectral
    ! sum (see man page on ccfft2d).

    complex,dimension(:,:,:)                          :: f
    complex,dimension(size(f,1),size(f,2),size(f,3))  :: fr
    integer,intent(in)                                :: dirn
    real                                              :: scale=1.0

    do n=1,size(f,3)
      if (dirn==-1) then
         scale=1.0
         call dfftw_execute_dft(plan_b,f(:,:,n),fr(:,:,n))
      elseif (dirn==1)  then
         scale=1.0/float(hres)**2
         call dfftw_execute_dft(plan_f,f(:,:,n),fr(:,:,n))
      endif
    enddo

    fr = scale*fr

  end function fft3bis

  !********************************************************************

end module fft_pak
