module qg_diagnostics

  !************************************************************************
  ! Energetics and other diagnostics
  !
  ! Routines:  Get_energetics, Get_spectra
  !
  ! Dependencies: qg_arrays, qg_params, io_tools, op_rules, transform_tools
  !
  !************************************************************************

  implicit none
  private
  save

  public :: Get_spectra, energy, enstrophy,diagnostics

contains

  !*************************************************************************
  subroutine diagnostics
    use qg_params,       only: d2frame,do_spectra,rossby
    use qg_arrays,       only: q,u_phi,u_a
    use io_tools,        only: Message

    real :: vvv
        ! 8/12/2021
        ! Computation of ke spectra
        if (do_spectra) d2frame = Get_spectra(d2frame)

        vvv=0.
        vvv=sqrt(sum(real(abs(q)**2)))
        call Message('rms ug=',r_tag=vvv)
        if (isnan(sum(real(abs(q))))) call Message('plantage Nan',fatal='y')

        if(rossby>0.) then
          vvv=0.
          vvv=sqrt(sum(real(abs(u_phi(:,:,1))**2 &
             +abs(u_phi(:,:,2))**2)))
          call Message('rms u_phi=',r_tag=vvv)
          vvv=0.
          vvv=sqrt(sum(real(abs(u_a(:,:,1))**2 &
             +abs(u_a(:,:,2))**2)))
          call Message('rms u_a=',r_tag=vvv)
        endif

    return
  end subroutine diagnostics


  !*************************************************************************

  real function energy(psi)

    ! 21/4/2022
    ! Get kinetic energy (geostrophic flow)

    use qg_arrays, only: ksqd_

    complex,dimension(:,:),intent(in) :: psi
    real                                :: ke=0.

    ke = real(sum(ksqd_*(abs(psi)**2)))
    !ke = 0.5*real(sum(ksqd_*(abs(psi)**2)))

    energy = ke

  end function energy

  !*************************************************************************

  real function enstrophy(q)

    ! 21/4/2022
    ! Get enstrophy (geostrophic flow)

    use qg_arrays, only: ksqd_

    complex,dimension(:,:),intent(in) :: q

    enstrophy = real(sum( (ksqd_*(q*conjg(q)) )))
    !enstrophy = 0.5*real(sum( (ksqd_*(q*conjg(q)) )))

  end function enstrophy

  !*************************************************************************

  function Get_spectra(framein) result(dframe)

    !************************************************************************
    ! Calculate the isotropic horizontal wavenumber vs vertical
    ! wavenumber spectra of modal and layered energetics
    !************************************************************************

    use io_tools,        only: Write_field
    use numerics_lib,    only: Ring_integral
    use qg_arrays,       only: ksqd_,kxv,kyv,psi
    use qg_params,       only: time,kmax, &
                               diag2_time_file,kes_file

    integer,intent(in)                     :: framein
    integer                                :: dframe
    real,dimension(:,:),allocatable      :: field
    real,dimension(:),allocatable        :: spec
    dframe = framein + 1
    call Write_field(time,diag2_time_file,dframe) ! Track diagnostic-writes

    allocate(spec(kmax),field(-kmax:kmax,0:kmax))

    ! ******************* 
    ! KE spectra
    ! 21/4/2022 

    field = ksqd_*abs(psi(:,:))**2
    !field = 0.5*ksqd_*abs(psi(:,:))**2       
    spec = Ring_integral(field,kxv,kyv,kmax)
    call Write_field(spec,kes_file,dframe)
    ! ******************* 


    deallocate(field,spec)

  end function Get_spectra

  !*********************************************************************

end module qg_diagnostics
