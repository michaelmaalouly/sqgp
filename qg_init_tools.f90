module qg_init_tools         !-*-f90-*-  <- tells emacs: use f90-mode

  !************************************************************************
  ! Contains all the routines for initializing fields
  ! IO done here is for reading input fields, error msgs, and writing
  ! preset fields.
  !
  ! Routines:   Init_topo, Init_tracer,
  !            Init_streamfunction, Init_filter, Normu
  !
  ! Dependencies:  qg_transform_tools, io_tools, qg_arrays, qg_params,
  !                qg_strat_tools,qg_run_tools, op_rules, numerics_lib
  !************************************************************************

  implicit none
  private
  save

  public :: Init_streamfunction,Init_filter,&
            Init_u_phi,Init_u_a,&
            init_traj

contains

  !************************************************************************
! 28/3/2022
! initialization of particle positions
  subroutine init_traj
    use io_tools, only: Message, Read_field
    use qg_params, only: pi,nx,ny,ntra,delta0,idump
    use qg_params, only: xtra_file,ytra_file,lagr_init_type,xtra_init_file,ytra_init_file,&
            ci,Parameters_ok,lagr_frame,lagr_start_frame,restarting,&
            xtra_restart_file,ytra_restart_file
    use qg_arrays, only: xtra,ytra,xtra_1,ytra_1,xtra_0,ytra_0
    use numerics_lib, only: ranp
    integer :: ib,jb,mt,mtra
    real :: dxtra,r,thetatra


    restart: if (restarting) then

      if (xtra_restart_file=='') then
          call Read_field(xtra,xtra_file,lagr_frame+1)
      else
          call Read_field(xtra,xtra_restart_file,lagr_start_frame+1)
      endif
      if (ytra_restart_file=='') then
          call Read_field(ytra,ytra_file,lagr_frame+1)
      else
          call Read_field(ytra,ytra_restart_file,lagr_start_frame+1)
      endif

    else

       call Message("Lagrangian tracers initial distribution:")

       select case(trim(lagr_init_type))

       case('uniform')
          ! mtra: number of particles per line
          r=sqrt(float(ntra))
          mtra=int(r)
          if(mod(mtra,2)/=0) then
            call Message("mtra=sqrt(ntra) has to be even")
            stop
          endif
          ! initial particle spacing
          dxtra=2.*pi/mtra
          ! initial particle distribution:
          ! mtra particles per line
          ! particle counter
          mt=0
          do jb=1,mtra
             do ib=1,mtra
                mt=mt+1
                xtra(mt)=(ib-1)*dxtra
                ytra(mt)=(jb-1)*dxtra
!                write(1,*) jb,ib,xtra(mt),ytra(mt)
!                call flush(1)
             enddo
          enddo
          call Message("uniform distribution")
          call Message("number of particles per line:",tag=mtra)
          call Message("initial particle spacing:",r_tag=dxtra)
          xtra_1=0.
          ytra_1=0.
          xtra_0=0.
          ytra_0=0.

       case('uniform_pairs_xy')
          ! mtra: number of pairs per line
          r=sqrt(float(ntra)/3.)
          mtra=int(r)
          if(mod(mtra,2)/=0) then
            call Message("mtra=sqrt(ntra/3) has to be even")
            stop
          endif
          ! initial triplet spacing
          dxtra=2.*pi/mtra
          ! initial particle distribution:
          ! mtra particles per line
          ! particle counter
          mt=0
          do jb=1,mtra
             do ib=1,mtra
                mt=mt+1
                xtra(mt)=(ib-1)*dxtra
                ytra(mt)=(jb-1)*dxtra
                mt=mt+1
                xtra(mt)=xtra(mt-1)+delta0
                ytra(mt)=ytra(mt-1)
                mt=mt+1
                xtra(mt)=xtra(mt-2)
                ytra(mt)=ytra(mt-2)+delta0
!                write(1,*) iz,jb,ib,xtra(mt),ytra(mt)
!                call flush(1)
             enddo
          enddo
          call Message("uniform distribution of pairs")
          call Message("number of triplets per line:",tag=mtra)
          call Message("initial spacing between triplets:",r_tag=dxtra)
          call Message("1 triplet: 1 pair along x and 1 pair along y")
          call Message("initial particle spacing:",r_tag=delta0)
          xtra_1=0.
          ytra_1=0.
          xtra_0=0.
          ytra_0=0.


! ***********************************************          
! 20/4/2022
       case('random')
          if(mod(ntra,2)/=0) then
            call Message("ntra has to be even")
            stop
          endif
          do ib=1,ntra-1,2
             xtra(ib)=2.*pi*ranp(idump)
             ytra(ib)=2.*pi*ranp(idump)
             thetatra=2.*pi*ranp(idump)
             xtra(ib+1)=xtra(ib)+delta0*cos(thetatra)
             ytra(ib+1)=ytra(ib)+delta0*sin(thetatra)
          enddo
          call Message("random distribution")
          call Message("number of particles:",tag=ntra)
          call Message("initial particle spacing:",r_tag=delta0)
          xtra_1=0.
          ytra_1=0.
          xtra_0=0.
          ytra_0=0.
! ***********************************************          


       case('read')
          if (trim(xtra_init_file)=='' .or. trim(ytra_init_file)=='') then
             call Message('Error: no input file for xtra ytra given')
             Parameters_ok=.false.
          endif
          call Message('Initial trajectories will be read from: '&
               &//trim(xtra_init_file)//' and '//trim(ytra_init_file)//&
               ', lagr. frame:', tag=lagr_start_frame)
          if (lagr_start_frame==ci) then
             call Message('Warning: lagr_start_frame & 
                  not initialized-setting to 1')
             lagr_start_frame = 1
          elseif (lagr_start_frame<=0) then
             call Message('Error: require lagr_start_frame>=0')
             Parameters_ok=.false.
          endif
          call Read_field(xtra,xtra_init_file,&
               frame=lagr_start_frame,exclude_dd=1)
          call Read_field(ytra,ytra_init_file,&
               frame=lagr_start_frame,exclude_dd=1)
          xtra_1=0.
          ytra_1=0.
          xtra_0=0.
          ytra_0=0.

       end select

    end if restart

  end subroutine init_traj

  !************************************************************************

  function Init_filter(filter_type,filter_exp,k_cut) result(filter)

    !************************************************************************
    ! Set dealiasing mask for isotropic truncation (semicircle is just tangent
    ! to line in '4/3' rule) and combine it with small scale spatial filter.
    ! Parameter 'filtdec' below is set so that filter decays to
    ! (1+4*pi/nx)**(-1) at kmax
    ! Initialize small scale filter/de-aliasing mask.  Legal filter types are:
    !
    !   hyperviscous : equivalent to RHS dissipation
    !                  nu*del^(2*filter_exp)*field, with nu set optimally.
    !                  Requires 'filter_exp'
    !
    !   exp_cutoff   : exponential cutoff filter (see code below)
    !                  Requires 'filter_exp' and 'k_cut'
    !
    !   none         : none
    !************************************************************************

    use io_tools,  only: Message
    use qg_params, only: kmax,nx,pi,parameters_ok,cr,nu, &
                         dt                                 ! 18/10/2021
    use qg_arrays, only: ksqd_

    character(*),intent(in)            :: filter_type
    real,intent(in)                    :: filter_exp,k_cut
    real,dimension(-kmax:kmax,0:kmax)  :: filter
    ! Local
    real                               :: filtdec
    real                               :: k_al

    filter = 1.                ! Set up de-aliasing part of filter first


    select case (trim(filter_type))
    case ('hyperviscous')

       call Message('Using hyperviscous filter')
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('...with (del**2)**',r_tag=filter_exp)
          call Message('...and with nu = ',r_tag=nu)
       endif

       where (filter>0.)filter=1./(1.+nu*(4.*pi/nx) &
            *(ksqd_/kmax**2)**filter_exp)

! *****************************************
    ! 18/10/2021
    case ('hyperviscous_sb1')

       call Message('Using hyperviscous sb1 filter')
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('...with (del**2)**',r_tag=filter_exp)
          call Message('...and with nu = ',r_tag=nu)
       endif

       where (filter>0.)filter=1./(1.+0.5*dt*nu &
            *ksqd_**filter_exp)

    case ('hyperviscous_sb2')

       call Message('Using hyperviscous sb2 filter')
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('...with (del**2)**',r_tag=filter_exp)
          call Message('...and with prefactor = ',r_tag=nu)
       endif

       where (filter>0.)filter=1./(1.+nu*0.5 &
            *(pi*ksqd_/kmax**2)**filter_exp)


    case ('hyperviscous_ab3')

       call Message('Using hyperviscous filter for Adams Bashforth')
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('...with (del**2)**',r_tag=filter_exp)
          call Message('...and with prefactor = ',r_tag=nu)
       endif

       filter=exp(-dt*nu*(ksqd_/kmax**2)**filter_exp)

! *****************************************

    case ('exp_cutoff')

       call Message('Using exponential cutoff filter')
       call Message('version modifiee pour SQG')
       if (k_cut==cr) then
          call Message('Error: cutoff scale k_cut not set')
          parameters_ok=.false.
       else
          call Message('Cutoff scale k_cut =',r_tag=k_cut)
       endif
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('Filter exponent filter_exp =',r_tag=filter_exp)
       endif
       k_al=sqrt(8./9.)*(kmax+1)
!       k_al=kmax

! Shafer met **2 ???
       filtdec = -log(1+4*pi/nx**1)/(k_al-k_cut)**filter_exp
       where (ksqd_>k_cut**2)
          filter = exp(filtdec*(sqrt(ksqd_)-k_cut)**filter_exp)
       end where

    case ('exp_cutoff_p')

       call Message('Using exponential cutoff filter_p')
       if (k_cut==cr) then
          call Message('Error: cutoff scale k_cut not set')
          parameters_ok=.false.
       else
          call Message('Cutoff scale k_cut =',r_tag=k_cut)
       endif
       if (filter_exp==cr) then
          call Message('Error: filter_exp not set')
          parameters_ok=.false.
       else
          call Message('Filter exponent filter_exp =',r_tag=filter_exp)
       endif
       Call Message('Dimensionless nu parameter =',r_tag=nu)

       k_al=sqrt(8./9.)*(kmax+1)
!       k_al=kmax

! this filter should work well to suppress the Gibbs phenomenon when computing
! the precipitation field
       filtdec = -log(1+4*pi/nx*nu)/(k_al-k_cut)**filter_exp
       where (ksqd_>k_cut**2)
          filter = exp(filtdec*(sqrt(ksqd_)-k_cut)**filter_exp)
       end where

! exemple
!       filtdec = -log(1+4*pi/nx*5)/(k_al-k_cut)**3
!       where (ksqd_>k_cut**2)
!          filter = filter*exp(filtdec*(sqrt(ksqd_)-k_cut)**3)
!       end where


    case ('none')

       call Message('No spatial filtering')

    case default ! or filter_type = 'none' -- make that the default in decln.

       call Message('Error: Must select filter_type.  Legal choices are &
                     &filter_type = hyperviscous|exp_cutoff|none')
       parameters_ok=.false.

    end select

! dealiasing
    where (ksqd_>= (8./9.)*(kmax+1)**2) filter = 0.

    filter(-kmax:0,0) = 0.

  end function Init_filter



  !************************************************************************


  function Init_streamfunction(psi_init_type,restarting) result(psi)

    !************************************************************************
    ! Read in or create initial streamfunction field in manner specified
    !   by variable 'psi_init_type', which can have the values:
    !
    ! spectral_m :  Gaussian spread of initial energy about isotropic horiz.
    !               wavenumber 'k_o' and with width 'delk', and all energy
    !               in vertical mode 'm_o'
    ! spectral_z :  Same as above in horizontal plane, but with all initial
    !               energy in level 'z_o'
    ! elliptical_vortex :  Elliptical gaussian bump in initial grid vorticity
    !               field, aspect ratio 'aspect_vort' and width 'del_vort',
    !               and contained in mode 'm_o'
    ! read :        Read in from 'psi_init_file' at frame 'start_frame'
    !
    ! All of the values in quotes can be set in input namelist
    !************************************************************************

    use io_tools,        only: Message, Read_field
    use qg_params,       only: kmax,z_o,k_o,m_o,delk,i,pi,idum,nkx,nky,&
                               psi_init_file,start_frame,nx,ny,aspect_vort,&
                               del_vort,parameters_ok,cr,ci,psi_file,psi_restart_file,frame,&
                               uag_file,uphi_file,&
                               rossby,isign
    use qg_arrays,       only: ksqd_,kx_,ky_
    use numerics_lib,    only: Ran
    use transform_tools, only: Grid2spec,Spec2grid_cc,ir_prod

    character(*),intent(inout)                 :: psi_init_type
    logical, intent(in)                        :: restarting
    complex,dimension(-kmax:kmax,0:kmax)    :: psi
    real,dimension(:,:),allocatable            :: espec, mu
    real,dimension(:,:),allocatable          :: zetag
    complex,dimension(:,:),allocatable       :: zeta
    real                                       :: radius2,x_0,y_0
    integer                                    :: ix, iy, midx, midy
    complex,dimension(-kmax:kmax,0:kmax,2) :: toto2,AA
    complex,dimension(-kmax:kmax,0:kmax,4) :: toto4
    complex,dimension(nx,ny,4)             :: totog4
    complex,dimension(nx,ny,2)             :: qxg, qyg


    if (trim(psi_init_type)=='spectral') psi_init_type = 'spectral_m'

    restart: if (restarting) then
! ************************************************
! 29/3/2021
       !call Read_field(psi,psi_file,frame+1)
       if (psi_restart_file=='') then
           call Read_field(psi,psi_file,frame+1)
       else
           call Read_field(psi,psi_restart_file,start_frame+1)
       endif
! ************************************************

    else

       psi = 0.
       select case (trim(psi_init_type))

       !********************
       ! 8/12/2021
       ! all scales with apmplitude uniform in k (and random phases)
       case('all_scales_u')

          call Message('Initial condition: random phases uniform amplitude in k')

          ! to have a flat kinetic energy spectrum
          psi(:,:) = ksqd_**(-3./4.)*exp(i*2*pi*Ran(idum,nkx,nky))
       
       !********************
       ! 8/12/2021
       ! from case('all_scales'), available in qg_init_tools.f90.old
       ! initial condition as in Hakim et al., JAS (2002)
       ! default values: k_o=14, m_o=25
       case('spectral_h')

          call Message('Initial condition: streamfunction with random phases, & 
               k.e. spectrum peak at k_o')

          if (k_o==cr) then
             call Message('Error: k_o must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (k_o<=0) then
             call Message('Error: require k_o > 0 - yours is:', r_tag=k_o)
             Parameters_ok=.false.
          else
             call Message('Initial energy centroid at isotropic wavenumber &
                  &k_o =', r_tag=k_o)
          endif
          if (m_o==cr) then
             call Message('Error: exponent m_o must be set to make streamfunction')
             Parameters_ok=.false.
          else
             call Message('Exponent m_o =', tag=m_o)
          endif

          psi(:,:) = sqrt(ksqd_)**(0.25*m_o-1.)/(sqrt(ksqd_)+k_o)**(0.5*m_o) &
               *exp(i*2*pi*Ran(idum,nkx,nky))

          if(rossby>0.)then
            ! calcul de la correction de temperature ageostrophique
            !   sans correction de vent ageostrophique

            ! (u d_u_dz, v, d_v_dz)
            toto4(:,:,1)=-i*ky_*psi
            toto4(:,:,2)=i*((ky_*sqrt(ksqd_))*psi)
            toto4(:,:,3)=i*kx_*psi
            toto4(:,:,4)=-i*((kx_*sqrt(ksqd_))*psi)
            totog4(:,:,1:4)=Spec2grid_cc(toto4(:,:,1:4))

            !  theta, d_theta_dz
            toto2(:,:,1)=-sqrt(ksqd_)*psi
            toto2(:,:,2)=ksqd_*psi
            qxg = Spec2grid_cc(toto2)

            ! -ug theta - F[ theta d_u_dz + ug d_theta_dz ]
            AA(:,:,1)= - Grid2spec(ir_prod(totog4(:,:,1),qxg(:,:,1)) ) &
               - Grid2spec( ir_prod(totog4(:,:,2),qxg(:,:,1)) &
               + ir_prod(totog4(:,:,1),qxg(:,:,2)) )/sqrt(ksqd_)

               AA(:,:,2)= - Grid2spec(ir_prod(totog4(:,:,3),qxg(:,:,1))) &
               - Grid2spec(ir_prod(totog4(:,:,4),qxg(:,:,1)) &
               +ir_prod(totog4(:,:,3),qxg(:,:,2)))/sqrt(ksqd_)

               qyg(:,:,1)=ir_prod(qxg(:,:,1),qxg(:,:,2))
               qyg(:,:,2)=ir_prod(qxg(:,:,1),qxg(:,:,1))

               ! *********************
               ! 18/12/2021
               !toto2=Grid2spec(qyg) ! error on ciclad (divide by zero)
               toto2(:,:,1)=Grid2spec(qyg(:,:,1))
               toto2(:,:,2)=Grid2spec(qyg(:,:,2))
               ! *********************

               psi=psi &
               -rossby*(toto2(:,:,1)+0.5*sqrt(ksqd_)*toto2(:,:,2) &
               +i*kx_*AA(:,:,2)-i*ky_*AA(:,:,1))/sqrt(ksqd_)

          endif

       !********************


       case('spectral_m')

          call Message('Initial streamfunction will be spectrally local')
          if (k_o==cr) then
             call Message('Error: k_o must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (k_o<=0) then
             call Message('Error: require k_o > 0 - yours is:', r_tag=k_o)
             Parameters_ok=.false.
          else
             call Message('Initial energy centroid at isotropic wavenumber &
                  &k_o =', r_tag=k_o)
          endif
          if (delk==cr) then
             call Message('Error: delk must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (delk<=0) then
             call Message('Error: need delk>0 - yours is:',r_tag=delk)
             Parameters_ok=.false.
          else
             call Message('Initial energy peak wavenumber width delk =',&
                  r_tag=delk)
          endif
          !m_o = 0 ! 30/3/2021 because Init_streamfuction called after setting
                   ! m_o in qg_driver

          allocate(espec(-kmax:kmax,0:kmax),mu(-kmax:kmax,0:kmax))
          espec = 1.
          if (delk/=0) espec = exp(-(sqrt(ksqd_)-k_o)**2/delk**2)
          mu = sqrt(ksqd_)        ! Total geostrophic wavenumber
          psi(:,:) = sqrt(espec)/mu*exp(i*2*pi*Ran(idum,nkx,nky))
          deallocate(espec,mu)

       case('spectral_z')

          call Message('Initial streamfunction will be spectral by layer')
          if (k_o==cr) then
             call Message('Error: k_o must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (k_o<=0) then
             call Message('Error: require k_o > 0 - yours is:', r_tag=k_o)
             Parameters_ok=.false.
          else
             call Message('Initial energy centroid at isotropic wavenumber &
                  &k_o =', r_tag=k_o)
          endif
          if (delk==cr) then
             call Message('Error: delk must be set to make streamfunction')
             Parameters_ok=.false.
          elseif (delk<0) then
             call Message('Error: require delk>=0 - yours is:',r_tag=delk)
             Parameters_ok=.false.
          else
             call Message('Initial energy peak wavenumber width delk =',&
                  r_tag=delk)
          endif

          z_o = 1

          allocate(espec(-kmax:kmax,0:kmax))
          espec = 1.
          if (delk/=0) espec = exp(-(sqrt(ksqd_)-k_o)**2/delk)
          psi(:,:) = sqrt(espec)/sqrt(ksqd_)*cexp(i*2*pi*Ran(idum,nkx,nky))
          deallocate(espec)

       case('elliptical_vortex')

          call Message('Initial vorticity field will be an elliptical vortex')
          if (del_vort==cr) then
             call Message('Error: del_vort must be set to make streamfunction')
             Parameters_ok=.false.
          else
             call Message('Initial vortex width del_vort =',r_tag=del_vort)
          endif
          if (aspect_vort==cr) then
             call Message('Error: aspect_vort must be set to make streamfuncn')
             Parameters_ok=.false.
          else
             call Message('Initial vortex aspect ratio aspect_vort =', &
                  r_tag=aspect_vort)
          endif
          m_o = 0


          allocate(zetag(nx,ny),zeta(-kmax:kmax,0:kmax))
          midx = kmax+1; midy = midx; zetag = 0.
          do iy = 1,ny
             do ix = 1,nx
                radius2 = ((ix-midx)**2+aspect_vort*(iy-midy)**2)
                radius2 = radius2*((2*pi)**2/(nx*ny))
                zetag(ix,iy) = exp(-radius2/del_vort**2)
             enddo
          enddo
          zeta = Grid2spec(zetag)
          psi = -zeta*(1./ksqd_)
          deallocate(zetag,zeta)




      case('dipole')


        allocate(zetag(nx,ny),zeta(-kmax:kmax,0:kmax))
    !          allocate(psig_0(nx,ny,1))

        midx = kmax+1; midy =kmax+1;
        zetag = 0.; zeta=0.

        call Message('Initial vorticity field will be two vortices')
        if (del_vort==cr) then
          call Message('Error: del_vort must be set to make streamfunction')
          Parameters_ok=.false.
        else
          call Message('Initial vortex width del_vort =',r_tag=del_vort)
        endif
        if (delk==cr) then
          call Message('Error: delk must be set to make streamfunction')
          Parameters_ok=.false.
        elseif (delk<0) then
          call Message('Error: require delk>=0 - yours is:',r_tag=delk)
          Parameters_ok=.false.
        else
          call Message('Initial distance between vortices delk =',&
                            r_tag=delk)
        endif
        if (isign==cr) then
          call Message('Error: isign must be set to make streamfunction')
          Parameters_ok=.false.
        else
          call Message('Initial sign of second vortex isign =',&
                            r_tag=isign)
        endif

        do iy = 1,ny
           do ix = 1,nx
              x_0=(ix-midx)*2*pi/nx-delk/2.
              y_0=(iy-midy)*2*pi/nx
              if(sqrt(x_0**2+y_0**2)<=del_vort) &
                       zetag(iy,ix)=1.-(x_0**2+y_0**2)/del_vort**2
          !                if(sqrt(x_0**2+y_0**2)<=del_vort)zetag(iy,ix)=1.

              x_0=(ix-midx)*2*pi/nx+delk/2.
              y_0=(iy-midy)*2*pi/nx
              if(sqrt(x_0**2+y_0**2)<=del_vort) &
                     zetag(iy,ix)=isign*(1.-(x_0**2+y_0**2)/del_vort**2)


          enddo
        enddo

        zeta = Grid2spec(zetag)
        psi = -zeta*(1./(ksqd_))
        deallocate(zetag,zeta)

     case('interaction')

        call Message('Initial configuration with 3 eddies interacting')
        if (delk==cr) then
          call Message('Error: delk must be set to make streamfunction')
          Parameters_ok=.false.
        else
          call Message('Initial dipole distance delk =',r_tag=delk)
        endif
        if (del_vort==cr) then
          call Message('Error: del_vort must be set to make streamfunction')
          Parameters_ok=.false.
        else
          call Message('Initial eddy radius del_vort =',r_tag=del_vort)
        endif

        ! trois tourbillons (papier JMR)
        allocate(zetag(nx,ny),zeta(-kmax:kmax,0:kmax))
        midx = kmax+1; midy =(kmax+1)/4-1; zetag = 0.
        do iy = 1,ny
          do ix = 1,nx

            ! 1er tourbillon
            x_0=(ix-midx)*2*pi/nx-delk/2.
            y_0=(iy-midy)*2*pi/nx
            if(sqrt(x_0**2+y_0**2)<=del_vort) &
                     zetag(iy,ix)=(1.-(x_0**2+y_0**2)/del_vort**2)

            ! 2eme tourbillon
            x_0=(ix-midx)*2*pi/nx+delk/2.
            y_0=(iy-midy)*2*pi/nx
            if(sqrt(x_0**2+y_0**2)<=del_vort) &
                     zetag(iy,ix)=-(1.-(x_0**2+y_0**2)/del_vort**2)


            ! 3eme tourbillon
            x_0=(ix-midx)*2*pi/nx+delk/2.*1.2
            y_0=(float(iy)-float(midx)*1.5)*2*pi/float(nx)
            if(sqrt(x_0**2+y_0**2)<=del_vort) &
                 zetag(iy,ix)=isign*(1.-(x_0**2+y_0**2)/del_vort**2)/2.

          enddo
        enddo

        zeta = Grid2spec(zetag)
        psi = 0.
        psi = -zeta*(1./sqrt(ksqd_))

        deallocate(zetag,zeta)

       case('read')

          if (trim(psi_init_file)=='') then
             call Message('Error: no input file for psi given')
             Parameters_ok=.false.
          endif
          call Message('Initial streamfunction will be read from: '&
               &//trim(psi_init_file)//', frame:', tag=start_frame)
          if (start_frame==ci) then
             call Message('Warning: start_frame not initialized-setting to 1')
             start_frame = 1
          elseif (start_frame<=0) then
             call Message('Error: require start_frame>=0')
             Parameters_ok=.false.
          endif

          call Read_field(psi,psi_init_file,frame=start_frame,exclude_dd=1)

      case('read_grid')

          if (trim(psi_init_file)=='') then
             call Message('Error: no input file for psi given')
             Parameters_ok=.false.
          endif
          call Message('Initial grid tracer will be read from: '&
               &//trim(psi_init_file)//', frame:', tag=start_frame)
          if (start_frame==ci) then
             call Message('Warning: start_frame not initialized-setting to 1')
             start_frame = 1
          elseif (start_frame<=0) then
             call Message('Error: require start_frame>=0')
             Parameters_ok=.false.
          endif
          allocate(zetag(nx,ny),zeta(-kmax:kmax,0:kmax))
          call Read_field(zetag,psi_init_file,frame=start_frame,exclude_dd=1)
          zeta = Grid2spec(zetag)
          psi = -zeta*(1./sqrt(ksqd_))
          deallocate(zetag,zeta)


       case default

          call Message('Error: must select psi_init_type = &
               &read|spectral_m|spectral_z|spectral|elliptical_vortex &
               &|all_scales - yours is:'//trim(psi_init_type))
          Parameters_ok=.false.

       end select

    end if restart

  end function Init_streamfunction

  !*********************************************************************
! 29/3/2021
  function Init_u_phi(restarting) result(u_phi)
    use io_tools,        only: Message, Read_field
    use qg_params,       only: kmax,uphi_file,frame,uphi_restart_file,start_frame

    complex,dimension(-kmax:kmax,0:kmax,2)   :: u_phi
    logical, intent(in)                      :: restarting

! ************************************************
! 29/3/2021
       if (uphi_restart_file=='') then
           call Read_field(u_phi,uphi_file,frame+1)
       else
           call Read_field(u_phi,uphi_restart_file,start_frame+1)
       endif
! ************************************************

  end function Init_u_phi

  !*********************************************************************
! 29/3/2021
  function Init_u_a(restarting) result(u_a)
    use io_tools,        only: Message, Read_field
    use qg_params,       only: kmax,uag_file,frame,uag_restart_file,start_frame

    complex,dimension(-kmax:kmax,0:kmax,2)   :: u_a
    logical, intent(in)                      :: restarting

! ************************************************
! 29/3/2021
       if (uag_restart_file=='') then
           call Read_field(u_a,uag_file,frame+1)
       else
           call Read_field(u_a,uag_restart_file,start_frame+1)
       endif
! ************************************************

  end function Init_u_a



    !************************************************************************
    ! Read in or create initial streamfunction field in manner specified
end module qg_init_tools
