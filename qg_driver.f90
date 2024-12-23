program qg_driver

! ***********************
! 20/9/2021
! included variables fro Markovian forcing
! force_o (from qg_arrays)
! kf_max,kf_min,forc_coef,forc_corr,norm_forcing (from qg_params)
! Markovian (from qg_run_tools)
! ***********************

  use qg_arrays,       only: q,psi,filter, &
                             rhs, &
                             u_a,u_phi,Setup_fields, &
                             force_o, &
                             xtra,ytra,xtra_1,ytra_1,xtra_0,ytra_0, &
                             uvd_total,u,v,u_1,v_1,u_0,v_0
  use qg_params,       only: kmax, cr,parameters_ok, &
                             pi,dt, forc_coef, use_forcing,&
                             cntr,cnt,total_counts,start, &
                             time, &
                             filter_type,filter_exp,k_cut, &
                             e_o, restarting,psi_init_type, &
                             frame,write_step, &
                             diag2_step,do_spectra,&
                             mean_temperature, &
                             Initialize_parameters,rossby, &
                             kf_max,kf_min,forc_coef,forc_corr,norm_forcing,&
                             cab,cab_1,cab_0,&                                      ! 28/3/2022 for ab3
                             lagrangian,write_lagr,lagr_frame,xtra_file,ytra_file,& ! 28/3/2022 to write Lagrangian trajectories 
                             comp_lag_grad                                          ! 29/3/2022
!                             lagr_time_file

  use qg_run_tools,    only: Get_pv,Invert_pv, &
                             Write_snapshots,get_rhs, &
                             Markovian, &
                             dynlaab,lag_grad                                       ! 29/3/2022
  use qg_init_tools,   only: Init_streamfunction, Init_filter,&
                             Init_u_phi,Init_u_a,&
                             init_traj
  use qg_diagnostics,  only: energy, enstrophy, Get_spectra,diagnostics
  use transform_tools, only: Init_transform
  use io_tools,        only: Message,Write_field

  implicit none
  integer :: icntr
  real  :: mtemp_sortie
  complex,dimension(:,:),allocatable :: rhs_0,rhs_1
  real :: t1,t2

  call Initialize_parameters
  call Message('Initializing model...')

  call Init_transform(kmax)
  call Setup_fields

  ! 24/3/2022
  ! Lagrangian tracers dynamics: 
  ! initialization of particle positions
  if (lagrangian) call init_traj

  ! filtre spatial (dissipation)
  filter = Init_filter(filter_type,filter_exp,k_cut)

  ! rossby number
  if(rossby==cr)then
    call Message('Info: rossby not set - setting to 0')
    rossby=0.
  else
    call Message('Rossby number set to rossby=',r_tag=rossby)
  endif

  psi =Init_streamfunction(psi_init_type,restarting)
! 29/3/2021
  if (rossby>0 .and. restarting) then
     u_phi = Init_u_phi(restarting)
     u_a = Init_u_a(restarting)
  endif

! 25/3/2021
! without filtering the initial condition of a restart
! the first frame written in the restart coincides
! with the last one written in the previous run
  if(.not.(restarting) .and. &
       (.not.(psi_init_type=='read')) .and. &
          (.not.(psi_init_type=='read_grid'))) psi=filter*psi

  if (.not.restarting.and. &
        (.not.(psi_init_type=='read')) .and. &
        (.not.(psi_init_type=='read_grid'))) then
     if (e_o==cr) then
        call Message('Info: e_o not set - setting to 1.')
        e_o=1.
     elseif (e_o<=0) then
        call Message('Error: must have e_o > 0')
        parameters_ok = .false.
     else
        call Message('Initial energy will be set to e_o =', r_tag=e_o)
     endif
     if (energy(psi)>0) then
        psi = sqrt(e_o/energy(psi))*psi
     elseif  (energy(psi)<=epsilon(e_o)) then
        call Message('Error: no amplitude in initial psi field:',&
             r_tag=energy(psi))
        parameters_ok = .false.
     endif
  endif



  allocate(rhs_1(-kmax:kmax,0:kmax))
  allocate(rhs_0(-kmax:kmax,0:kmax))



! 29/3/2021
  if (.not.restarting) call get_rhs

  if (.not.parameters_ok) then
     call Message('The listed errors pertain to values set in your input')
     call Message('namelist file - correct the entries and try again.', &
                   fatal='y')
  endif


  q = Get_pv(psi)

  ! *********** Main time loop *************************

  call Message('Beginning calculation')
  call cpu_time(t1)

  call get_rhs

  ! 24/3/2022
  rhs_0=0.
  rhs_1=0.

! ***************

  do icntr = cnt, total_counts

     cntr=icntr
     start = (cntr==1)

     ! *********************
     ! 24/3/2022 for 2nd step of AB3
     ! cab, cab_1, cab_0 in qg_params.f90 (stessi per Eul. e Lagr.)
     if (cntr==cnt) then
        cab=1.
        cab_1=0.
        cab_0=0.
     elseif (cntr==cnt+1) then
        cab=1.5
        cab_1=0.5
        cab_0=0.
     else
        cab=23./12.
        cab_1=16./12.
        cab_0=5./12.
     endif
     ! *********************

     ! *************************
     ! 29/3/2022
     if (mod(cntr,diag2_step)==0.or.start) then
        call diagnostics
     endif

     ! *************************
     ! 29/3/2022
     if (mod(cntr,write_step)==0.or.start) then
        call Message('mean_temperature:',r_tag=mean_temperature)
        frame = Write_snapshots(frame)
     endif


     ! 28/3/2022
     ! Write trajectories
     if (lagrangian) then
        if (mod(cntr,write_lagr)==0.or.start)then
           lagr_frame=lagr_frame+1
           call Write_field(xtra,xtra_file,lagr_frame)
           call Write_field(ytra,ytra_file,lagr_frame)
!           call Write_field(time,lagr_time_file,lagr_frame)
           call Message('Wrote Lagrangian frame',tag=lagr_frame)
           call Message('time=',r_tag=time)
           ! 29/3/2022
           ! compute gradients (divergence,vorticity)
           ! at particle positions
           if (comp_lag_grad) then 
              call lag_grad(psi,u_phi,u_a)
           endif
        endif
     endif


! advance in time
     ! ********************************
     ! 20/9/2021 Markovian forcing [step with sqrt(delta_t)] /  seems ok (convergence vs delta_t)
     ! to use with forc_corr=0 (delta-correlated in time)
     if (use_forcing) then 
         q=q+sqrt(dt)*Markovian(kf_min,kf_max,forc_coef,forc_corr,force_o,norm_forcing,psi)
         if ((mod(cntr,write_step)==0).and.(cntr>0)) then 
            call Message('gen enstrophy input:', & 
                 r_tag=abs(real(-sum((conjg(q)*force_o)))))
            call Message('gen energy input:', & 
                 r_tag=abs(real(-sum((conjg(psi)*force_o)))))
         endif
     endif
! ********************************
      ! 24/3/2022 AB3
      q=filter*( q + &
           dt*( cab*rhs  - cab_1*(filter*rhs_1)  + cab_0*(filter*filter*rhs_0) ))

     ! 29/3/2022 dynlaab
     ! u,v computed in Get_rhs
     ! either geostrophic or full velocity
     ! selected with flag advect_g
     if (lagrangian) then
        call dynlaab
        u_0=u_1
        v_0=v_1
        u_1=u
        v_1=v
        xtra_0=xtra_1
        ytra_0=ytra_1
        xtra_1=xtra
        ytra_1=ytra
     endif

     rhs_0=rhs_1
     rhs_1=rhs

     psi = Invert_pv(q)
     call get_rhs

     time = time + dt

  enddo  ! End of main time loop

  call cpu_time(t2)
  call Message('Calculation done')
  call Message('total cpu time: ',r_tag=t2-t1)

  !******************* END OF PROGRAM *********************************

end program qg_driver
