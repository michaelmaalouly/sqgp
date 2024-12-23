module qg_params                   !-*-f90-*-

  !************************************************************************
  ! Contains all global parameters for program, and routines for
  ! processing their I/O
  !
  ! Routines: Initialize_parameters, Write_parameters, Check_parameters
  !
  ! Dependencies: IO_tools, Syscalls
  !************************************************************************

  implicit none
  public
  save

  integer,parameter       :: ci=-9      ! For checking init status of variables
  real,parameter          :: cr=-9.

  ! Parameters in namelist input -- initialize some parameterss to
  ! negative values (cr and ci) in order to check they are intentionally
  ! initialized in Parameters_ok (below)

  ! Resolution

  integer                 :: kmax        = ci          ! Spectral resolution

!temps
  real                    :: dt          = cr          ! Time step
  integer                 :: write_step  = ci          ! Frame snapshot step
  integer                 :: diag1_step  = ci          ! Diagnostics 1 step
  integer                 :: diag2_step  = ci          ! Diagnostics 2 step
  integer                 :: total_counts= ci          ! Total timesteps to do
  integer                 :: start_frame = ci          ! Frame to start from
  integer                 :: cntr        = 1           ! Timestep counter value
  integer                 :: frame       = 0           ! Field snapshot frame
  integer                 :: d1frame     = 0           ! Diagnostics 1 frame
  integer                 :: d2frame     = 0           ! Diagnostics 2 frame
  real                    :: time        = 0
  character(8)            :: end_date                  ! End date of sim
  character(10)           :: end_time                  ! End time of sim
! *********************************************
! 28/3/20222
! coefficients for Adams-Bashforth order 3 (values set in qg_driver.f90) 
  real                    :: cab         = 0.           ! coeff of term at t
  real                    :: cab_1       = 0.           ! coeff of term at t-1
  real                    :: cab_0       = 0.           ! coeff of term at t-2
! *********************************************

!parametre du modele
  real                    :: rossby = cr  !Rossby number

! *********************************************
! 28/3/20222
! parameters for Lagrangian dynamics in SQG+1 
  integer                 :: lagr_ifr_run     = 0               ! run number
  integer                 :: ntra             = ci              ! number of trajectories
  integer                 :: lagr_frame       = 0               ! Lagrangian trajecrories frame
  integer                 :: write_lagr       = ci              ! Lagrangian frame snapshot step
  integer                 :: lagr_start_frame = ci              ! Frame to start from (Lagrangian)
  real                    :: delta0      = cr                   ! distance between particles in a pair
  character(20)           :: lagr_init_type  = ''               ! Init lagr type
  character(70)           :: xtra_init_file = ''                ! xtra input field
  character(70)           :: ytra_init_file = ''                ! ytra input field
! *********************************************



!initialisation
  logical                 :: restarting     = .false.     ! Is this a restart?
  logical                 :: linear         = .false.     ! Omit non-linear terms
  logical                 :: lagrangian     = .false.     ! Must be True for Lagrangian dynamics  (28/3/2022)
  logical                 :: advect_g       = .false.     ! Must be True for advection with geostrophic flow  (29/3/2022)
  logical                 :: comp_lag_grad  = .false.     ! gradients at lagrangian positions (29/3/3033)
! ***********
! 25/3/2021
  logical                 :: split_runs  = .false.     ! Must be True to split runs
  integer                 :: ifr_run     = 0           ! run number
! ***********
  character(20)           :: psi_init_type  = ''       ! Init streamfn type
  real                    :: e_o         = cr          ! Initial energy
  real                    :: k_o         = cr          ! Initial peak k
  real                    :: delk        = cr          ! Initial spread in k
  real                    :: aspect_vort = cr          ! Initl vort aspct ratio
  real                    :: del_vort    = cr          ! Initial vortex width
  integer                 :: m_o         = ci          ! Initial modal peak
  integer                 :: z_o         = ci          ! Initial level
  character(70)           :: psi_init_file = ''        ! Psi input field
! ***********
  real                 :: isign = cr
! ***********

  real                    :: mean_temperature   = 0.    ! mean temperature
  real                    :: mean_divtemp = 0.  ! average of div x temperature

!forcage
  logical                 :: use_forcing = .false.     ! Use rndm markov frcing
  logical                 :: norm_forcing= .false.     ! Norm gen rate from RMF
  real                    :: forc_coef   = cr          ! BT forcing coefficient
  real                    :: forc_corr   = cr          ! BT forcing correlation
  real                    :: kf_min      = cr          ! min k^2 for BT frc
  real                    :: kf_max      = cr          ! max k^2 for BT frc

!spectres diagno
  logical                 :: do_spectra  = .true.      ! Calc/save spectra

!dissipation
  real                    :: nu          =  1.         ! viscous parameter
  character(20)           :: filter_type = ''          ! Spatial filter
  real                    :: filter_exp  = cr          ! Filter exponent
  real                    :: k_cut       = cr          ! Exp cutoff scale
  real                    :: therm_drag  = 0.          ! Thermal drag
  integer                 :: n_hypo      = ci          ! hypofriction exponent

  ! DIP Switches and tuning factors - factory settings.
  ! All of these are included in namelist input as well, but they
  ! are not checked for initialization in Parameters_ok so that you
  ! dont have to include them in the input file (but can if you need to).

  integer                 :: recunit     = 8            ! For direct access IO
  integer                 :: idum                       ! Random generator seed
  integer                 :: idump                      ! Random generator seed (particles' init. pos.)

  ! Numerical stability tuning parameters


  real                    :: rmf_norm_min= 1.e-5        ! RMF genn min 4 normn

  ! **************End of namelist input parameters*************
  !
  ! Parameters for global internal use - NOT included in any namelist input
  !
  ! Output file and directory names - diagnostics outputs are defined in
  ! respective diagnostics module.

  character(80)           :: datadir       = '.'
  character(32)           :: inputfile     = 'input.nml'
  character(32)           :: restartfile   = 'restart.nml'
! ******************
! 29/3/2021
  character(32)           :: psi_file              = 'psi'       ! Psi snapshots 
  character(32)           :: psi_restart_file      = ''          ! Psi snapshots for restart (for write on separate files)
  character(32)           :: uag_file              = 'uag'
  character(32)           :: uag_restart_file      = ''
  character(32)           :: uphi_file             = 'uphi'
  character(32)           :: uphi_restart_file     = ''
  character(32)           :: write_time_file       = 'write_time'
  character(32)           :: energy_file           = 'energy'    ! 8/12/2021 
  character(32)           :: enstrophy_file        = 'enstrophy' ! 8/12/2021
! 8/12/2021
  character(32)           :: diag2_time_file       = 'diag2_time'
  character(32)           :: kes_file              = 'kes'
! ******************
! 29/3/2022
! For Lagrangian SQG+1
  character(32)           :: xtra_restart_file     = ''          ! xtra snapshots for restart (for write on separate files) 
  character(32)           :: ytra_restart_file     = ''          ! ytra snapshots for restart (for write on separate files) 
!  character(32)           :: lagr_time_file       = 'lagr_time'
! output files for Lagrangian SQG+1
  character(32)           :: xtra_file          = 'xtra'
  character(32)           :: ytra_file          = 'ytra'
  character(32)           :: Div_tra_file       = 'Div_tra'     ! div at particle positions file
  character(32)           :: Vort_g_tra_file    = 'Vortg_tra'   ! geostrophic vort at particle positions file
  character(32)           :: Vort_phi_tra_file  = 'Vortphi_tra' ! vort phi at particle positions file
  character(32)           :: Vort_a_tra_file    = 'Vorta_tra'   ! vort a at particle positions file
! ******************

  ! Resolution and tuning parameters (set as functions of kmax)

  integer                 :: numk, nkx, nky, nx, ny, nv, nmask

  ! Internal flags and counters

  logical                 :: free_surface = .false.
  logical                 :: parameters_ok = .true.
  logical                 :: start
  integer                 :: cnt = 1, call_q = 0, call_b = 0, call_t = 0

  ! Cabalistic numbers

  real,parameter          :: pi          = 3.14159265358979
  complex,parameter       :: i           = (0.,1.)

  ! Namelist declarations


! declaration variables spatiales
  namelist/run_params/kmax
! ****************
! 25/3/2021
  namelist/run_params/recunit,idum,idump
! ****************

! declaration variables temporelles
  namelist/run_params/dt
  namelist/run_params/total_counts,start_frame,cntr,frame,d1frame,d2frame,time
  namelist/run_params/end_date,end_time
  namelist/run_params/rmf_norm_min

!declaration initialisation
  namelist/run_params/restarting,linear,psi_init_type
! ****************
! 25/3/2021
  namelist/run_params/split_runs,ifr_run
! ****************
  namelist/run_params/k_o,delk,aspect_vort,del_vort,m_o,z_o,e_o,isign
! ****************
! 29/3/2021
  namelist/run_params/psi_init_file,psi_restart_file
  namelist/run_params/uphi_restart_file,uag_restart_file
  namelist/run_params/mean_temperature
! ****************

!declaration parametres du modele
  namelist/run_params/rossby

! *********************************************
! 28/3/2022
! parameters for Lagrangian dynamics
  namelist/run_params/lagrangian,lagr_init_type,ntra,delta0
  namelist/run_params/lagr_ifr_run,lagr_frame,write_lagr,lagr_start_frame
  namelist/run_params/xtra_init_file,ytra_init_file,xtra_restart_file,ytra_restart_file
  namelist/run_params/advect_g,comp_lag_grad
! *********************************************

!declaration forcage
  namelist/run_params/use_forcing,norm_forcing
  namelist/run_params/forc_coef,forc_corr,kf_min,kf_max

!declaration dissipation
  namelist/run_params/filter_type,nu
  namelist/run_params/filter_exp,k_cut
  namelist/run_params/therm_drag
  namelist/run_params/n_hypo

!declaration diagnostics spectres
  namelist/run_params/do_spectra
  namelist/run_params/write_step,diag1_step,diag2_step


  !************** Routines for parameter I/O*********************

contains

  subroutine Initialize_parameters

    !************************************************************************
    ! Read in command line arguments 'datadir' and 'inputfile'
    ! then read namelist in 'inputfile' for input parameters,
    ! and pass some parameters to 'io_tools' which are needed for I/O.
    ! If no first arg supplied, then datadir = ./, and if no
    ! 2nd arg, then inputfile = input.nml (in datadir)
    ! Any I/O errors in this part occur before program
    ! knows where to write 'run.log' to (I/O params not passed
    ! to io_mod yet), so errors are written to screen and to
    ! error.log in executable directory for non-interactive runs.
    !************************************************************************

    use io_tools, only: Message, Pass_params
    use syscalls, only: Get_arg
!    external : getarg
    character(80)         :: fnamein='', temp=''
    integer               :: fin=7, nchars, iock
    logical               :: restart_exists=.false.

    call Get_arg(1,datadir,nchars,iock)
    if (iock/=0) call Message('1st getarg failed:'//trim(datadir),&
         iock,fatal='y')
    if (datadir(nchars:nchars)/='/') datadir=trim(datadir)//'/'
    call Pass_params(datadir,recunit) ! Send to io_mod

    ! If there is a restart file in the data directory, then use
    ! this for parameter input, unless there is a second cmnd line arg,
    ! which allows selection of input nml file.

    inquire(file=trim(datadir)//trim(restartfile),exist=restart_exists)
    if (restart_exists) inputfile = restartfile

    call Get_arg(2,temp,nchars,iock)
    if (iock/=0) then
       call Message('2nd getarg failed:'//trim(temp),tag=iock,fatal='y')
    else
       if (trim(temp)/='') then
          inputfile = temp
       endif
    endif

    fnamein = trim(datadir)//trim(inputfile)

    open(unit=fin,file=fnamein,status='old',delim="apostrophe",iostat=iock)
    if (iock/=0) call Message('Open input namelist error; file:'//fnamein, &
         tag=iock,fatal='y')
    read (unit=fin,nml=run_params,iostat=iock)

    if (iock/=0) call Message('Read input namelist error; file, iock ='&
         &//fnamein, tag=iock,fatal='y')
    close(fin)

    parameters_ok = Check_parameters()  ! See below - check and set counters

    nx = 2*(kmax+1)            ! Physical resolution in x
    ny = 2*(kmax+1)            ! Physical resolution in y
    nkx = 2*kmax+1             ! Number of zonal wavenumbers
    nky = kmax+1               ! Number of meridional wavenumbers
    numk = 2*kmax*(kmax+1)     ! Total number of (k,l) points in horiz plane
    idum = -abs(idum)          ! Make sure random num gen is set to start right
    idump = -abs(idump)        

    cnt = cntr

  end subroutine Initialize_parameters

  !**************************************************************

  subroutine Write_parameters

    !************************************************************************
    ! Write parameters to restart namelist file and selected
    ! params to a special file for use by SQG Matlab routines.
    !************************************************************************

    use io_tools, only: Message, Open_file

    character(70) :: restartname, exitcheckname
    integer       :: outf = 40,iock

    restartname = trim(datadir)//trim(restartfile)

    ! Set certain parameters the way they should be set for restarting,
    ! regardless of what goes on in the program.  Store temp vals as nec.

    restarting = .true.
    start_frame = frame
    if (lagrangian) lagr_start_frame = lagr_frame
    call date_and_time(DATE=end_date,TIME=end_time)

    open (unit=outf,file=restartname,status='unknown',delim="apostrophe",&
         iostat=iock)
    if (iock /= 0) call Message('write params nml open error; &
         &iock =',tag=iock,fatal='y')
    write (unit=outf,nml=run_params,iostat=iock)
    if (iock /= 0) call Message('write params nml write error; &
         &iock =',tag=iock,fatal='y')
    close(outf)

    ! Write params to bin file for reading with MATLAB SQG function getparams.m

    call Open_file(outf,'parameters','unknown',9)
    write(unit=outf,rec=1) kmax,1,d1frame,d2frame,frame
    close(outf)


  end subroutine Write_parameters

  !**************************************************************

  logical function Check_parameters()

    !**************************************************************
    ! This function will test that certain input namelist params are consistent
    ! and print some messages with basic info for run
    !**************************************************************

    use io_tools, only: Message

    Check_parameters=.true.

    call Message('')
! **********************
! 25/3/2021
    call Message('SQG plus 1 model')
! **********************
   if (lagrangian) call Message('Lagrangian tracers dynamics')    ! 28/3/2022
    call Message('')
    call Message('Checking parameters for consistency')
    call Message('')
    call Message('Input file = '//trim(datadir)//inputfile)
    call Message('Output directory = '//trim(datadir))

    ! Resolution

    if (kmax==ci) then
       call Message('Error: kmax not set')
       Check_parameters=.false.
    elseif (mod(kmax+1,2)/=0) then
       call Message('Error: kmax must be odd - yours is: ',tag=kmax)
       Check_parameters=.false.
    else
       call Message('Horizontal spatial resolution =',tag=2*(kmax+1))
    endif

    ! 26/1/2022  
    ! large-scale dissipation
    if (therm_drag==cr)  then
          call Message('Error: therm_drag not set')
          Check_parameters=.false.
       else
          call Message('therm_drag =',r_tag=therm_drag)
          if(n_hypo==ci) then
             call Message('Error: n_hypo not set')
             Check_parameters=.false.
          else
             call Message('n_hypo =', tag=n_hypo)
          endif
    endif

! *******************************************************
    ! 28/3/2022
    ! Check parameters for Lagrangian SQG+1
    if (lagrangian) then
       if (ntra==ci)  then
          call Message('Error: Number of particles per depth must be set for Lagrangian SQG+1')
          Check_parameters=.false.
       else
          call Message('Number of particles per depth =',tag=ntra)
       endif
       if (.not.advect_g) then 
          call Message('Lagrangian advection by full velocity field')
       else
          call Message('Lagrangian advection by geostrophic flow')
       endif
       if (comp_lag_grad) then
          call Message('Velocity gradients at particle positions are computed')
       else
          call Message('Velocity gradients at particle positions are not computed')
       endif
    endif

! *******************************************************
! 25/3/2021
    ! Check if runs must be splitted 
    if (split_runs) then
       write(psi_file,"('psi.',i2.2)") ifr_run
       write(uag_file,"('uag.',i2.2)") ifr_run
       write(uphi_file,"('uphi.',i2.2)") ifr_run
       write(write_time_file,"('write_time.',i2.2)") ifr_run
       write(energy_file,"('energy.',i2.2)") ifr_run
       write(enstrophy_file,"('enstrophy.',i2.2)") ifr_run
       write(diag2_time_file,"('diag2_time.',i2.2)") ifr_run
       write(kes_file,"('kes.',i2.2)") ifr_run
       ! 28/3/2022
       ! For Lagrangian SQG+1
       if (lagrangian) then
          write(xtra_file,"('xtra.',i2.2)") lagr_ifr_run
          write(ytra_file,"('ytra.',i2.2)") lagr_ifr_run
!          write(lagr_time_file,"('lagr_time.',i2.2)") lagr_ifr_run
          if (comp_lag_grad) then
             if (.not.advect_g) then  
                write(Div_tra_file,"('Div_tra.',i2.2)") lagr_ifr_run
                write(Vort_phi_tra_file,"('Vortphi_tra.',i2.2)") lagr_ifr_run
                write(Vort_a_tra_file,"('Vorta_tra.',i2.2)") lagr_ifr_run
             endif
             write(Vort_g_tra_file,"('Vortg_tra.',i2.2)") lagr_ifr_run
          endif
       endif
    endif

! *******************************************************
    ! Check whether restarting, and if not, check psi init type

    restart: if (restarting) then

       call Message('This is a restart run')
       if (start_frame==ci) then
          call Message('Error: start_frame must be set when restarting')
          Check_parameters=.false.
       elseif (start_frame<0.or.start_frame>frame) then
          call Message('Error: require 0 <= start_frame <= frame')
          Check_parameters=.false.
       elseif (start_frame<=frame) then      ! Set counters for restart
! 25/3/2021
!          frame = start_frame
          ! Check if runs must be splitted 
          if (.not.split_runs) then
             frame = start_frame
             if (lagrangian) lagr_frame = lagr_start_frame !28/3/2022
             if (frame==0) then
                cntr = 1
                d1frame = 0
                d2frame = 0
             else
               ! ******************************        
                ! 25/3/2021
                ! just in case write_step=1
                if (write_step==1) then
                   cntr = (frame+1)*write_step
                else
                   cntr = frame*write_step
                endif
                ! ******************************    
                d1frame = (cntr-1-mod(cntr-1,diag1_step))/diag1_step + 1
                d2frame = (cntr-1-mod(cntr-1,diag2_step))/diag2_step + 1
             endif
          else
          ! 25/3/2021 
          ! runs are splitted 
             frame=0
             ! ******************************        
             ! just in case write_step=1
             if (write_step==1 .and. ifr_run==0) then
                cntr = (start_frame+1)*write_step*(ifr_run+1)
             else
                cntr = start_frame*write_step*(ifr_run+1)
             endif
             ! ******************************        
             d1frame = 0
             d2frame = 0
             psi_restart_file = psi_file
             uphi_restart_file = uphi_file
             uag_restart_file = uag_file
             write(psi_file,"('psi.',i2.2)") ifr_run+1
             write(uag_file,"('uag.',i2.2)") ifr_run+1
             write(uphi_file,"('uphi.',i2.2)") ifr_run+1
             write(write_time_file,"('write_time.',i2.2)") ifr_run+1
             write(energy_file,"('energy.',i2.2)") ifr_run+1
             write(enstrophy_file,"('enstrophy.',i2.2)") ifr_run+1
             write(diag2_time_file,"('diag2_time.',i2.2)") ifr_run+1
             write(kes_file,"('kes.',i2.2)") ifr_run+1
             ifr_run=ifr_run+1  ! Update for next restart
             ! 26/3/2019
             ! For Lagrangian SQG+1
             if (lagrangian) then
                if (lagr_start_frame==ci) then
                   call Message('Error: lagr_start_frame must be set when restarting Lagrangian dynamics')
                   Check_parameters=.false.
                else
                   lagr_frame=0
                   xtra_restart_file = xtra_file
                   ytra_restart_file = ytra_file
                   write(xtra_file,"('xtra.',i2.2)") lagr_ifr_run+1
                   write(ytra_file,"('ytra.',i2.2)") lagr_ifr_run+1
                   if (comp_lag_grad) then
                      if (.not.advect_g) then 
                         write(Div_tra_file,"('Div_tra.',i2.2)") lagr_ifr_run+1
                         write(Vort_a_tra_file,"('Vorta_tra.',i2.2)") lagr_ifr_run+1
                         write(Vort_phi_tra_file,"('Vortphi_tra.',i2.2)") lagr_ifr_run+1
                      endif
                      write(Vort_g_tra_file,"('Vortg_tra.',i2.2)") lagr_ifr_run+1
                   endif
                   lagr_ifr_run=lagr_ifr_run+1  ! Update for next restart
                endif
             endif
          endif
       endif
       if (psi_restart_file=='') then
          call Message('Will read streamfunction from: '//trim(psi_file)//&
           &', frame:',tag=start_frame+1)
       else
          call Message('Will read streamfunction from: '//trim(psi_restart_file)//&
           &', frame:',tag=start_frame+1)
       endif
       ! ******************************       
       ! 29/3/2021 
       if (uphi_restart_file=='') then
          call Message('Will read u_phi from: '//trim(uphi_file)//&
           &', frame:',tag=start_frame+1)
       else
          call Message('Will read u_phi from: '//trim(uphi_restart_file)//&
           &', frame:',tag=start_frame+1)
       endif
       if (uag_restart_file=='') then
          call Message('Will read u_a from: '//trim(uag_file)//&
           &', frame:',tag=start_frame+1)
       else
          call Message('Will read u_a from: '//trim(uag_restart_file)//&
           &', frame:',tag=start_frame+1)
       endif
       ! ******************************       
       ! 26/3/2019
       ! For Lagrangian SQG+1
       if (lagrangian) then
          if (xtra_restart_file=='') then
             call Message('Will read Lagrangian coordinates x from: '//trim(xtra_file)//&
               &', frame:',tag=lagr_start_frame+1)
          else
             call Message('Will read Lagrangian coordinates x from: '//trim(xtra_restart_file)//&
               &', frame:',tag=lagr_start_frame+1)
          endif
          if (ytra_restart_file=='') then
             call Message('Will read Lagrangian coordinates y from: '//trim(ytra_file)//&
               &', frame:',tag=lagr_start_frame+1)
          else
             call Message('Will read Lagrangian coordinates y from: '//trim(ytra_restart_file)//&
               &', frame:',tag=lagr_start_frame+1)
          endif
       endif
       ! ******************************       
    elseif (.not.restarting) then

       frame = 0; d1frame = 0; d2frame = 0; cntr = 1

       call Message('Checking counters...')  ! Check internal counters

       if (total_counts==ci) then
          call Message('Error: total_counts not set')
          Check_parameters=.false.
       else
          call Message('Total number of timesteps model will run =', &
                        tag=total_counts)
       endif
       if (diag1_step==ci) then
          call Message('Info: timeseries write interval not set - &
               &setting diag1_step = 50')
          diag1_step = 50
       endif
       call Message('Timeseries will be written at timestep interval =', &
                     tag=diag1_step)
       if (diag2_step==ci.and.do_spectra) then
          call Message('Info: spectra write interval not set - &
               &setting diag2_step = 100')
          diag2_step = 100
       endif
       call Message('Spectra will be written at timestep interval =', &
                     tag=diag2_step)
       if (write_step==ci) then
          call Message('Error: Full field snapshot write interval not set - &
               &choose a value for write_step')
          Check_parameters=.false.
       else
          call Message('Snapshots will be written at timestep interval =', &
                     tag=write_step)
       endif
       if (write_lagr==ci) then
          call Message('Error: Lagrangian trajectories snapshot write interval not set - &
               &choose a value for write_lagr')
          Check_parameters=.false.
       else
          call Message('Lagrangian trajectories snapshots will be written at timestep interval =', &
                     tag=write_lagr)
       endif

    endif restart

! *******************************************************
    ! 26/5/2021 (Markovian forcing)
    if (use_forcing) then
       call Message('Random Markovian forcing on')
       if (forc_coef==cr) then
          call Message('Error: forc_coef must be set for RM forcing')
          Check_parameters=.false.
       else
          call Message('Forcing with coefficient forc_coef =',r_tag=forc_coef)
       endif
       if (forc_corr==cr) then
          call Message('Info: forc_corr not set - setting to .5')
          forc_corr = .5
       endif
       if ((forc_corr<0).or.(forc_corr>1)) then
          call Message('Error: require 0<=forc_corr<=1 - &
               &yours is:',&
               r_tag=forc_corr)
          Check_parameters=.false.
       endif
       if (kf_min==cr)  then
          call Message('Error: kf_min must be set for RM forcing')
          Check_parameters=.false.
       else
          call Message('Minimum forcing wavenumber kf_min =',tag=int(kf_min))
       endif
       if (kf_max==cr)  then
          call Message('Error: kf_max must be set for RM forcing')
          Check_parameters=.false.
       else
          call Message('Maximum forcing wavenumber kf_max =',tag=int(kf_max))
       endif
       if (norm_forcing) call Message('Info: eddy generation &
            &rate set to forc_coef')
    endif
! *******************************************************

  end function Check_parameters

!*************************************************************************

end module qg_params
