module qg_run_tools

  ! Contains core tools for running the model:
  ! PV inversion, time-stepping and forcing.
  !
  ! Routines: Get_pv, Invert_pv, March, Markovian, Write_snapshots
  !

  !
  ! Dependencies: io_tools, qg_arrays, qg_params, qg_strat_tools,
  !               transform_tools, numerics_lib

  implicit none
  save

contains


  !*******************************************************************

  function Get_pv(psi) result(q)

    !**************************************************************
    ! Calculate PV field from streamfunction
    !**************************************************************

    use qg_params, only: kmax
    use qg_arrays, only: ksqd_

    complex,dimension(-kmax:kmax,0:kmax),intent(in) :: psi
    complex,dimension(-kmax:kmax,0:kmax)            :: q


    q = -sqrt(ksqd_)*psi

  end function Get_pv

  !*******************************************************************

  function Invert_pv(q) result(psi)

    !**************************************************************
    ! Invert PV (q) in manner depending on surface_bc, which can be
    ! either 'rigid_lid', 'free_surface', or 'periodic'.  These
    ! types are checked in qg_params/Check_parameters right now.
    !**************************************************************

    use qg_params,    only: kmax
    use qg_arrays,    only: ksqd_
    complex,dimension(-kmax:kmax,0:kmax) :: psi, q



    psi = - (1./sqrt(ksqd_))*q

    psi(-kmax:0,0) = 0.

  end function Invert_pv



  !*********************************************************************

    subroutine Get_rhs
      use qg_arrays,       only: rhs,thetag,kx_,ky_,psi,filter,q,&
                                 u_a,u_phi,div, uvd_total,filter,ksqd_,ktdrag_,force_o,&
                                 u,v                           ! 29/3/2022
      use qg_params,       only: kmax,nx,ny,linear,i,rossby,&
                                  therm_drag,&
                                  kf_max,kf_min,&
                                  mean_divtemp,dt,&
                                  advect_g                     ! 29/3/2022
      use transform_tools, only: Spec2grid_cc, Grid2spec, ir_prod
      implicit none
      complex,dimension(-kmax:kmax,0:kmax,3) :: toto3
      complex,dimension(nx,ny,3)             :: totog3
      complex,dimension(-kmax:kmax,0:kmax,2) :: uv

      rhs=0.
      mean_divtemp=0.

      thetag=0.
      toto3=0.
      totog3=0.
      uvd_total=0.
      uv=0.

      if (.not.linear) then

        thetag = Spec2grid_cc(q)

        ! vitesses geostrophiques
        uv(:,:,1) = -i*(ky_*psi)
        uv(:,:,2) =  i*(kx_*psi)

        if (advect_g) then 
           u=Spec2grid_cc(uv(:,:,1))
           v=Spec2grid_cc(uv(:,:,2))
        endif

        if(rossby>0.)then

          call get_ageostrophic

          div=0.

          ! calcul des vitesses et de la divergence ageostrophique
          !    u_1 ne contribue pas ???
          uv(:,:,1) = uv(:,:,1) &
              + rossby*(u_phi(:,:,1) + u_a(:,:,1))

          uv(:,:,2)= uv(:,:,2) &
              + rossby*(u_phi(:,:,2) + u_a(:,:,2))

          ! calcul divergence
          div=rossby*(i*(kx_*u_a(:,:,1) + ky_*u_a(:,:,2)))
          uvd_total(:,:,3)=Spec2grid_cc(div)

          ! calcul div*theta
          totog3(:,:,3)=ir_prod(uvd_total(:,:,3), thetag)
          !  moyenne de div*theta
          mean_divtemp=real(sum(totog3(:,:,3)))/float(nx*ny)

          toto3(:,:,3)=Grid2spec(totog3(:,:,3))


        endif

        uvd_total(:,:,1)=Spec2grid_cc(uv(:,:,1))
        uvd_total(:,:,2)=Spec2grid_cc(uv(:,:,2))

        if (.not.advect_g) then 
           u=uvd_total(:,:,1)
           v=uvd_total(:,:,2)
        endif

        ! flux total de temperature
        totog3(:,:,1)=ir_prod(uvd_total(:,:,1), thetag)
        totog3(:,:,2)=ir_prod(uvd_total(:,:,2), thetag)


        toto3(:,:,1)=Grid2spec(totog3(:,:,1))
        toto3(:,:,2)=Grid2spec(totog3(:,:,2))

        ! divergence du flux
        rhs(:,:)= -i*(kx_*toto3(:,:,1) + ky_*toto3(:,:,2))

        if(rossby>0.)rhs = rhs + toto3(:,:,3)


      endif

     ! ******************************************
     ! 26/1/2022
     ! hypofriction
     if (therm_drag/=0) rhs = rhs + therm_drag*(ktdrag_*psi)
     ! ******************************************

    end subroutine Get_rhs

    !*********************************************************************

    subroutine Get_ageostrophic
      use qg_arrays,       only: thetag,kx_,ky_,psi,filter, &
                                u_a,u_phi,ksqd_,q
      use qg_params,       only: kmax,nx,ny,i
      use transform_tools, only: Spec2grid_cc, Grid2spec, ir_prod
      implicit none
      complex,dimension(-kmax:kmax,0:kmax) :: tmp
      complex,dimension(-kmax:kmax,0:kmax,2) :: theta2z,uvtz,uvzt,toto2
      complex,dimension(nx,ny,2)             :: uvg,uvzg
      complex,dimension(nx,ny)               :: dthetagdz,tmpg

      u_phi=0.
      u_a=0.

      tmp=0.
      dthetagdz=0.
      tmpg=0.
      theta2z=0.
      uvg=0.
      uvzg=0.
      uvtz=0.
      uvzt=0.
      toto2=0.


      !    dthetadz
      tmp = -sqrt(ksqd_)*q
      dthetagdz=spec2grid_cc(tmp)

      ! calcul de psi_1_tilde= - F[ - theta * d_theta_dz] /K
      ! theta^2/2 ne contribue pas au flux de theta, donc on pourrait l'oublier
      tmpg=ir_prod(thetag,thetag)
      theta2z(:,:,1)=grid2spec(tmpg)

      tmpg=ir_prod(thetag,dthetagdz)
      theta2z(:,:,2)=grid2spec(tmpg)

      tmp = filter*(0.5*theta2z(:,:,1)+(1./sqrt(ksqd_))*theta2z(:,:,2))

      u_phi(:,:,1)=-i*(ky_*tmp)
      u_phi(:,:,2)= i*(kx_*tmp)


      ! calcul de (ug, vg)  et (d_u_dz d_v_dz)
      tmp = -i*(ky_*psi)
      uvg(:,:,1)=Spec2grid_cc(tmp)
      tmp = i*(kx_*psi)
      uvg(:,:,2)=Spec2grid_cc(tmp)

      tmp = i*((ky_*sqrt(ksqd_))*psi)
      uvzg(:,:,1)=Spec2grid_cc(tmp)
      tmp = -i*((kx_*sqrt(ksqd_))*psi)
      uvzg(:,:,2)=Spec2grid_cc(tmp)

      !  (u d_theta_dz,   d_u_dz theta)
      tmpg=ir_prod(uvg(:,:,1),dthetagdz)
      uvtz(:,:,1)=Grid2spec(tmpg)

      tmpg=ir_prod(uvg(:,:,2),dthetagdz)
      uvtz(:,:,2)=Grid2spec(tmpg)

      tmpg=ir_prod(uvzg(:,:,1),thetag)
      uvzt(:,:,1)=Grid2spec(tmpg)
      tmpg=ir_prod(uvzg(:,:,2),thetag)
      uvzt(:,:,2)=Grid2spec(tmpg)


      !  (ug theta, vg theta)
      tmpg=ir_prod(uvg(:,:,1),thetag)
      toto2(:,:,1)=Grid2spec(tmpg)

      tmpg=ir_prod(uvg(:,:,2),thetag)
      toto2(:,:,2)=Grid2spec(tmpg)

      ! en theorie, toto4(:,:,3:4) ne contribue pas aux flux
      u_a(:,:,1)=filter*( uvtz(:,:,1)+uvzt(:,:,1) &
             +sqrt(ksqd_)*toto2(:,:,1) )
      u_a(:,:,2)=filter*( uvzt(:,:,2)+uvtz(:,:,2) &
             +sqrt(ksqd_)*toto2(:,:,2) )


    end subroutine Get_ageostrophic


    !*********************************************************************

  function Markovian(kf_min,kf_max,amp,lambda,frc_o,norm_forcing,field) &
       result(forc)

    !**************************************************************
    ! Random Markovian forcing function.  If norm_forcing = T, function
    ! will normalize the forcing such that the total generation = amp
    !**************************************************************

    use qg_arrays,    only: ksqd_
    use qg_params,    only: nkx, nky, kmax, idum, i, pi, rmf_norm_min
    use numerics_lib, only: Ran
    use io_tools,     only: Message

    complex,dimension(-kmax:kmax,0:kmax)         :: forc
    complex,dimension(-kmax:kmax,0:kmax),intent(inout)  :: frc_o
    real,intent(in)                                 :: kf_min,kf_max,amp,lambda
    logical,intent(in),optional                     :: norm_forcing
    complex,dimension(:,:),intent(in),optional    :: field
    ! Local
    real                                            :: gamma=1.

    if (present(norm_forcing)) then
       if ((norm_forcing).and.(.not.(present(field)))) then
          call Message('Error:Markovian: need field with norm_forcing',&
                        fatal='y')
       endif
    endif

    where((ksqd_ > kf_min**2).and.(ksqd_ <= kf_max**2))
       frc_o = lambda*frc_o &
             + amp*sqrt(1-lambda**2)*cexp(i*2*pi*Ran(idum,nkx,nky))
    endwhere
    forc(:,:) = frc_o  ! Force top layer only
    if (norm_forcing) then
       gamma = real(-2*sum((conjg(field)*forc)))/amp
       if (abs(gamma)>rmf_norm_min) then
           forc = forc/gamma
       else
          call Message('RandMarkForc: normalization factor too small, =',&
                       r_tag = gamma)
       endif
    endif
    forc(-kmax:0,0) = 0.
    frc_o = forc(:,:)

  end function Markovian

  !*********************************************************************

  function Write_snapshots(framein) result(frameout)

    !**************************************************************
    ! Write full snapshot of dynamic field and make restart files
    !**************************************************************

    use io_tools,    only: Message, Write_field
    use qg_params,   only: psi_file,time,write_time_file,&
                           energy_file,enstrophy_file,&
                           rossby,Write_parameters,dt, &
                           uphi_file,uag_file
    use qg_diagnostics, only : energy, enstrophy
    use qg_arrays,   only: psi,q,u_phi,u_a

    integer,intent(in) :: framein
    integer            :: frameout
    real               :: scra

    frameout = framein + 1                ! Update field frame counter

    call Write_field(psi,psi_file,frameout)
    ! **************************************
    ! 8/12/2021
    ! kinetic energy and enstrophy 
    ! of the geostrophic flow
    scra=0.
    scra=energy(psi)
    call message('kinetic energy (geostrophic): ',r_tag=scra)
    call Write_field(scra,energy_file,frameout)
    scra=0.
    scra=enstrophy(q)
    call message('enstrophy (geostrophic): ',r_tag=scra)
    call Write_field(scra,enstrophy_file,frameout)
    ! **************************************

    if(rossby>0.)then
        call Write_field(u_phi,uphi_file,frameout)
        call Write_field(u_a,uag_file,frameout)
    endif

    call Write_parameters                 ! Write params for restart.nml
    call Write_field(time,write_time_file,frameout)

    call Message('Wrote frame: ',tag=frameout)

  end function Write_snapshots

  !*********************************************************************

  !*********************************************************************  
! 24/3/2022: Lagrangian tracers dynamics (Adams-Bashforth order 3)
  subroutine dynlaab

    use qg_params, only: nx,ny,pi,dt,ntra,cab,cab_1,cab_0
    use qg_arrays, only: u,v,u_1,v_1,u_0,v_0,xtra,ytra,xtra_1,ytra_1,xtra_0,ytra_0

    real,dimension(ntra,3) :: up,vp
    real,dimension(ntra) :: xt,yt
    real,dimension(ntra) :: uu,vv
    integer,dimension(ntra) :: ix,jy
    real,dimension(ntra,4) :: eg,ef
    integer,dimension(ntra,4) :: ic,jc
    real,dimension(ntra,4,4) :: au,av
    real :: dx, r
    integer :: l,i,k,ii,jj,j


    dx=2.*pi/float(nx)

!----------------------------------------------------
!  adams-bashforth integration
!----------------------------------------------------
! x(n+3)=x(n+2)+dt*(cab*u-cab_1*u_1+cab_0*u_0)
! u(x(n+2),t(n+2))
! u_1(x(n+1),t(n+1))
! u_0(x(n),t(n))
!----------------------------------------------------

       up=0.
       vp=0.

!----------------------------------------------------
       do  l=1,3
!----------------------------------------------------
!   positions for interpolation 
!   and their normalization on the grid 

         if (l==1) then
           xt(:)=xtra(:)/dx+1.
           yt(:)=ytra(:)/dx+1.
         endif
         if (l==2) then
           xt(:)=xtra_1(:)/dx+1.
           yt(:)=ytra_1(:)/dx+1.
         endif
         if (l==3) then
           xt(:)=xtra_0(:)/dx+1.
           yt(:)=ytra_0(:)/dx+1.
         endif

!-----------------------------------------------------
!   interpolation of order 3  of velocity in (xt,yt)
!-----------------------------------------------------
! ix(k) et jy(k) sont les coordonnees entieres du coin inferieur 
! gauche de la maille qui contient le point a interpoler

          ix(:)=int(xt(:)+(sign(1.,xt(:))-1.)/2.)
          jy(:)=int(yt(:)+(sign(1.,yt(:))-1.)/2.)

!---------------------------------------------------
! calcul des coordonnees des 16 points de grille entourant
!     le point ou l'on interpole la vitesse
! ces coordonnees sont ramenees entre 1 et nx
!---------------------------------------------------

          do  i=1,4
             do  k=1,ntra
                ii=ix(k)+i-2
                jj=jy(k)+i-2
                if (ii.gt.0) then
                   ic(k,i)=1+mod(ii-1,nx)
                else
                   ic(k,i)=mod(ii,nx)+nx
                endif
                if (jj.gt.0) then
                   jc(k,i)=1+mod(jj-1,ny)
                else
                   jc(k,i)=mod(jj,ny)+ny
                endif
             enddo
          enddo

!------------------------------------------------
! calcul des coefficients ef(k,i) et eg(k,j)
!  de la formule d'interpolation.
!------------------------------------------------
          ef(:,1)=-(xt(:)-float(ix(:)))*(xt(:)-float(ix(:))-1.)*(xt(:)-float(ix(:))-2.)/6.
          ef(:,2)=(xt(:)-float(ix(:))+1.)*(xt(:)-float(ix(:))-1.)*(xt(:)-float(ix(:))-2.)/2.
          ef(:,3)=-(xt(:)-float(ix(:))+1.)*(xt(:)-float(ix(:)))*(xt(:)-float(ix(:))-2.)/2.
          ef(:,4)=(xt(:)-float(ix(:))+1.)*(xt(:)-float(ix(:)))*(xt(:)-float(ix(:))-1.)/6.
          eg(:,1)=-(yt(:)-float(jy(:)))*(yt(:)-float(jy(:))-1.)*(yt(:)-float(jy(:))-2.)/6.
          eg(:,2)=(yt(:)-float(jy(:))+1.)*(yt(:)-float(jy(:))-1.)*(yt(:)-float(jy(:))-2.)/2.
          eg(:,3)=-(yt(:)-float(jy(:))+1.)*(yt(:)-float(jy(:)))*(yt(:)-float(jy(:))-2.)/2.
          eg(:,4)=(yt(:)-float(jy(:))+1.)*(yt(:)-float(jy(:)))*(yt(:)-float(jy(:))-1.)/6.


! on prend les vitesses a t
          if (l.eq.1) then

          do  i=1,4
             do  j=1,4
                do  k=1,ntra
                      au(k,i,j)=u(ic(k,i),jc(k,j))
                      av(k,i,j)=v(ic(k,i),jc(k,j))
                   enddo
                enddo
             enddo

          endif

! calcul du champ des vitesses a t-dt
          if (l.eq.2) then

          do  i=1,4
             do  j=1,4
                do  k=1,ntra
                      au(k,i,j)=u_1(ic(k,i),jc(k,j))
                      av(k,i,j)=v_1(ic(k,i),jc(k,j))
                   enddo
                enddo
             enddo

          endif

! calcul du champ des vitesses a t-2dt
          if (l.eq.3) then

          do  i=1,4
             do  j=1,4
                do  k=1,ntra
                      au(k,i,j)=u_0(ic(k,i),jc(k,j))
                      av(k,i,j)=v_0(ic(k,i),jc(k,j))
                   enddo
                enddo
             enddo

          endif

!----------------------------------------------------
          do j=1,4
             do i=1,4
               do k=1,ntra
                  r=ef(k,i)*eg(k,j)
                  up(k,l)=r*au(k,i,j)+up(k,l)
                  vp(k,l)=r*av(k,i,j)+vp(k,l)
               enddo
             enddo
          enddo

!----------------------------------------------------
!---------- end of interpolation ------------------
!----------------------------------------------------


!----------------------------------------------------
       enddo
!--------------------------------------------------

!      displacement 
!      dx=cab*up(1)-cab_1*up(2)+cab_0*up(3)
!      dy=cab*vp(1)-cab_1*vp(2)+cab_0*vp(3)
!--------------------------------------------------

       uu=cab*up(:,1)-cab_1*up(:,2)+cab_0*up(:,3)
       vv=cab*vp(:,1)-cab_1*vp(:,2)+cab_0*vp(:,3)

!--------------------------------------------------------
!   end of adams-bashforth integration 
!--------------------------------------------------------

!-----------------------------------------------------
!   advance Lagrangian tracers
!-----------------------------------------------------

       xtra(:)=xtra(:)+dt*uu(:)
       ytra(:)=ytra(:)+dt*vv(:)

!--------------------------------------------------------

  end subroutine dynlaab

  !*********************************************************************  
  ! 29/3/2022
  ! Velocity gradients (div, vort) at Lagrangian positions
  subroutine lag_grad(ppsi,uu_phi,uu_a)

  use qg_params, only: kmax,nx,ny,ntra,pi,rossby,lagr_frame, &
                       Div_tra_file,Vort_phi_tra_file,Vort_a_tra_file, &
                       Vort_g_tra_file,advect_g
  use qg_arrays, only: xtra,ytra,kx_,ky_
  use transform_tools, only: Spec2grid_cc
  use io_tools, only: Write_field

  complex,parameter :: ij=(0.,1.)
  complex,dimension(-kmax:kmax,0:kmax) :: ppsi,scra
  complex,dimension(-kmax:kmax,0:kmax,2) :: uu_phi,uu_a
  real,dimension(nx,ny) :: Div,Vort_phi,Vort_a,Vort_g
  real,dimension(ntra) :: xt,yt,Div_tra,Vort_phi_tra,Vort_a_tra,Vort_g_tra 
  integer,dimension(ntra) :: ix,jy
  real,dimension(ntra,4) :: eg,ef
  integer,dimension(ntra,4) :: ic,jc
  real :: dx
  integer :: i,j,k,ii,jj

  if (.not.advect_g) then 
     ! compute divergence (div(u_phi)=0)
     scra(:,:)=ij*kx_*(uu_a(:,:,1)+uu_phi(:,:,1))+ij*ky_*(uu_a(:,:,2)+uu_phi(:,:,2))
     Div=rossby*Spec2grid_cc(scra)
     Div_tra=0.
     ! compute vorticity
     scra(:,:)=ij*kx_*uu_phi(:,:,2)-ij*ky_*(uu_phi(:,:,1))
     Vort_phi=rossby*Spec2grid_cc(scra)
     scra(:,:)=ij*kx_*uu_a(:,:,2)-ij*ky_*(uu_a(:,:,1))
     Vort_a=rossby*Spec2grid_cc(scra)
     Vort_phi_tra=0.
     Vort_a_tra=0.
  endif
  ! geostrophic vorticity
  scra=-(kx_*kx_+ky_*ky_)*ppsi
  Vort_g=Spec2grid_cc(scra)
  Vort_g_tra=0.

  ! interpolation
  dx=2.*pi/float(nx)

!----------------------------------------------------
! on normalise les coordonnees des flotteurs
! (elles peuvent etre negatives)
!----------------------------------------------------
   xt(:)=xtra(:)/dx+1.
   yt(:)=ytra(:)/dx+1.

!-----------------------------------------------------
!   interpolation d'ordre 3  des champs en (xt,yt)
!-----------------------------------------------------
! ix(k) et jy(k) sont les coordonnees entieres du coin inferieur 
! gauche de la maille qui contient le point a interpoler

   ix(:)=int(xt(:)+(sign(1.,xt(:))-1.)/2.)
   jy(:)=int(yt(:)+(sign(1.,yt(:))-1.)/2.)

!---------------------------------------------------
! calcul des coordonnees des 16 points de grille entourant
! le point ou l'on interpole le champ
! ces coordonnees sont ramenees entre 1 et nx
!---------------------------------------------------

   do  i=1,4
      do  k=1,ntra
         ii=ix(k)+i-2
         jj=jy(k)+i-2
         if (ii.gt.0) then
            ic(k,i)=1+mod(ii-1,nx)
         else
            ic(k,i)=mod(ii,nx)+nx
         endif
         if (jj.gt.0) then
            jc(k,i)=1+mod(jj-1,ny)
         else
            jc(k,i)=mod(jj,ny)+ny
         endif
      enddo
   enddo

!------------------------------------------------
! calcul des coefficients ef(k,i) et eg(k,j)
!  de la formule d'interpolation.
!------------------------------------------------
   ef(:,1)=-(xt(:)-float(ix(:)))*(xt(:)-float(ix(:))-1.)*(xt(:)-float(ix(:))-2.)/6.
   ef(:,2)=(xt(:)-float(ix(:))+1.)*(xt(:)-float(ix(:))-1.)*(xt(:)-float(ix(:))-2.)/2.
   ef(:,3)=-(xt(:)-float(ix(:))+1.)*(xt(:)-float(ix(:)))*(xt(:)-float(ix(:))-2.)/2.
   ef(:,4)=(xt(:)-float(ix(:))+1.)*(xt(:)-float(ix(:)))*(xt(:)-float(ix(:))-1.)/6.
   eg(:,1)=-(yt(:)-float(jy(:)))*(yt(:)-float(jy(:))-1.)*(yt(:)-float(jy(:))-2.)/6.
   eg(:,2)=(yt(:)-float(jy(:))+1.)*(yt(:)-float(jy(:))-1.)*(yt(:)-float(jy(:))-2.)/2.
   eg(:,3)=-(yt(:)-float(jy(:))+1.)*(yt(:)-float(jy(:)))*(yt(:)-float(jy(:))-2.)/2.
   eg(:,4)=(yt(:)-float(jy(:))+1.)*(yt(:)-float(jy(:)))*(yt(:)-float(jy(:))-1.)/6.

!----------------------------------------------------
  if (.not.advect_g) then
     do j=1,4
        do i=1,4
           do k=1,ntra
              Div_tra(k)=ef(k,i)*eg(k,j)*Div(ic(k,i),jc(k,j))+Div_tra(k)
              Vort_phi_tra(k)=ef(k,i)*eg(k,j)*Vort_phi(ic(k,i),jc(k,j))+Vort_phi_tra(k)
              Vort_a_tra(k)=ef(k,i)*eg(k,j)*Vort_a(ic(k,i),jc(k,j))+Vort_a_tra(k)
           enddo
        enddo
     enddo
     call Write_field(Div_tra,Div_tra_file,lagr_frame) 
     call Write_field(Vort_phi_tra,Vort_phi_tra_file,lagr_frame) 
     call Write_field(Vort_a_tra,Vort_a_tra_file,lagr_frame) 
  endif  
  ! geostrophic vortcity
  do j=1,4
     do i=1,4
        do k=1,ntra
           Vort_g_tra(k)=ef(k,i)*eg(k,j)*Vort_g(ic(k,i),jc(k,j))+Vort_g_tra(k)
        enddo
     enddo
  enddo
  call Write_field(Vort_g_tra,Vort_g_tra_file,lagr_frame) 


!----------------------------------------------------
!---------- fin de l'interpolation ------------------
!----------------------------------------------------

  end subroutine lag_grad

  !*************************************************************************

end module qg_run_tools
