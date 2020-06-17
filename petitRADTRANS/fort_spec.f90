!!$*****************************************************************************
!!$*****************************************************************************
!!$*****************************************************************************
!!$ fort_spec.f90: utility functions to calculate cloud opacities, optical
!!$ depths, spectra, and spectral contribution functions for the petitRADTRANS
!!$ radiative transfer package.
!!$
!!$ Copyright 2016-2018, Paul Molliere
!!$ Maintained by Paul Molliere, molliere@strw.leidenunivl.nl
!!$ Status: under development
!!$*****************************************************************************
!!$*****************************************************************************
!!$*****************************************************************************

!!$ Natural constants block

module constants_block
  implicit none
  DOUBLE PRECISION,parameter      :: AU = 1.49597871d13, R_sun = 6.955d10, R_jup=6.9911d9
  DOUBLE PRECISION,parameter      :: pi = 3.14159265359d0, sig=5.670372622593201d-5, c_l=2.99792458d10
  DOUBLE PRECISION,parameter      :: G = 6.674d-8, M_jup = 1.89813e30, deg = Pi/1.8d2
  DOUBLE PRECISION,parameter      :: kB=1.3806488d-16, hplanck=6.62606957d-27, amu = 1.66053892d-24
  DOUBLE PRECISION,parameter      :: sneep_ubachs_n = 25.47d18, L0 = 2.68676d19
end module constants_block

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate tau with 2nd order accuracy

subroutine calc_tau_g_tot_ck(gravity,press,total_kappa,struc_len,freq_len,g_len,N_species,tau)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                          :: struc_len, freq_len, g_len, N_species
  DOUBLE PRECISION, intent(in)                 :: total_kappa(g_len,freq_len,N_species,struc_len)
  DOUBLE PRECISION, intent(in)                 :: gravity, press(struc_len)
  DOUBLE PRECISION, intent(out)                :: tau(g_len,freq_len,N_species,struc_len)
  ! internal
  integer                                      :: i_struc, i_freq, i_g, i_spec
  DOUBLE PRECISION                             :: del_tau_lower_ord, &
       gamma_second(g_len,freq_len,N_species), f_second, kappa_i(g_len,freq_len,N_species), &
       kappa_im(g_len,freq_len,N_species), kappa_ip(g_len,freq_len,N_species)
  LOGICAL                                      :: second_order
  !~~~~~~~~~~~~~

  tau = 0d0
  second_order = .FALSE.

  if (second_order) then
     do i_struc = 2, struc_len
        if (i_struc .EQ. struc_len) then
           tau(:,:,:,i_struc) = tau(:,:,:,i_struc-1) + &
                (total_kappa(:,:,:,i_struc)+total_kappa(:,:,:,i_struc-1)) &
                /2d0/gravity*(press(i_struc)-press(i_struc-1))
        else
           f_second = (press(i_struc+1)-press(i_struc))/(press(i_struc)-press(i_struc-1))
           kappa_i = total_kappa(:,:,:,i_struc)
           kappa_im = total_kappa(:,:,:,i_struc-1)
           kappa_ip = total_kappa(:,:,:,i_struc+1)
           gamma_second = (kappa_ip-(1d0+f_second)*kappa_i+f_second*kappa_im) / &
                (f_second*(1d0+f_second))
           tau(:,:,:,i_struc) = tau(:,:,:,i_struc-1) + &
                ((kappa_i+kappa_im)/2d0-gamma_second/6d0) &
                /gravity*(press(i_struc)-press(i_struc-1))
           do i_spec = 1, N_species
              do i_freq = 1, freq_len
                 do i_g = 1, g_len
                    if (tau(i_g,i_freq,i_spec,i_struc) < tau(i_g,i_freq,i_spec,i_struc-1)) then
                       if (i_struc .EQ. 2) then
                          tau(i_g,i_freq,i_spec,i_struc) = &
                               tau(i_g,i_freq,i_spec,i_struc-1)*1.01d0
                       else
                          tau(i_g,i_freq,i_spec,i_struc) = &
                               tau(i_g,i_freq,i_spec,i_struc-1) + &
                               (tau(i_g,i_freq,i_spec,i_struc-1)- &
                               tau(i_g,i_freq,i_spec,i_struc-2))*0.01d0
                       end if
                    end if
                    del_tau_lower_ord = (kappa_i(i_g,i_freq,i_spec)+ &
                         kappa_im(i_g,i_freq,i_spec))/2d0/gravity* &
                         (press(i_struc)-press(i_struc-1))
                    if ((tau(i_g,i_freq,i_spec,i_struc) - &
                         tau(i_g,i_freq,i_spec,i_struc-1)) > del_tau_lower_ord) then
                       tau(i_g,i_freq,i_spec,i_struc) = &
                            tau(i_g,i_freq,i_spec,i_struc-1) + del_tau_lower_ord
                    end if
                 end do
              end do
           end do
        end if
     end do
  else
     do i_struc = 2, struc_len
        tau(:,:,:,i_struc) = tau(:,:,:,i_struc-1) + &
             (total_kappa(:,:,:,i_struc)+total_kappa(:,:,:,i_struc-1)) &
             /2d0/gravity*(press(i_struc)-press(i_struc-1))
     end do
  end if

end subroutine calc_tau_g_tot_ck

!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to do the radiative transport, using the mean transmission method

subroutine flux_ck(freq,tau,temp,mu,w_gauss_mu, &
     w_gauss,contribution,freq_len,struc_len,N_mu,g_len,N_species,flux,contr_em)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                         :: freq_len, struc_len,g_len, N_species
  DOUBLE PRECISION, intent(in)                :: freq(freq_len)
  DOUBLE PRECISION, intent(in)                :: temp(struc_len) !, press(struc_len)
  DOUBLE PRECISION, intent(in)                :: tau(g_len,freq_len,N_species,struc_len)

  INTEGER, intent(in)                         :: N_mu
  DOUBLE PRECISION, intent(in)                :: mu(N_mu) !, gravity
  DOUBLE PRECISION, intent(in)                :: w_gauss_mu(N_mu)
  DOUBLE PRECISION, intent(in)                :: w_gauss(g_len)
  LOGICAL, intent(in)                         :: contribution
  DOUBLE PRECISION, intent(out)               :: flux(freq_len)
  DOUBLE PRECISION, intent(out)               :: contr_em(struc_len,freq_len)

  ! Internal
  INTEGER                                     :: i_mu,i_freq,i_str,i_spec
  DOUBLE PRECISION                            :: r(struc_len)
  DOUBLE PRECISION                            :: transm_mu(g_len,freq_len,N_species,struc_len), &
       mean_transm(freq_len,N_species,struc_len), transm_all(freq_len,struc_len), &
       transm_all_loc(struc_len), flux_mu(freq_len)

  flux = 0d0

  if (contribution) then
     contr_em = 0d0
  end if

  do i_mu = 1, N_mu

     ! will contain species' product of g-space integrated transmissions
     transm_all = 1d0
     ! Transmissions for a given incidence angle
     transm_mu = exp(-tau/mu(i_mu))
     ! Flux contribution from given mu-angle
     flux_mu = 0d0

     do i_str = 1, struc_len
        do i_spec = 1, N_species
           do i_freq = 1, freq_len
              ! Integrate transmission over g-space
              mean_transm(i_freq,i_spec,i_str) = sum(transm_mu(:,i_freq,i_spec,i_str)*w_gauss)
           end do
        end do
     end do

     ! Multiply transmissions of infdiv. species
     do i_spec = 1, N_species
        transm_all = transm_all*mean_transm(:,i_spec,:)
     end do

     ! Do the actual radiative transport
     do i_freq = 1, freq_len
        ! Get source function
        r = 0
        call planck_f(struc_len,temp,freq(i_freq),r)
        ! Spatial transmissions at given wavelength
        transm_all_loc = transm_all(i_freq,:)
        ! Calc Eq. 9 of manuscript (em_deriv.pdf) 
        do i_str = 1, struc_len-1
           flux_mu(i_freq) = flux_mu(i_freq)+ &
                (r(i_str)+r(i_str+1))*(transm_all_loc(i_str)-transm_all_loc(i_str+1))/2d0
           if (contribution) then
              contr_em(i_str,i_freq) = contr_em(i_str,i_freq) + (r(i_str)+r(i_str+1)) * &
                   (transm_all_loc(i_str)-transm_all_loc(i_str+1)) &
                   *mu(i_mu)*w_gauss_mu(i_mu)              
           end if
        end do
        flux_mu(i_freq) = flux_mu(i_freq) + r(struc_len)*transm_all_loc(struc_len)
        if (contribution) then
           contr_em(struc_len,i_freq) = contr_em(struc_len,i_freq) + &
                2d0*r(struc_len)*transm_all_loc(struc_len)*mu(i_mu)*w_gauss_mu(i_mu)
        end if
     end do
     ! angle integ, factor 1/2 needed for flux calc. from upward pointing intensity
     flux = flux + flux_mu/2d0*mu(i_mu)*w_gauss_mu(i_mu)

  end do
  ! Normalization
  flux = flux*4d0*pi

  if (contribution) then
     do i_freq = 1, freq_len
        contr_em(:,i_freq) = contr_em(:,i_freq)/SUM(contr_em(:,i_freq))
     end do
  end if

end subroutine flux_ck

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate the Planck source function

subroutine planck_f(struc_len,T,nu,B_nu)

  use constants_block
  implicit none
  INTEGER                         :: struc_len
  DOUBLE PRECISION                :: T(struc_len),B_nu(struc_len), nu
  DOUBLE PRECISION                :: buffer

  !~~~~~~~~~~~~~

  B_nu = 0d0
  buffer = 2d0*hplanck*nu**3d0/c_l**2d0
  B_nu = buffer / (exp(hplanck*nu/kB/T)-1d0)
  
end subroutine planck_f

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate the transmission spectrum

subroutine calc_transm_spec(freq,total_kappa_in,temp,press,gravity,mmw,P0_bar,R_pl, &
     w_gauss,scat,continuum_opa_scat,var_grav,transm,freq_len,struc_len,g_len,N_species)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                         :: freq_len, struc_len, g_len, N_species
  DOUBLE PRECISION, intent(in)                :: freq(freq_len), P0_bar, R_pl
  DOUBLE PRECISION, intent(in)                :: temp(struc_len), press(struc_len), mmw(struc_len)
  DOUBLE PRECISION, intent(in)                :: total_kappa_in(g_len,freq_len,N_species,struc_len)

  DOUBLE PRECISION, intent(in)                :: gravity
  DOUBLE PRECISION, intent(in)                :: w_gauss(g_len), continuum_opa_scat(freq_len,struc_len)
  LOGICAL, intent(in)                         :: scat !, contribution
  DOUBLE PRECISION, intent(out)               :: transm(freq_len) !, contr_tr(struc_len,freq_len)

  ! Internal
  DOUBLE PRECISION                            :: P0_cgs, rho(struc_len), radius(struc_len), &
       radius_var(struc_len), total_kappa(g_len,freq_len,N_species,struc_len)
  INTEGER                                     :: i_str, i_freq, i_g, i_spec, j_str
  LOGICAL                                     :: var_grav
  DOUBLE PRECISION                            :: alpha_t2(g_len,freq_len,N_species,struc_len-1)
  DOUBLE PRECISION                            :: t_graze(g_len,freq_len,N_species,struc_len), s_1, s_2, &
       t_graze_wlen_int(struc_len,freq_len), &
       alpha_t2_scat(freq_len,struc_len-1), t_graze_scat(freq_len,struc_len)

  total_kappa = total_kappa_in
  ! Some cloud opas can be < 0 sometimes, apparently.
  do i_str = 1, struc_len
     do i_spec = 1, N_species
        do i_freq = 1, freq_len
           do i_g = 1, g_len
              if (total_kappa(i_g,i_freq,i_spec,i_str) < 0d0) then
                 total_kappa(i_g,i_freq,i_spec,i_str) = 0d0
              end if
           end do
        end do
     end do
  end do
        
  transm = 0d0
  t_graze = 0d0
  t_graze_scat = 0d0

  ! Convert reference pressure to cgs
  P0_cgs = P0_bar*1d6
  ! Calculate density
  rho = mmw*amu*press/kB/temp
  ! Calculate planetary radius (in cm), assuming hydrostatic equilibrium
  call calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,var_grav,radius)

  ! Calc. mean free paths across grazing distances
  do i_str = 1, struc_len-1
     alpha_t2(:,:,:,i_str) = (total_kappa(:,:,:,i_str)*rho(i_str)+total_kappa(:,:,:,i_str+1)*rho(i_str+1))
  end do

  if (scat) then
     do i_str = 1, struc_len-1
        alpha_t2_scat(:,i_str) = (continuum_opa_scat(:,i_str)*rho(i_str)+ &
             continuum_opa_scat(:,i_str+1)*rho(i_str+1))
     end do
  end if
  
  ! Cacuclate grazing rays optical depths
  do i_str = 2, struc_len
     s_1 = sqrt(radius(1)**2d0-radius(i_str)**2d0)
     do j_str = 1, i_str-1
        if (j_str > 1) then
           s_1 = s_2
        end if
        s_2 = sqrt(radius(j_str+1)**2d0-radius(i_str)**2d0)
        t_graze(:,:,:,i_str) = t_graze(:,:,:,i_str)+alpha_t2(:,:,:,j_str)*(s_1-s_2)
     end do
  end do
  if (scat) then
     do i_str = 2, struc_len
        s_1 = sqrt(radius(1)**2d0-radius(i_str)**2d0)
        do j_str = 1, i_str-1
           if (j_str > 1) then
              s_1 = s_2
           end if
           s_2 = sqrt(radius(j_str+1)**2d0-radius(i_str)**2d0)
           t_graze_scat(:,i_str) = t_graze_scat(:,i_str)+alpha_t2_scat(:,j_str)*(s_1-s_2)
        end do
     end do
  end if

  ! Calculate transmissions, update tau array to store these
  t_graze = exp(-t_graze)
  if (scat) then
     t_graze_scat = exp(-t_graze_scat)
  end if
  
  t_graze_wlen_int = 1d0
  ! Wlen (in g-space) integrate transmissions
  do i_str = 2, struc_len ! i_str=1 t_grazes are 1 anyways
     do i_spec = 1, N_species
        do i_freq = 1, freq_len
           t_graze_wlen_int(i_str,i_freq) = t_graze_wlen_int(i_str,i_freq)* &
                SUM(t_graze(:,i_freq,i_spec,i_str)*w_gauss)
           if (scat .AND. (i_spec .EQ. 1)) then
              t_graze_wlen_int(i_str,i_freq) = t_graze_wlen_int(i_str,i_freq)* &
                   t_graze_scat(i_freq,i_str)
           end if
        end do
     end do
  end do

  ! Get effective area fraction from transmission
  t_graze_wlen_int = 1d0-t_graze_wlen_int

  ! Caculate planets effectice area (leaving out pi, because we want the radius in the end)
  do i_freq = 1, freq_len
     do i_str = 2, struc_len
        transm(i_freq) = transm(i_freq)+(t_graze_wlen_int(i_str-1,i_freq)*radius(i_str-1)+ &
             t_graze_wlen_int(i_str,i_freq)*radius(i_str))*(radius(i_str-1)-radius(i_str))
     end do
  end do
  ! Get radius
  transm = sqrt(transm+radius(struc_len)**2d0)

!!$  if (contribution) then
!!$     contr_tr = t_graze_wlen_int
!!$  end if
!!$
!!$  call calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,.FALSE.,radius)
!!$  call calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,.TRUE.,radius_var)
!!$  open(unit=10,file='rad_test.dat')
!!$  do i_str = 1, struc_len
!!$     write(10,*) press(i_str)*1d-6, radius(i_str)/R_jup, radius_var(i_str)/R_jup, rho(i_str)
!!$  end do
!!$  close(10)
  
end subroutine calc_transm_spec

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate the radius from the pressure grid

subroutine calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,var_grav,radius)

  implicit none
  ! I/O
  INTEGER, intent(in)                         :: struc_len
  DOUBLE PRECISION, intent(in)                :: P0_cgs
  DOUBLE PRECISION, intent(in)                :: temp(struc_len), press(struc_len), mmw(struc_len), &
       rho(struc_len)
  DOUBLE PRECISION, intent(in)                :: gravity, R_pl
  LOGICAL, intent(in)                         :: var_grav
  DOUBLE PRECISION, intent(out)               :: radius(struc_len)

  ! Internal
  INTEGER                                     :: i_str
  DOUBLE PRECISION                            :: R0, inv_rho(struc_len), integ_parab

  inv_rho = 1d0/rho
  
  radius = 0d0

  if (var_grav) then

     !write(*,*) '####################### VARIABLE GRAVITY'
     !write(*,*) '####################### VARIABLE GRAVITY'
     !write(*,*) '####################### VARIABLE GRAVITY'
     !write(*,*) '####################### VARIABLE GRAVITY'
     
     ! Calculate radius with vertically varying gravity, set up such that at P=P0, i.e. R=R_pl
     ! the planet has the predefined scalar gravity value
     do i_str = struc_len-1, 1, -1
        if ((press(i_str+1) > P0_cgs) .AND. (press(i_str) <= P0_cgs)) then
           if (i_str <= struc_len-2) then
              R0 = radius(i_str+1) - integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                   inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),P0_cgs,press(i_str+1))/gravity &
                   /R_pl**2d0
           else
              R0 = radius(i_str+1)-(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                   (press(i_str+1)-P0_cgs)/R_pl**2d0
           end if
        end if
        if (i_str <= struc_len-2) then
           radius(i_str) = radius(i_str+1) - integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),press(i_str),press(i_str+1))/gravity &
                /R_pl**2d0
        else
           radius(i_str) = radius(i_str+1)-(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                (press(i_str+1)-press(i_str))/R_pl**2d0
        end if
     end do
     R0 = 1d0/R_pl-R0
     radius = radius + R0
     radius = 1d0/radius
     
  else

     !write(*,*) '####################### CONSTANT GRAVITY'
     !write(*,*) '####################### CONSTANT GRAVITY'
     !write(*,*) '####################### CONSTANT GRAVITY'
     !write(*,*) '####################### CONSTANT GRAVITY'

     
     ! Calculate radius with vertically constant gravity
     do i_str = struc_len-1, 1, -1
        if ((press(i_str+1) > P0_cgs) .AND. (press(i_str) <= P0_cgs)) then
           if (i_str <= struc_len-2) then
              R0 = radius(i_str+1) + integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                   inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),P0_cgs,press(i_str+1))/gravity
           else
              R0 = radius(i_str+1)+(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                   (press(i_str+1)-P0_cgs)
           end if
        end if
        if (i_str <= struc_len-2) then
           radius(i_str) = radius(i_str+1) + integ_parab(press(i_str),press(i_str+1),press(i_str+2), &
                inv_rho(i_str),inv_rho(i_str+1),inv_rho(i_str+2),press(i_str),press(i_str+1))/gravity
        else
           radius(i_str) = radius(i_str+1)+(1d0/rho(i_str)+1d0/rho(i_str+1))/(2d0*gravity)* &
                (press(i_str+1)-press(i_str))
        end if
     end do

     R0 = R_pl-R0
     radius = radius + R0
!!$     write(*,*) R0, P0_cgs, gravity, R_pl, press(20), rho(20)

  end if

end subroutine calc_radius

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to add Rayleigh scattering

subroutine add_rayleigh(spec,abund,lambda_angstroem,MMW,temp,press,rayleigh_kappa,struc_len,freq_len)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                         :: freq_len, struc_len
  CHARACTER*20, intent(in)                    :: spec
  DOUBLE PRECISION, intent(in)                :: lambda_angstroem(freq_len), &
       abund(struc_len), &
       MMW(struc_len), temp(struc_len), press(struc_len)
  DOUBLE PRECISION, intent(out)               :: rayleigh_kappa(freq_len,struc_len)

  ! Internal
  INTEGER                                     :: i_str, i_freq
  DOUBLE PRECISION                            :: lambda_cm(freq_len), &
       lamb_inv(freq_len), alpha_pol, lamb_inv_use
  DOUBLE PRECISION                            :: a0, a1, a2, a3, a4, a5, &
       a6, a7, luv, lir, l(freq_len), &
       d(struc_len), T(struc_len), retVal, retValMin, retValMax, mass_h2o, &
       nm1, fk, scale, mass_co2, &
       mass_o2, mass_n2, A, B, C, mass_co, nfr_co, &
       mass_ch4, nfr_ch4
  
  rayleigh_kappa = 0d0
  
  if (trim(adjustl(spec)) .EQ. 'H2') then

     ! H2 Rayleigh according to dalgarno & williams (1962)
     do i_str = 1, struc_len

        if (abund(i_str) > 1d-60) then
           rayleigh_kappa(:,i_str) = rayleigh_kappa(:,i_str) + &
                (8.14d-13/lambda_angstroem**4+1.28d-6/lambda_angstroem**6+1.61d0/lambda_angstroem**8)/2d0 &
                /1.66053892d-24*abund(i_str)
        end if
        
     end do
     
  else if (trim(adjustl(spec)) .EQ. 'He') then
     
     ! He Rayleigh scattering according to Chan & Dalgarno alphas (1965)
     lambda_cm = lambda_angstroem*1d-8
     
     do i_str = 1, struc_len
        if (abund(i_str) > 1d-60) then
           do i_freq = 1, freq_len

              if (lambda_cm(i_freq) >= 0.9110d-4) then
                 alpha_pol = 1.379
              else
                 alpha_pol = 2.265983321d0 - 3.721350022d0*lambda_cm(i_freq)/1d-4 &
                      + 3.016150391d0*(lambda_cm(i_freq)/1d-4)**2d0
              end if

              rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                   +128d0*pi**5d0/3d0/lambda_cm(i_freq)**4d0*(alpha_pol*1.482d-25)**2d0/4d0 &
                   /1.66053892d-24*abund(i_str)
           end do
        end if
     end do

  else if (trim(adjustl(spec)) .EQ. 'H2O') then

     ! For H2O Rayleigh scattering according to Harvey et al. (1998)
     a0 = 0.244257733
     a1 = 9.74634476d-3
     a2 = -3.73234996d-3
     a3 = 2.68678472d-4
     a4 = 1.58920570d-3
     a5 = 2.45934259d-3
     a6 = 0.900704920
     a7 = -1.66626219d-2
     luv = 0.2292020d0
     lir = 5.432937d0
     mass_h2o = 18d0*amu

     lambda_cm = lambda_angstroem*1d-8
     lamb_inv = 1d0/lambda_cm

     l = lambda_cm/1d-4/0.589d0
     d = MMW*amu*press/kB/temp*abund
     T = temp/273.15d0

     do i_str = 1, struc_len
        if (abund(i_str) > 1d-60) then
           do i_freq = 1, freq_len

              retVal = (a0+a1*d(i_str)+a2*T(i_str)+a3*l(i_freq)**2d0*T(i_str)+a4/l(i_freq)**2d0 &
                   + a5/(l(i_freq)**2d0-luv**2d0) + a6/(l(i_freq)**2d0-lir**2d0) + &
                   a7*d(i_str)**2d0)*d(i_str)

              retValMin = (a0+a1*d(i_str)+a2*T(i_str)+a3*(0.2d0/0.589d0)**2d0*T(i_str)+a4/(0.2d0/0.589d0)**2d0 &
                   + a5/((0.2d0/0.589d0)**2d0-luv**2d0) + a6/((0.2d0/0.589d0)**2d0-lir**2d0) + &
                   a7*d(i_str)**2d0)*d(i_str)

              retValMax = (a0+a1*d(i_str)+a2*T(i_str)+a3*(1.1d0/0.589d0)**2d0*T(i_str)+a4/(1.1d0/0.589d0)**2d0 &
                   + a5/((1.1d0/0.589d0)**2d0-luv**2d0) + a6/((1.1d0/0.589d0)**2d0-lir**2d0) + &
                   a7*d(i_str)**2d0)*d(i_str)

              if ((lambda_cm(i_freq)/1d-4 > 0.2d0) .AND. (lambda_cm(i_freq)/1d-4 < 1.1d0)) then
                 nm1 = sqrt((1d0+2d0*retVal)/(1d0-retVal))
              else if (lambda_cm(i_freq)/1d-4 >= 1.1d0) then
                 nm1 = sqrt((1.+2.*retValMax)/(1.-retValMax))
              else
                 nm1 = sqrt((1.+2.*retValMin)/(1.-retValMin))
              end if

              nm1 = nm1 - 1d0
              fk = 1.0
              rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                   + 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(d(i_str)/18d0/amu)**2d0* &
                   (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_h2o * &
                   abund(i_str)

           end do
        end if
     end do

  else if (trim(adjustl(spec)) .EQ. 'CO2') then

     ! CO2 Rayleigh scattering according to Sneep & Ubachs (2004)
     d = MMW*amu*press/kB/temp*abund

     lambda_cm = lambda_angstroem*1d-8
     lamb_inv = 1d0/lambda_cm
     mass_co2 = 44d0*amu

     do i_str = 1, struc_len
        if (abund(i_str) > 1d-60) then
           scale = d(i_str)/44d0/amu/sneep_ubachs_n
           do i_freq = 1, freq_len

              nm1 = 1d-3*1.1427d6*( 5799.25d0/max(20d0**2d0,128908.9d0**2d0-lamb_inv(i_freq)**2d0) + &
                   120.05d0/max(20d0**2d0,89223.8d0**2d0-lamb_inv(i_freq)**2d0) + &
                   5.3334d0/max(20d0**2d0,75037.5d0**2d0-lamb_inv(i_freq)**2d0) + &
                   4.3244/max(20d0**2d0,67837.7d0**2d0-lamb_inv(i_freq)**2d0) + &
                   0.1218145d-4/max(20d0**2d0,2418.136d0**2d0-lamb_inv(i_freq)**2d0))
              nm1 = nm1 * scale
              fk = 1.1364+25.3d-12*lamb_inv(i_freq)**2d0
              rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                   + 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(scale*sneep_ubachs_n)**2d0* &
                   (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_co2 * &
                   abund(i_str)

           end do
        end if
     end do

  else if (trim(adjustl(spec)) .EQ. 'O2') then

     ! O2 Rayleigh scattering according to Thalman et al. (2014).
     ! Also see their erratum!
     d = MMW*amu*press/kB/temp*abund

     lambda_cm = lambda_angstroem*1d-8
     lamb_inv = 1d0/lambda_cm
     mass_o2 = 32d0*amu

     do i_str = 1, struc_len
        if (abund(i_str) > 1d-60) then
           scale = d(i_str)/mass_o2/2.68678d19

           do i_freq = 1, freq_len

              if (lamb_inv(i_freq) > 18315d0) then
                 A = 20564.8d0
                 B = 2.480899d13
              else
                 A = 21351.1d0
                 B = 2.18567d13
              end if
              C = 4.09d9

              nm1 = 1d-8*(A+B/(C-lamb_inv(i_freq)**2d0))
              nm1 = nm1 !* scale
              fk = 1.096d0+1.385d-11*lamb_inv(i_freq)**2d0+1.448d-20*lamb_inv(i_freq)**4d0
              rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                   + 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(2.68678d19)**2d0* & !(d(i_str)/mass_o2)**2d0* &
                   (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_o2 * &
                   abund(i_str)

           end do
        end if
     end do

  else if (trim(adjustl(spec)) .EQ. 'N2') then

     ! N2 Rayleigh scattering according to Thalman et al. (2014).
     ! Also see their erratum!
     d = MMW*amu*press/kB/temp*abund

     lambda_cm = lambda_angstroem*1d-8
     lamb_inv = 1d0/lambda_cm
     mass_n2 = 34d0*amu

     do i_str = 1, struc_len
        if (abund(i_str) > 1d-60) then
           scale = d(i_str)/mass_n2/2.546899d19

           do i_freq = 1, freq_len

              if (lamb_inv(i_freq) > 4860d0) then 

                 if (lamb_inv(i_freq) > 21360d0) then
                    A = 5677.465d0
                    B = 318.81874d12
                    C = 14.4d9
                 else
                    A = 6498.2d0
                    B = 307.43305d12
                    C = 14.4d9
                 end if

                 nm1 = 1d-8*(A+B/(C-lamb_inv(i_freq)**2d0))
                 nm1 = nm1 !* scale
                 fk = 1.034d0+3.17d-12*lamb_inv(i_freq)**2d0
                 rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                      + 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(2.546899d19)**2d0* & !(d(i_str)/mass_n2)**2d0* &
                      (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_n2 * &
                      abund(i_str)

              end if

           end do
        end if
     end do

  else if (trim(adjustl(spec)) .EQ. 'CO') then

     ! CO Rayleigh scattering according to Sneep & Ubachs (2004)

     lambda_cm = lambda_angstroem*1d-8
     lamb_inv = 1d0/lambda_cm

     d = MMW*amu*press/kB/temp*abund

     do i_str = 1, struc_len

        if (abund(i_str) > 1d-60) then

           scale = d(i_str)/28d0/amu/sneep_ubachs_n
           nfr_co = d(i_str)/28d0/amu
           mass_co = 28d0*amu

           do i_freq = 1, freq_len

              lamb_inv_use = lamb_inv(i_freq)
              if (lambda_cm(i_freq)/1e-4 < 0.168d0) then
                 lamb_inv_use = 1d0/0.168d-4
              end if
              nm1 = (22851d0 + 0.456d12/(71427d0**2d0-lamb_inv_use**2d0))*1d-8
              nm1 = nm1 * scale
              fk = 1.016d0
              
              rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                   + 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(nfr_co)**2d0* &
                   (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_co * &
                   abund(i_str)

           end do
        end if

     end do

  else if (trim(adjustl(spec)) .EQ. 'CH4') then

     ! CH4 Rayleigh scattering according to Sneep & Ubachs (2004)

     lambda_cm = lambda_angstroem*1d-8
     lamb_inv = 1d0/lambda_cm

     d = MMW*amu*press/kB/temp*abund

     do i_str = 1, struc_len

        if (abund(i_str) > 1d-60) then

           scale = d(i_str)/16d0/amu/sneep_ubachs_n
           nfr_ch4 = d(i_str)/16d0/amu
           mass_ch4 = 16d0*amu

           do i_freq = 1, freq_len
              
              nm1 = (46662d0 + 4.02d-6*lamb_inv(i_freq)**2d0)*1d-8
              nm1 = nm1 * scale
              fk = 1.0
              rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                   + 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(nfr_ch4)**2d0* &
                   (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_ch4 * &
                   abund(i_str)

           end do
        end if

     end do


  end if

end subroutine add_rayleigh

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate the contribution function of the transmission spectrum

subroutine calc_transm_spec_contr(freq,total_kappa,temp,press,gravity,mmw,P0_bar,R_pl, &
     w_gauss,transm_in,scat,continuum_opa_scat,var_grav,contr_tr,freq_len,struc_len,g_len,N_species)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                         :: freq_len, struc_len, g_len, N_species
  DOUBLE PRECISION, intent(in)                :: freq(freq_len), P0_bar, R_pl
  DOUBLE PRECISION, intent(in)                :: temp(struc_len), press(struc_len), mmw(struc_len)
  DOUBLE PRECISION, intent(in)                :: total_kappa(g_len,freq_len,N_species,struc_len)

  DOUBLE PRECISION, intent(in)                :: gravity
  DOUBLE PRECISION, intent(in)                :: w_gauss(g_len), continuum_opa_scat(freq_len,struc_len)
  LOGICAL, intent(in)                         :: scat
  DOUBLE PRECISION, intent(in)                :: transm_in(freq_len)
  DOUBLE PRECISION, intent(out)               :: contr_tr(struc_len,freq_len)

  ! Internal
  DOUBLE PRECISION                            :: P0_cgs, rho(struc_len), radius(struc_len),radius_var(struc_len)
  INTEGER                                     :: i_str, i_freq, i_g, i_spec, j_str, i_leave_str
  LOGICAL                                     :: var_grav
  DOUBLE PRECISION                            :: alpha_t2(g_len,freq_len,N_species,struc_len-1)
  DOUBLE PRECISION                            :: t_graze(g_len,freq_len,N_species,struc_len), s_1, s_2, &
       t_graze_wlen_int(struc_len,freq_len), alpha_t2_scat(freq_len,struc_len-1), t_graze_scat(freq_len,struc_len), &
       total_kappa_use(g_len,freq_len,N_species,struc_len), continuum_opa_scat_use(freq_len,struc_len), &
       transm(freq_len)

  ! Convert reference pressure to cgs
  P0_cgs = P0_bar*1d6
  ! Calculate density
  rho = mmw*amu*press/kB/temp
  ! Calculate planetary radius (in cm), assuming hydrostatic equilibrium
  call calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,var_grav,radius)

  do i_leave_str = 1, struc_len

     transm = 0d0
     t_graze = 0d0
     t_graze_scat = 0d0

     continuum_opa_scat_use = continuum_opa_scat
     total_kappa_use = total_kappa
     total_kappa_use(:,:,:,i_leave_str) = 0d0
     continuum_opa_scat_use(:,i_leave_str) = 0d0
     
     ! Calc. mean free paths across grazing distances
     do i_str = 1, struc_len-1
        alpha_t2(:,:,:,i_str) = (total_kappa_use(:,:,:,i_str)*rho(i_str)+total_kappa_use(:,:,:,i_str+1)*rho(i_str+1))
     end do
     if (scat) then
        do i_str = 1, struc_len-1
           alpha_t2_scat(:,i_str) = (continuum_opa_scat_use(:,i_str)*rho(i_str)+ &
                continuum_opa_scat_use(:,i_str+1)*rho(i_str+1))
        end do
     end if

     ! Cacuclate grazing rays optical depths
     do i_str = 2, struc_len
        s_1 = sqrt(radius(1)**2d0-radius(i_str)**2d0)
        do j_str = 1, i_str-1
           if (j_str > 1) then
              s_1 = s_2
           end if
           s_2 = sqrt(radius(j_str+1)**2d0-radius(i_str)**2d0)
           t_graze(:,:,:,i_str) = t_graze(:,:,:,i_str)+alpha_t2(:,:,:,j_str)*(s_1-s_2)
        end do
     end do
     if (scat) then
        do i_str = 2, struc_len
           s_1 = sqrt(radius(1)**2d0-radius(i_str)**2d0)
           do j_str = 1, i_str-1
              if (j_str > 1) then
                 s_1 = s_2
              end if
              s_2 = sqrt(radius(j_str+1)**2d0-radius(i_str)**2d0)
              t_graze_scat(:,i_str) = t_graze_scat(:,i_str)+alpha_t2_scat(:,j_str)*(s_1-s_2)
           end do
        end do
     end if

     ! Calculate transmissions, update tau array to store these
     t_graze = exp(-t_graze)
     if (scat) then
        t_graze_scat = exp(-t_graze_scat)
     end if

     t_graze_wlen_int = 1d0
     ! Wlen (in g-space) integrate transmissions
     do i_str = 2, struc_len ! i_str=1 t_grazes are 1 anyways
        do i_spec = 1, N_species
           do i_freq = 1, freq_len
              t_graze_wlen_int(i_str,i_freq) = t_graze_wlen_int(i_str,i_freq)* &
                   SUM(t_graze(:,i_freq,i_spec,i_str)*w_gauss)
              if (scat .AND. (i_spec .EQ. 1)) then
                 t_graze_wlen_int(i_str,i_freq) = t_graze_wlen_int(i_str,i_freq)* &
                      t_graze_scat(i_freq,i_str)
              end if
           end do
        end do
     end do

     ! Get effective area fraction from transmission
     t_graze_wlen_int = 1d0-t_graze_wlen_int

     ! Caculate planets effectice area (leaving out pi, because we want the radius in the end)
     do i_freq = 1, freq_len
        do i_str = 2, struc_len
           transm(i_freq) = transm(i_freq)+(t_graze_wlen_int(i_str-1,i_freq)*radius(i_str-1)+ &
                t_graze_wlen_int(i_str,i_freq)*radius(i_str))*(radius(i_str-1)-radius(i_str))
        end do
     end do
     ! Get radius
     transm = transm+radius(struc_len)**2d0
     contr_tr(i_leave_str,:) = transm_in - transm

     write(*,*) i_leave_str
     
  end do

  do i_freq = 1, freq_len
     contr_tr(:,i_freq) = contr_tr(:,i_freq)/SUM(contr_tr(:,i_freq))
  end do
  
end subroutine calc_transm_spec_contr

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Function to calc higher order integ.

function integ_parab(x,y,z,fx,fy,fz,a,b)

  implicit none
  ! I/O
  double precision :: x,y,z,fx,fy,fz,a,b
  double precision :: integ_parab
  ! Internal
  double precision :: c1,c2,c3

  c3 = ((fz-fy)/(z-y)-(fz-fx)/(z-x))/(y-x)
  c2 = (fz-fx)/(z-x)-c3*(z+x)
  c1 = fx-c2*x-c3*x**2d0

  integ_parab = c1*(b-a)+c2*(b**2d0-a**2d0)/2d0+c3*(b**3d0-a**3d0)/3d0

end function integ_parab
  
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################


!!$ Subroutine to calculate cloud opacities

subroutine calc_cloud_opas(rho,rho_p,cloud_mass_fracs,r_g,sigma_n,cloud_rad_bins,cloud_radii,cloud_lambdas, &
     cloud_specs_abs_opa,cloud_specs_scat_opa,cloud_aniso,cloud_abs_opa_TOT,cloud_scat_opa_TOT, &
     cloud_red_fac_aniso_TOT,struc_len,N_cloud_spec,N_cloud_rad_bins, N_cloud_lambda_bins)
  
  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in) :: struc_len, N_cloud_spec, N_cloud_rad_bins, N_cloud_lambda_bins
  DOUBLE PRECISION, intent(in) :: rho(struc_len), rho_p(N_cloud_spec)
  DOUBLE PRECISION, intent(in) :: cloud_mass_fracs(struc_len,N_cloud_spec),r_g(struc_len,N_cloud_spec) 
  DOUBLE PRECISION, intent(in) :: sigma_n
  DOUBLE PRECISION, intent(in) :: cloud_rad_bins(N_cloud_rad_bins+1), cloud_radii(N_cloud_rad_bins), &
       cloud_lambdas(N_cloud_lambda_bins)
  DOUBLE PRECISION, intent(in) :: cloud_specs_abs_opa(N_cloud_rad_bins,N_cloud_lambda_bins,N_cloud_spec), &
       cloud_specs_scat_opa(N_cloud_rad_bins,N_cloud_lambda_bins,N_cloud_spec), &
       cloud_aniso(N_cloud_rad_bins,N_cloud_lambda_bins,N_cloud_spec)
  DOUBLE PRECISION, intent(out) :: cloud_abs_opa_TOT(N_cloud_lambda_bins,struc_len), &
       cloud_scat_opa_TOT(N_cloud_lambda_bins,struc_len), &
       cloud_red_fac_aniso_TOT(N_cloud_lambda_bins,struc_len)

  ! internal
  INTEGER :: i_struc, i_spec, i_rad, i_lamb
  DOUBLE PRECISION :: N, dndr(N_cloud_rad_bins), integrand_abs(N_cloud_rad_bins), &
       integrand_scat(N_cloud_rad_bins), add_abs, add_scat, integrand_aniso(N_cloud_rad_bins), add_aniso
  !~~~~~~~~~~~~~~~~

  cloud_abs_opa_TOT = 0d0
  cloud_scat_opa_TOT = 0d0
  cloud_red_fac_aniso_TOT = 0d0

  do i_struc = 1, struc_len
     do i_spec = 1, N_cloud_spec

           do i_lamb = 1, N_cloud_lambda_bins

              N = 3d0*cloud_mass_fracs(i_struc,i_spec)*rho(i_struc)/4d0/pi/rho_p(i_spec)/ &
                   r_g(i_struc,i_spec)**3d0*exp(-9d0/2d0*log(sigma_n)**2d0)

              dndr = N/(cloud_radii*sqrt(2d0*pi)*log(sigma_n))* &
                   exp(-log(cloud_radii/r_g(i_struc,i_spec))**2d0/(2d0*log(sigma_n)**2d0))

              integrand_abs = 4d0*pi/3d0*cloud_radii**3d0*rho_p(i_spec)*dndr* &
                   cloud_specs_abs_opa(:,i_lamb,i_spec)
              integrand_scat = 4d0*pi/3d0*cloud_radii**3d0*rho_p(i_spec)*dndr* &
                   cloud_specs_scat_opa(:,i_lamb,i_spec)
              integrand_aniso = integrand_scat*(1d0-cloud_aniso(:,i_lamb,i_spec))

              add_abs = sum(integrand_abs*(cloud_rad_bins(2:N_cloud_rad_bins+1)- &
                   cloud_rad_bins(1:N_cloud_rad_bins)))
              cloud_abs_opa_TOT(i_lamb,i_struc) = cloud_abs_opa_TOT(i_lamb,i_struc) + &
                   add_abs
                   
              add_scat = sum(integrand_scat*(cloud_rad_bins(2:N_cloud_rad_bins+1)- &
                   cloud_rad_bins(1:N_cloud_rad_bins)))
              cloud_scat_opa_TOT(i_lamb,i_struc) = cloud_scat_opa_TOT(i_lamb,i_struc) + &
                   add_scat

              add_aniso = sum(integrand_aniso*(cloud_rad_bins(2:N_cloud_rad_bins+1)- &
                   cloud_rad_bins(1:N_cloud_rad_bins)))
              cloud_red_fac_aniso_TOT(i_lamb,i_struc) = cloud_red_fac_aniso_TOT(i_lamb,i_struc) + &
                   add_aniso
              
           end do
                   
     end do

     do i_lamb = 1, N_cloud_lambda_bins
        if (cloud_scat_opa_TOT(i_lamb,i_struc) > 1d-200) then
           cloud_red_fac_aniso_TOT(i_lamb,i_struc) = cloud_red_fac_aniso_TOT(i_lamb,i_struc)/ &
                     cloud_scat_opa_TOT(i_lamb,i_struc)
        else
           cloud_red_fac_aniso_TOT(i_lamb,i_struc) = 0d0
        end if
     end do
     
     cloud_abs_opa_TOT(:,i_struc) = cloud_abs_opa_TOT(:,i_struc)/rho(i_struc)
     cloud_scat_opa_TOT(:,i_struc) = cloud_scat_opa_TOT(:,i_struc)/rho(i_struc)
     
  end do

end subroutine calc_cloud_opas

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Interpolate cloud opacities to actual radiative transfer wavelength grid

subroutine interp_integ_cloud_opas(cloud_abs_opa_TOT,cloud_scat_opa_TOT, &
     cloud_red_fac_aniso_TOT,cloud_lambdas,HIT_border_freqs,HIT_kappa_tot_g_approx, &
     HIT_kappa_tot_g_approx_scat,red_fac_aniso_final, HIT_kappa_tot_g_approx_scat_unred, &
     N_cloud_lambda_bins,struc_len,HIT_coarse_borders)

  use constants_block
  implicit none
  ! I/O
  INTEGER, intent(in)           :: N_cloud_lambda_bins,struc_len,HIT_coarse_borders
  DOUBLE PRECISION, intent(in)  :: cloud_abs_opa_TOT(N_cloud_lambda_bins,struc_len), &
       cloud_scat_opa_TOT(N_cloud_lambda_bins,struc_len), &
       cloud_red_fac_aniso_TOT(N_cloud_lambda_bins,struc_len), cloud_lambdas(N_cloud_lambda_bins), &
       HIT_border_freqs(HIT_coarse_borders)
  DOUBLE PRECISION, intent(out) :: HIT_kappa_tot_g_approx(HIT_coarse_borders-1,struc_len), &
       HIT_kappa_tot_g_approx_scat(HIT_coarse_borders-1,struc_len), &
       red_fac_aniso_final(HIT_coarse_borders-1,struc_len), &
       HIT_kappa_tot_g_approx_scat_unred(HIT_coarse_borders-1,struc_len)
  
  ! internal
  DOUBLE PRECISION :: kappa_integ(struc_len), kappa_scat_integ(struc_len), red_fac_aniso_integ(struc_len), &
       kappa_tot_integ(HIT_coarse_borders-1,struc_len), kappa_tot_scat_integ(HIT_coarse_borders-1,struc_len)
  INTEGER          :: HIT_i_lamb
  DOUBLE PRECISION :: HIT_border_lamb(HIT_coarse_borders)
  INTEGER          :: intp_index_small_min, intp_index_small_max, i_lamb, i_struc, &
       new_small_ind, i_ng

  HIT_kappa_tot_g_approx = 0d0
  HIT_kappa_tot_g_approx_scat = 0d0
  HIT_kappa_tot_g_approx_scat_unred = 0d0
  
  
  HIT_border_lamb = c_l/HIT_border_freqs
  red_fac_aniso_final = 0d0

  kappa_tot_integ = 0d0
  kappa_tot_scat_integ = 0d0

  do HIT_i_lamb = 1, HIT_coarse_borders-1

     intp_index_small_min = MIN(MAX(INT((log10(HIT_border_lamb(HIT_i_lamb))-log10(cloud_lambdas(1))) / &
          log10(cloud_lambdas(N_cloud_lambda_bins)/cloud_lambdas(1))*DBLE(N_cloud_lambda_bins-1) &
          +1d0),1),N_cloud_lambda_bins-1)
     
     intp_index_small_max = MIN(MAX(INT((log10(HIT_border_lamb(HIT_i_lamb+1))-log10(cloud_lambdas(1))) / &
          log10(cloud_lambdas(N_cloud_lambda_bins)/cloud_lambdas(1))*DBLE(N_cloud_lambda_bins-1) &
          +1d0),1),N_cloud_lambda_bins-1)

     kappa_integ = 0d0
     kappa_scat_integ = 0d0
     red_fac_aniso_integ = 0d0
     
     if ((intp_index_small_max-intp_index_small_min) .EQ. 0) then

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_abs_opa_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),HIT_border_lamb(HIT_i_lamb+1),kappa_integ)

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_scat_opa_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),HIT_border_lamb(HIT_i_lamb+1),kappa_scat_integ)

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_red_fac_aniso_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),HIT_border_lamb(HIT_i_lamb+1),red_fac_aniso_integ)
                
     else if ((intp_index_small_max-intp_index_small_min) .EQ. 1) then

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_abs_opa_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),cloud_lambdas(intp_index_small_min+1),kappa_integ)
        
        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_scat_opa_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),cloud_lambdas(intp_index_small_min+1),kappa_scat_integ)

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_red_fac_aniso_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),cloud_lambdas(intp_index_small_min+1),red_fac_aniso_integ)

        call integ_kaps(intp_index_small_max,N_cloud_lambda_bins,struc_len,cloud_abs_opa_TOT, &
             cloud_lambdas,cloud_lambdas(intp_index_small_max),HIT_border_lamb(HIT_i_lamb+1),kappa_integ)
        
        call integ_kaps(intp_index_small_max,N_cloud_lambda_bins,struc_len,cloud_scat_opa_TOT, &
             cloud_lambdas,cloud_lambdas(intp_index_small_max),HIT_border_lamb(HIT_i_lamb+1),kappa_scat_integ)

        call integ_kaps(intp_index_small_max,N_cloud_lambda_bins,struc_len,cloud_red_fac_aniso_TOT, &
             cloud_lambdas,cloud_lambdas(intp_index_small_max),HIT_border_lamb(HIT_i_lamb+1),red_fac_aniso_integ)
                
     else

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_abs_opa_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),cloud_lambdas(intp_index_small_min+1),kappa_integ)
        
        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_scat_opa_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),cloud_lambdas(intp_index_small_min+1),kappa_scat_integ)

        call integ_kaps(intp_index_small_min,N_cloud_lambda_bins,struc_len,cloud_red_fac_aniso_TOT, &
             cloud_lambdas,HIT_border_lamb(HIT_i_lamb),cloud_lambdas(intp_index_small_min+1),red_fac_aniso_integ)

        new_small_ind = intp_index_small_min+1
        do while (intp_index_small_max-new_small_ind .NE. 0)
           
           call integ_kaps(new_small_ind,N_cloud_lambda_bins,struc_len,cloud_abs_opa_TOT, &
                cloud_lambdas,cloud_lambdas(new_small_ind),cloud_lambdas(new_small_ind+1),kappa_integ)
        
           call integ_kaps(new_small_ind,N_cloud_lambda_bins,struc_len,cloud_scat_opa_TOT, &
                cloud_lambdas,cloud_lambdas(new_small_ind),cloud_lambdas(new_small_ind+1),kappa_scat_integ)

           call integ_kaps(new_small_ind,N_cloud_lambda_bins,struc_len,cloud_red_fac_aniso_TOT, &
                cloud_lambdas,cloud_lambdas(new_small_ind),cloud_lambdas(new_small_ind+1),red_fac_aniso_integ)

           new_small_ind = new_small_ind+1
           
        end do

        call integ_kaps(intp_index_small_max,N_cloud_lambda_bins,struc_len,cloud_abs_opa_TOT, &
             cloud_lambdas,cloud_lambdas(intp_index_small_max),HIT_border_lamb(HIT_i_lamb+1),kappa_integ)
        
        call integ_kaps(intp_index_small_max,N_cloud_lambda_bins,struc_len,cloud_scat_opa_TOT, &
             cloud_lambdas,cloud_lambdas(intp_index_small_max),HIT_border_lamb(HIT_i_lamb+1),kappa_scat_integ)

        call integ_kaps(intp_index_small_max,N_cloud_lambda_bins,struc_len,cloud_red_fac_aniso_TOT, &
             cloud_lambdas,cloud_lambdas(intp_index_small_max),HIT_border_lamb(HIT_i_lamb+1),red_fac_aniso_integ)
                
     end if

     kappa_integ = kappa_integ/(HIT_border_lamb(HIT_i_lamb+1)-HIT_border_lamb(HIT_i_lamb))
     kappa_scat_integ = kappa_scat_integ/(HIT_border_lamb(HIT_i_lamb+1)-HIT_border_lamb(HIT_i_lamb))
     red_fac_aniso_integ = red_fac_aniso_integ/(HIT_border_lamb(HIT_i_lamb+1)-HIT_border_lamb(HIT_i_lamb))
     
     kappa_tot_integ(HIT_i_lamb,:) = kappa_integ
     kappa_tot_scat_integ(HIT_i_lamb,:) = kappa_scat_integ

     HIT_kappa_tot_g_approx(HIT_i_lamb,:) = HIT_kappa_tot_g_approx(HIT_i_lamb,:) + &
          kappa_integ
     HIT_kappa_tot_g_approx_scat(HIT_i_lamb,:) = HIT_kappa_tot_g_approx_scat(HIT_i_lamb,:) + &
          kappa_integ + kappa_scat_integ*red_fac_aniso_integ
     HIT_kappa_tot_g_approx_scat_unred(HIT_i_lamb,:) = HIT_kappa_tot_g_approx_scat_unred(HIT_i_lamb,:) + &
          kappa_integ + kappa_scat_integ

     red_fac_aniso_final(HIT_i_lamb,:) = red_fac_aniso_final(HIT_i_lamb,:) + red_fac_aniso_integ
     
  end do

end subroutine interp_integ_cloud_opas

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

subroutine integ_kaps(intp_ind,N_cloud_lambda_bins,struc_len,kappa,lambda,l_bord1,l_bord2,kappa_integ)
  implicit none
  INTEGER, intent(in) :: intp_ind,N_cloud_lambda_bins,struc_len
  DOUBLE PRECISION, intent(in) :: lambda(N_cloud_lambda_bins), kappa(N_cloud_lambda_bins,struc_len)
  DOUBLE PRECISION, intent(in) :: l_bord1,l_bord2
  DOUBLE PRECISION, intent(out) :: kappa_integ(struc_len)

  ! This subroutine calculates the integral of a linearly interpolated function kappa.
  
  kappa_integ = kappa_integ + kappa(intp_ind,:)*(l_bord2-l_bord1) + (kappa(intp_ind+1,:)-kappa(intp_ind,:))/ &
       (lambda(intp_ind+1)-lambda(intp_ind))* &
       0.5d0*((l_bord2-lambda(intp_ind))**2d0-(l_bord1-lambda(intp_ind))**2d0)

end subroutine integ_kaps

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine get_rg_N(gravity,rho,rho_p,temp,MMW,frain,cloud_mass_fracs, &
     sigma_n,Kzz,r_g,struc_len,N_cloud_spec)

  use constants_block
  implicit none
  ! I/O
  INTEGER, intent(in)  :: struc_len, N_cloud_spec
  DOUBLE PRECISION, intent(in) :: gravity, rho(struc_len), rho_p(N_cloud_spec), temp(struc_len), &
       MMW(struc_len), frain, cloud_mass_fracs(struc_len,N_cloud_spec), &
       sigma_n, Kzz(struc_len)
  DOUBLE PRECISION, intent(out) :: r_g(struc_len,N_cloud_spec)
  ! Internal
  INTEGER, parameter :: N_fit = 100
  INTEGER          :: i_str, i_size, i_spec, i_rad
  DOUBLE PRECISION :: turbulent_settling_speed_ret1, turbulent_settling_speed_ret2, &
       bisect_particle_rad
  DOUBLE PRECISION :: w_star(struc_len), H(struc_len)
  DOUBLE PRECISION :: r_w(struc_len,N_cloud_spec), alpha(struc_len,N_cloud_spec)
  DOUBLE PRECISION :: rad(N_fit), vel(N_fit), f_fill(N_cloud_spec)
  DOUBLE PRECISION :: a, b

  H = kB*temp/(MMW*amu*gravity)
  w_star = Kzz/H
  
  f_fill = 1d0
  
  do i_str = 1, struc_len
     do i_spec = 1, N_cloud_spec
        r_w(i_str,i_spec) = bisect_particle_rad(1d-16,1d2,gravity,rho(i_str), &
             rho_p(i_spec),temp(i_str),MMW(i_str),w_star(i_str))
        if (r_w(i_str,i_spec) > 1d-16) then
           if (frain > 1d0) then
              do i_rad = 1, N_fit
                 rad(i_rad) = r_w(i_str,i_spec)/max(sigma_n,1.1d0) + &
                      (r_w(i_str,i_spec)-r_w(i_str,i_spec)/max(sigma_n,1.1d0))* &
                      DBLE(i_rad-1)/DBLE(N_fit-1)
                 call turbulent_settling_speed(rad(i_rad),gravity,rho(i_str),rho_p(i_spec),temp(i_str), &
                      MMW(i_str),vel(i_rad))
              end do
           else
              do i_rad = 1, N_fit
                 rad(i_rad) = r_w(i_str,i_spec) + (r_w(i_str,i_spec)*max(sigma_n,1.1d0)- &
                      r_w(i_str,i_spec))* &
                      DBLE(i_rad-1)/DBLE(N_fit-1)
                 call turbulent_settling_speed(rad(i_rad),gravity,rho(i_str),rho_p(i_spec),temp(i_str), &
                      MMW(i_str),vel(i_rad))
              end do
           end if

           call fit_linear(log(rad), log(vel/w_star(i_str)), N_fit, a, b)
           
           alpha(i_str,i_spec) = b
           r_w(i_str,i_spec) = exp(-a/b)
           r_g(i_str,i_spec) = r_w(i_str,i_spec) * frain**(1d0/alpha(i_str,i_spec))* &
                exp(-(alpha(i_str,i_spec)+6d0)/2d0*log(sigma_n)**2d0)
        else
           r_g(i_str,i_spec) = 1d-17
           alpha(i_str,i_spec) = 1d0
        end if
     end do

  end do
  
end subroutine get_rg_N

subroutine turbulent_settling_speed(x,gravity,rho,rho_p,temp,MMW,turbulent_settling_speed_ret)

  use constants_block
  implicit none
  DOUBLE PRECISION    :: turbulent_settling_speed_ret
  DOUBLE PRECISION    :: x,gravity,rho,rho_p,temp,MMW
  DOUBLE PRECISION, parameter :: d = 2.827d-8, epsilon = 59.7*kB
  DOUBLE PRECISION    :: N_Knudsen, psi, eta, CdNreSq, Nre, Cd, v_settling_visc

  
  N_Knudsen = MMW*amu/(pi*rho*d**2d0*x)
  psi = 1d0 + N_Knudsen*(1.249d0+0.42d0*exp(-0.87d0*N_Knudsen))
  eta = 15d0/16d0*sqrt(pi*2d0*amu*kB*temp)/(pi*d**2d0)*(kB*temp/epsilon)**0.16d0/1.22d0
  CdNreSq = 32d0*rho*gravity*x**3d0*(rho_p-rho)/(3d0*eta**2d0)
  Nre = exp(-2.7905d0+0.9209d0*log(CdNreSq)-0.0135d0*log(CdNreSq)**2d0)
  if (Nre < 1d0) then
     Cd = 24d0
  else if (Nre > 1d3) then
     Cd = 0.45d0
  else
     Cd = CdNreSq/Nre**2d0
  end if
  v_settling_visc = 2d0*x**2d0*(rho_p-rho)*psi*gravity/(9d0*eta)
  turbulent_settling_speed_ret = psi*sqrt(8d0*gravity*x*(rho_p-rho)/(3d0*Cd*rho))
  if ((Nre < 1d0) .AND. (v_settling_visc < turbulent_settling_speed_ret)) THEN
     turbulent_settling_speed_ret = v_settling_visc
  end if
  
end subroutine turbulent_settling_speed

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Function to find the particle radius, using a simple bisection method.

function bisect_particle_rad(x1,x2,gravity,rho,rho_p,temp,MMW,w_star)

  implicit none
  INTEGER, parameter :: ITMAX = 1000
  DOUBLE PRECISION :: gravity,rho,rho_p,temp,MMW,w_star  
  DOUBLE PRECISION :: bisect_particle_rad,x1,x2
  INTEGER :: iter
  DOUBLE PRECISION :: a,b,c,fa,fb,fc,del

  a=x1
  b=x2
  call turbulent_settling_speed(a,gravity,rho,rho_p,temp,MMW,fa)
  fa = fa - w_star
  call turbulent_settling_speed(b,gravity,rho,rho_p,temp,MMW,fb)
  fb = fb - w_star
  
  if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
     write(*,*) 'Warning: root not bracketed.'
     bisect_particle_rad = 1d-17
     return 
  end if

  do iter=1,ITMAX

     if (abs(log10(a/b)) > 1d0) then
        c = 1e1**(log10(a*b)/2d0)
     else
        c = (a+b)/2d0
     end if
     
     call turbulent_settling_speed(c,gravity,rho,rho_p,temp,MMW,fc)
     fc = fc - w_star
     
     if (((fc > 0d0) .and. (fa > 0d0)) .OR. ((fc < 0d0) .and. (fa < 0d0))) then
        del = 2d0*abs(a-c)/(a+b)
        a = c
        fa = fc
     else
        del = 2d0*abs(b-c)/(a+b)
        b = c
        fb = fc
     end if

     if (abs(del) .lt. 1d-9) then
        exit
     end if

  end do

  if (iter == ITMAX) then
     write(*,*) 'warning: maximum number of bisection root iterations reached!'
  end if

  bisect_particle_rad = c
  return

end function bisect_particle_rad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Subroutine to calculate slope and y-axis intercept of x,y data,
! assuming zero error on data.

SUBROUTINE fit_linear(x, y, ndata, a, b)

  implicit none
  INTEGER :: ndata
  DOUBLE PRECISION :: x(ndata), y(ndata)
  DOUBLE PRECISION :: a, b

  b = (sum(x)*sum(y)/dble(ndata) - sum(x*y))/ &
       (sum(x)**2d0/dble(ndata) - sum(x**2d0))
  a = sum(y-b*x)/dble(ndata)

end SUBROUTINE fit_linear

