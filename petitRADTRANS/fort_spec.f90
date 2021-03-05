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
  DOUBLE PRECISION,parameter      :: G = 6.674d-8, M_jup = 1.89813d30, deg = Pi/1.8d2
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

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate tau_scat with 2nd order accuracy

subroutine calc_tau_g_tot_ck_scat(gravity,press,total_kappa_in,do_scat_emis, &
     continuum_opa_scat_emis,struc_len,freq_len,g_len,tau,photon_destruction_prob)

  use constants_block
  implicit none

  ! I/O
  INTEGER, parameter                           :: N_species = 1
  INTEGER, intent(in)                          :: struc_len, freq_len, g_len
  DOUBLE PRECISION, intent(in)                 :: total_kappa_in(g_len,freq_len,N_species,struc_len)
  DOUBLE PRECISION, intent(in)                 :: gravity, press(struc_len)
  LOGICAL, intent(in)                          :: do_scat_emis
  DOUBLE PRECISION, intent(in)                 :: continuum_opa_scat_emis(freq_len,struc_len)
  DOUBLE PRECISION, intent(out)                :: tau(g_len,freq_len,N_species,struc_len), &
       photon_destruction_prob(g_len,freq_len,struc_len)
  ! internal
  integer                                      :: i_struc, i_freq, i_g, i_spec
  DOUBLE PRECISION                             :: del_tau_lower_ord, &
       gamma_second(g_len,freq_len,N_species), f_second, kappa_i(g_len,freq_len,N_species), &
       kappa_im(g_len,freq_len,N_species), kappa_ip(g_len,freq_len,N_species)
  DOUBLE PRECISION                             :: total_kappa(g_len,freq_len,N_species,struc_len)
  LOGICAL                                      :: second_order
  !~~~~~~~~~~~~~

  tau = 0d0
  second_order = .FALSE.

  total_kappa = total_kappa_in

  if (do_scat_emis) then
     do i_g = 1, g_len
        total_kappa(i_g,:,1,:) = total_kappa(i_g,:,1,:) + &
             continuum_opa_scat_emis(:,:)
        photon_destruction_prob(i_g,:,:) = continuum_opa_scat_emis(:,:) / &
             total_kappa(i_g,:,1,:)
     end do
     photon_destruction_prob = 1d0 - photon_destruction_prob
  else
     photon_destruction_prob = 1d0
  end if

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

end subroutine calc_tau_g_tot_ck_scat

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

subroutine calc_kappa_rosseland(total_kappa, temp, w_gauss, border_freqs, &
     do_scat_emis, continuum_opa_scat_emis, &
     g_len, freq_len, struc_len, freq_len_p_1, kappa_rosse)

  implicit none

  integer,          intent(in)  :: g_len, freq_len, struc_len, freq_len_p_1
  double precision, intent(in)  :: total_kappa(g_len, freq_len, struc_len)
  double precision, intent(in)  :: border_freqs(freq_len_p_1)
  double precision, intent(in)  :: temp(struc_len), w_gauss(g_len)
  LOGICAL, intent(in)           :: do_scat_emis
  DOUBLE PRECISION, intent(in)  :: continuum_opa_scat_emis(freq_len,struc_len)
  double precision, intent(out) :: kappa_rosse(struc_len)

  double precision              :: total_kappa_use(g_len, freq_len, struc_len)

  integer                       :: i_struc, i_g

  if (do_scat_emis) then
     do i_g = 1, g_len
        total_kappa_use(i_g,:,:) = total_kappa(i_g,:,:) + continuum_opa_scat_emis
     end do
  else
     total_kappa_use = total_kappa
  end if

  do i_struc = 1, struc_len
     call calc_rosse_opa(total_kappa_use(:,:,i_struc), border_freqs, temp(i_struc), &
          g_len, freq_len+1, &
          kappa_rosse(i_struc), w_gauss)
  end do

end subroutine calc_kappa_rosseland

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

subroutine calc_kappa_planck(total_kappa, temp, w_gauss, border_freqs, &
     do_scat_emis, continuum_opa_scat_emis, &
     g_len, freq_len, struc_len, freq_len_p_1, kappa_planck)

  implicit none

  integer,          intent(in)  :: g_len, freq_len, struc_len, freq_len_p_1
  double precision, intent(in)  :: total_kappa(g_len, freq_len, struc_len)
  double precision, intent(in)  :: border_freqs(freq_len_p_1)
  double precision, intent(in)  :: temp(struc_len), w_gauss(g_len)
  LOGICAL, intent(in)           :: do_scat_emis
  DOUBLE PRECISION, intent(in)  :: continuum_opa_scat_emis(freq_len,struc_len)
  double precision, intent(out) :: kappa_planck(struc_len)

  double precision              :: total_kappa_use(g_len, freq_len, struc_len)

  integer                       :: i_struc, i_g

  if (do_scat_emis) then
     do i_g = 1, g_len
        total_kappa_use(i_g,:,:) = total_kappa(i_g,:,:) + continuum_opa_scat_emis
     end do
  else
     total_kappa_use = total_kappa
  end if

  do i_struc = 1, struc_len
     call calc_planck_opa(total_kappa_use(:,:,i_struc), border_freqs, temp(i_struc), &
          g_len, freq_len+1, &
          kappa_planck(i_struc), w_gauss)
  end do

end subroutine calc_kappa_planck

!!$ #########################################################################
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
              contr_em(i_str,i_freq) = contr_em(i_str,i_freq) &
                  + (r(i_str)+r(i_str+1)) * &
                   (transm_all_loc(i_str)-transm_all_loc(i_str+1)) &
                   *mu(i_mu)*w_gauss_mu(i_mu)
           end if
        end do
        flux_mu(i_freq) = flux_mu(i_freq) + r(struc_len)*transm_all_loc(struc_len)
        if (contribution) then
           contr_em(struc_len,i_freq) = contr_em(struc_len,i_freq) + 2d0*r(struc_len)* &
                transm_all_loc(struc_len)*mu(i_mu)*w_gauss_mu(i_mu)
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

subroutine calc_transm_spec(total_kappa_in,temp,press,gravity,mmw,P0_bar,R_pl, &
     w_gauss,scat,continuum_opa_scat,var_grav,transm,freq_len,struc_len,g_len,N_species)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                         :: freq_len, struc_len, g_len, N_species
  DOUBLE PRECISION, intent(in)                :: P0_bar, R_pl
  DOUBLE PRECISION, intent(in)                :: temp(struc_len), press(struc_len), mmw(struc_len)
  DOUBLE PRECISION, intent(in)                :: total_kappa_in(g_len,freq_len,N_species,struc_len)

  DOUBLE PRECISION, intent(in)                :: gravity
  DOUBLE PRECISION, intent(in)                :: w_gauss(g_len), continuum_opa_scat(freq_len,struc_len)
  LOGICAL, intent(in)                         :: scat !, contribution
  LOGICAL, intent(in)                         :: var_grav

  DOUBLE PRECISION, intent(out)               :: transm(freq_len) !, contr_tr(struc_len,freq_len)

  ! Internal
  DOUBLE PRECISION                            :: P0_cgs, rho(struc_len), radius(struc_len), &
        total_kappa(g_len,freq_len,N_species,struc_len)
  INTEGER                                     :: i_str, i_freq, i_g, i_spec, j_str
  LOGICAL                                     :: rad_neg
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
  call calc_radius(struc_len,press,gravity,rho,P0_cgs,R_pl,var_grav,radius)

  rad_neg = .FALSE.
  do i_str = struc_len, 1, -1
     if (radius(i_str) < 0d0) then
        rad_neg = .TRUE.
        radius(i_str) = radius(i_str+1)
     end if
  end do
  if (rad_neg) then
     write(*,*) 'pRT: negative radius corretion applied!'
  end if

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

!!$  call calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,.FALSE.,radius)
!!$  call calc_radius(struc_len,temp,press,gravity,mmw,rho,P0_cgs,R_pl,.TRUE., radius_var)
!!$  open(unit=10,file='rad_test.dat')
!!$  do i_str = 1, struc_len
!!$     write(10,*) press(i_str)*1d-6, radius(i_str)/R_jup, radius_var(i_str)/R_jup
!!$  end do
!!$  close(10)

end subroutine calc_transm_spec

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to calculate the radius from the pressure grid

subroutine calc_radius(struc_len,press,gravity,rho,P0_cgs, &
     R_pl,var_grav,radius)

  implicit none
  ! I/O
  INTEGER, intent(in)                         :: struc_len
  DOUBLE PRECISION, intent(in)                :: P0_cgs
  DOUBLE PRECISION, intent(in)                :: press(struc_len), &
       rho(struc_len)
  DOUBLE PRECISION, intent(in)                :: gravity, R_pl
  LOGICAL, intent(in)                         :: var_grav
  DOUBLE PRECISION, intent(out)               :: radius(struc_len)

  ! Internal
  INTEGER                                     :: i_str
  DOUBLE PRECISION                            :: R0, inv_rho(struc_len), integ_parab

  inv_rho = 1d0/rho

  radius = 0d0
  R0=0d0
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
  DOUBLE PRECISION, intent(in)                :: lambda_angstroem(freq_len), abund(struc_len), &
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

              retVal = 24d0*pi**3d0*lamb_inv(i_freq)**4d0/(d(i_str)/18d0/amu)**2d0* &
                   (((nm1+1d0)**2d0-1d0)/((nm1+1d0)**2d0+2d0))**2d0*fk / mass_h2o * &
                   abund(i_str)

              if (.NOT. (retVal .NE. retVal)) then
                 rayleigh_kappa(i_freq,i_str) = rayleigh_kappa(i_freq,i_str) &
                      + retVal
              end if

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

subroutine calc_transm_spec_contr(total_kappa,temp,press,gravity,mmw,P0_bar,R_pl, &
     w_gauss,transm_in,scat,continuum_opa_scat,var_grav,contr_tr,freq_len,struc_len,g_len,N_species)

  use constants_block
  implicit none

  ! I/O
  INTEGER, intent(in)                         :: freq_len, struc_len, g_len, N_species
  DOUBLE PRECISION, intent(in)                :: P0_bar, R_pl
  DOUBLE PRECISION, intent(in)                :: temp(struc_len), press(struc_len), mmw(struc_len)
  DOUBLE PRECISION, intent(in)                :: total_kappa(g_len,freq_len,N_species,struc_len)

  DOUBLE PRECISION, intent(in)                :: gravity
  DOUBLE PRECISION, intent(in)                :: w_gauss(g_len), continuum_opa_scat(freq_len,struc_len)
  LOGICAL, intent(in)                         :: scat
  DOUBLE PRECISION, intent(in)                :: transm_in(freq_len)
  LOGICAL, intent(in)                         :: var_grav
  DOUBLE PRECISION, intent(out)               :: contr_tr(struc_len,freq_len)

  ! Internal
  DOUBLE PRECISION                            :: P0_cgs, rho(struc_len), radius(struc_len)
  INTEGER                                     :: i_str, i_freq,  i_spec, j_str, i_leave_str
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
  call calc_radius(struc_len,press,gravity,rho,P0_cgs,R_pl,var_grav,radius)

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
  INTEGER :: i_struc, i_spec, i_lamb
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
  INTEGER          :: intp_index_small_min, intp_index_small_max, &
       new_small_ind

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
       MMW(struc_len), frain(N_cloud_spec), cloud_mass_fracs(struc_len,N_cloud_spec), &
       sigma_n, Kzz(struc_len)
  DOUBLE PRECISION, intent(out) :: r_g(struc_len,N_cloud_spec)
  ! Internal
  INTEGER, parameter :: N_fit = 100
  INTEGER          :: i_str, i_spec, i_rad
  DOUBLE PRECISION :: bisect_particle_rad
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
           if (frain(i_spec) > 1d0) then
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
           r_g(i_str,i_spec) = r_w(i_str,i_spec) * frain(i_spec)**(1d0/alpha(i_str,i_spec))* &
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
     !write(*,*) 'warning: root must be bracketed for zbrent'
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Subroutine to randomly correlate the opacities

subroutine combine_opas_sample_ck(line_struc_kappas, g_gauss, weights, &
     nsample, g_len, freq_len, N_species, struc_len, line_struc_kappas_out)

  implicit none

  INTEGER, INTENT(IN)          :: nsample, g_len, freq_len, N_species, struc_len
  DOUBLE PRECISION, INTENT(IN) :: line_struc_kappas(g_len, freq_len, &
       N_species, struc_len), g_gauss(g_len), weights(g_len)
  DOUBLE PRECISION, INTENT(OUT) :: line_struc_kappas_out(g_len, freq_len, &
       struc_len)
  ! Internal
!!$    INTEGER          :: i_freq, i_spec, i_struc, inds_avail(48), &
  INTEGER          :: i_freq, i_spec, i_struc, inds_avail(32), &
!!$  INTEGER          :: i_freq, i_spec, i_struc, inds_avail(16), &
       ind_use(nsample), i_samp, intpint(g_len), i_g
!  DOUBLE PRECISION :: r_index(nsample,freq_len, &
  !       N_species, struc_len), weights_use(g_len), g_sample(nsample)
  DOUBLE PRECISION :: r_index(nsample), weights_use(g_len), g_sample(nsample)
  DOUBLE PRECISION :: sampled_opa_weights(nsample, 2, freq_len, struc_len), &
       cum_sum, k_min(freq_len, struc_len), k_max(freq_len, struc_len), &
       g_final(nsample+2), k_final(nsample+2)

!  DOUBLE PRECISION :: time_test, t1, t2, t0
  DOUBLE PRECISION :: threshold(freq_len, struc_len)
  INTEGER          :: take_spec(freq_len, struc_len), take_spec_ind(freq_len, struc_len)


  inds_avail = (/ 1, 2, 3, 4, 5, 6, 7, 8, &
       1, 2, 3, 4, 5, 6, 7, 8, &
       1, 2, 3, 4, 5, 6, 7, 8, &
       9, 10, 11, 12, 13, 14, 15, 16 /)

!!$  inds_avail = (/ 1, 2, 3, 4, 5, 6, 7, 8, &
!!$       9, 10, 11, 12, 13, 14, 15, 16 /)

!!$  inds_avail = (/ 1, 2, 3, 4, 5, 6, 7, 8, &
!!$       1, 2, 3, 4, 5, 6, 7, 8, &
!!$       1, 2, 3, 4, 5, 6, 7, 8, &
!!$       1, 2, 3, 4, 5, 6, 7, 8, &
!!$       1, 2, 3, 4, 5, 6, 7, 8, &
!!$       9, 10, 11, 12, 13, 14, 15, 16 /)

  sampled_opa_weights(:, 1, :, :) = 0d0
  sampled_opa_weights(:, 2, :, :) = 1d0
  k_min = 0d0
  k_max = 0d0
  weights_use = weights
  weights_use(1:8) = weights_use(1:8)/3d0
  take_spec = 0
  take_spec_ind = 1

!!$  weights_use(1:8) = weights_use(1:8)/5d0

  call init_random_seed()
  !call random_number(r_index)

  !time_test = TIME()
  !t1 = time_test

  ! Find threshold
  !time_test = TIME()
  !t0 = time_test

  ! In every layer and frequency bin:
  ! find the species with the largest kappa(g=0) value,
  ! save that value.
  do i_struc = 1, struc_len
     do i_freq = 1, freq_len
        threshold(i_freq, i_struc) = MAXVAL(line_struc_kappas(1, i_freq, :, i_struc))
     end do
  end do

  do i_struc = 1, struc_len
     do i_spec = 1, N_species
        do i_freq = 1, freq_len

           ! Only consider a species if kappa(g=1) > 0.01 * treshold
           if (line_struc_kappas(g_len, i_freq, i_spec, i_struc) < &
                threshold(i_freq, i_struc)*1d-2) then
              cycle
           end if

           take_spec(i_freq, i_struc) = take_spec(i_freq, i_struc)+1
           take_spec_ind(i_freq, i_struc) = i_spec

        end do
     end do
  end do

  !time_test = TIME()
  !t0 = time_test - t0

  do i_struc = 1, struc_len
     do i_spec = 1, N_species
        do i_freq = 1, freq_len

           ! Only do the sampling if more than one species is to be considered.
           if (take_spec(i_freq, i_struc) < 2) then
              cycle
           end if

           ! Check again: really sample the current species?
           if (line_struc_kappas(g_len, i_freq, i_spec, i_struc) < &
                threshold(i_freq, i_struc)*1d-2) then
              cycle
           end if

!!$           ind_use = inds_avail( &
!!$                int(r_index(:, i_freq, i_spec, i_struc)*(8*6))+1)

           call random_number(r_index)
           !ind_use = inds_avail( &
           !     int(r_index(:, i_freq, i_spec, i_struc)*(8*4))+1)
           ind_use = inds_avail(int(r_index*(8*4))+1)


!!$           ind_use = inds_avail( &
!!$                int(r_index(:, i_freq, i_spec, i_struc)*(8*2))+1)

           sampled_opa_weights(:, 1, i_freq, i_struc) = &
                sampled_opa_weights(:, 1, i_freq, i_struc) + &
                line_struc_kappas(ind_use, i_freq, i_spec, i_struc)

           sampled_opa_weights(:, 2, i_freq, i_struc) = &
                sampled_opa_weights(:, 2, i_freq, i_struc) * &
                weights_use(ind_use)

           k_min(i_freq, i_struc) = k_min(i_freq, i_struc) + &
                MINVAL(line_struc_kappas(:, i_freq, i_spec, i_struc))

           k_max(i_freq, i_struc) = k_max(i_freq, i_struc) + &
                MAXVAL(line_struc_kappas(:, i_freq, i_spec, i_struc))

        end do
     end do
     !write(*,*) take_spec(:, i_struc)
  end do

  !time_test = TIME()
  !t1 = time_test - t1

  !time_test = TIME()
  !t2 = time_test

  do i_struc = 1, struc_len
     do i_freq = 1, freq_len

        ! Interpolate new corr-k table if more than one species is to be considered
        if (take_spec(i_freq, i_struc) > 1) then

           call wrap_quicksort_swap(nsample, sampled_opa_weights(:, :, i_freq, i_struc))

           sampled_opa_weights(:, 2, i_freq, i_struc) = &
                sampled_opa_weights(:, 2, i_freq, i_struc) / &
                SUM(sampled_opa_weights(:, 2, i_freq, i_struc))

           g_sample = 0d0
           cum_sum = 0d0
           do i_samp = 1, nsample
              g_sample(i_samp) = &
                   sampled_opa_weights(i_samp, 2, i_freq, i_struc)/2d0 + &
                   cum_sum
              cum_sum = cum_sum + &
                   sampled_opa_weights(i_samp, 2, i_freq, i_struc)
           end do

           g_final(1) = 0d0
           g_final(2:nsample+1) = g_sample
           g_final(nsample+2) = 1d0

           k_final(1) = k_min(i_freq, i_struc)
           k_final(2:nsample+1) = sampled_opa_weights(:, 1, i_freq, i_struc)
           k_final(nsample+2) = k_max(i_freq, i_struc)

           call search_intp_ind(g_final, nsample+2, g_gauss, g_len, intpint)
           !if ((i_struc == 1) .AND. (i_freq == 1)) then
           do i_g = 1, g_len
!!$              write(*,*) g_final(intpint(i_g)), g_gauss(i_g), &
!!$                   g_final(intpint(i_g)+1), ((g_final(intpint(i_g)) <= &
!!$                   g_gauss(i_g)) .AND. &
!!$                   (g_gauss(i_g) <= g_final(intpint(i_g)+1)))

              line_struc_kappas_out(i_g, i_freq, i_struc) = &
                   k_final(intpint(i_g)) + &
                   (k_final(intpint(i_g)+1) - k_final(intpint(i_g))) / &
                   (g_final(intpint(i_g)+1) - g_final(intpint(i_g))) * &
                   (g_gauss(i_g) - g_final(intpint(i_g)))

           end do
           !end if

        ! Otherwise: just take the opacity of the only species as the full combined k-table
        else

           line_struc_kappas_out(:, i_freq, i_struc) = &
                line_struc_kappas(:, i_freq, take_spec_ind(i_freq, i_struc), i_struc)

        end if

     end do
  end do

  !time_test = TIME()
  !t2 = time_test - t2

  !write(*,*) 'Time spend where', t0/(t0+t1+t2), t1/(t0+t1+t2), t2/(t0+t1+t2)

end subroutine combine_opas_sample_ck

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Self-written? Too long ago... Check if not rather from numrep...
subroutine search_intp_ind(binbord,binbordlen,arr,arrlen,intpint)

  implicit none

  INTEGER            :: binbordlen, arrlen, intpint(arrlen)
  DOUBLE PRECISION   :: binbord(binbordlen),arr(arrlen)
  INTEGER            :: i_arr
  INTEGER            :: pivot, k0, km

  ! carry out a binary search for the interpolation bin borders
  do i_arr = 1, arrlen

     if (arr(i_arr) >= binbord(binbordlen)) then
        intpint(i_arr) = binbordlen - 1
     else if (arr(i_arr) <= binbord(1)) then
        intpint(i_arr) = 1
!!$        write(*,*) 'yes', arr(i_arr),binbord(1)
     else

        k0 = 1
        km = binbordlen
        pivot = (km+k0)/2

        do while(km-k0>1)

           if (arr(i_arr) >= binbord(pivot)) then
              k0 = pivot
              pivot = (km+k0)/2
           else
              km = pivot
              pivot = (km+k0)/2
           end if

        end do

        intpint(i_arr) = k0

     end if

  end do

end subroutine search_intp_ind

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine feautrier_rad_trans(border_freqs, &
     tau_approx_scat, &
     temp, &
     mu, &
     w_gauss_mu, &
     w_gauss_ck, &
     photon_destruct_in, &
     contribution, &
     surf_refl, &
     surf_emi, &
     I_star_0, &
     geom, &
     mu_star, &
     flux, &
     contr_em, &
     freq_len_p_1, &
     struc_len, &
     N_mu, &
     N_g)

  use constants_block
  implicit none

  ! I/O
  INTEGER, INTENT(IN)             :: freq_len_p_1, struc_len, N_mu, N_g
  DOUBLE PRECISION, INTENT(IN)    :: mu_star
  DOUBLE PRECISION, INTENT(IN)    :: surf_refl(freq_len_p_1-1),surf_emi(freq_len_p_1-1) !ELALEI
  DOUBLE PRECISION, INTENT(IN)    :: I_star_0(freq_len_p_1-1) !ELALEI
  DOUBLE PRECISION, INTENT(IN)    :: border_freqs(freq_len_p_1)
  DOUBLE PRECISION, INTENT(IN)    :: tau_approx_scat(N_g,freq_len_p_1-1,struc_len)
  DOUBLE PRECISION, INTENT(IN)    :: temp(struc_len)
  DOUBLE PRECISION, INTENT(IN)    :: mu(N_mu)
  DOUBLE PRECISION, INTENT(IN)    :: w_gauss_mu(N_mu), w_gauss_ck(N_g)
  DOUBLE PRECISION, INTENT(IN)    :: photon_destruct_in(N_g,freq_len_p_1-1,struc_len)
  LOGICAL, INTENT(IN)             :: contribution
  DOUBLE PRECISION, INTENT(OUT)   :: flux(freq_len_p_1-1)
  DOUBLE PRECISION, INTENT(OUT)   :: contr_em(struc_len,freq_len_p_1-1)
  CHARACTER*20, intent(in)        :: geom

  ! Internal
  INTEGER                         :: j,i,k,l
  DOUBLE PRECISION                :: I_J(struc_len,N_mu), I_H(struc_len,N_mu)
  DOUBLE PRECISION                :: source(N_g,freq_len_p_1-1,struc_len), &
       J_planet_scat(N_g,freq_len_p_1-1,struc_len), &
       photon_destruct(N_g,freq_len_p_1-1,struc_len), &
       source_planet_scat_n(N_g,freq_len_p_1-1,struc_len), &
       source_planet_scat_n1(N_g,freq_len_p_1-1,struc_len), &
       source_planet_scat_n2(N_g,freq_len_p_1-1,struc_len), &
       source_planet_scat_n3(N_g,freq_len_p_1-1,struc_len)
   DOUBLE PRECISION                :: J_star_ini(N_g,freq_len_p_1-1,struc_len)
   DOUBLE PRECISION                :: I_star_calc(N_g,N_mu,struc_len,freq_len_p_1-1)
   DOUBLE PRECISION                :: flux_old(freq_len_p_1-1), conv_val
  ! tridag variables
  DOUBLE PRECISION                :: a(struc_len),b(struc_len),c(struc_len),r(struc_len), &
       planck(struc_len)
  DOUBLE PRECISION                :: f1,f2,f3, deriv1, deriv2, I_plus, I_minus

  ! quantities for P-T structure iteration
  DOUBLE PRECISION                :: J_bol(struc_len)
  DOUBLE PRECISION                :: J_bol_a(struc_len)
  DOUBLE PRECISION                :: J_bol_g(struc_len)

  ! ALI
  DOUBLE PRECISION                :: lambda_loc(N_g,freq_len_p_1-1,struc_len)

  ! control
  DOUBLE PRECISION                :: inv_del_tau_min
  INTEGER                         :: iter_scat, i_iter_scat

  ! GCM spec calc
  LOGICAL                         :: GCM_read
  DOUBLE PRECISION                :: I_GCM(N_mu,freq_len_p_1-1)

  ! Variables for the contribution function calculation
  INTEGER :: i_mu, i_str, i_freq
  DOUBLE PRECISION :: transm_mu(N_g,freq_len_p_1-1,struc_len), &
                     transm_all(freq_len_p_1-1,struc_len), transm_all_loc(struc_len)

  ! PAUL NEW
  ! Variables for surface scattering
  DOUBLE PRECISION                :: I_plus_surface(N_mu, N_g, freq_len_p_1-1)


  I_plus_surface = 0d0
  I_minus = 0d0
  ! END PAUL NEW

  GCM_read = .FALSE.
  iter_scat = 100
  source = 0d0
  flux_old = 0d0
  flux = 0d0

  source_planet_scat_n = 0d0
  source_planet_scat_n1 = 0d0
  source_planet_scat_n2 = 0d0
  source_planet_scat_n3 = 0d0

  photon_destruct = photon_destruct_in

  ! DO THE STELLAR ATTENUATION CALCULATION

  J_star_ini = 0d0


  do i = 1, freq_len_p_1-1
    ! Irradiation treatment
    ! Dayside ave: multiply flux by 1/2.
    ! Planet ave: multiply flux by 1/4

    do i_mu = 1, N_mu
      if (trim(adjustl(geom)) .EQ. 'dayside_ave') then
           I_star_calc(:,i_mu,:,i) = 0.5* abs(I_star_0(i))*exp(-tau_approx_scat(:,i,:)/mu(i_mu))
           J_star_ini(:,i,:) = J_star_ini(:,i,:)+0.5d0*I_star_calc(:,i_mu,:,i)*w_gauss_mu(i_mu)
      else if (trim(adjustl(geom)) .EQ. 'planetary_ave') then
           I_star_calc(:,i_mu,:,i) = 0.25* abs(I_star_0(i))*exp(-tau_approx_scat(:,i,:)/mu(i_mu))
           J_star_ini(:,i,:) = J_star_ini(:,i,:)+0.5d0*I_star_calc(:,i_mu,:,i)*w_gauss_mu(i_mu)
      else if (trim(adjustl(geom)) .EQ. 'non-isotropic') then
           J_star_ini(:,i,:) = abs(I_star_0(i)/4.*exp(-tau_approx_scat(:,i,:)/mu_star))
      else
          write(*,*) 'Invalid geometry'
     end if
   end do
  end do


  do i_iter_scat = 1, iter_scat

    flux_old = flux

    lambda_loc = 0d0

    J_planet_scat = 0d0

    inv_del_tau_min = 1d10
    J_bol(1) = 0d0
    I_GCM = 0d0

    do i = 1, freq_len_p_1-1

       flux(i) = 0d0
       J_bol_a = 0d0

       r = 0

       call planck_f_lr(struc_len,temp(1:struc_len),border_freqs(i),border_freqs(i+1),r)
       planck = r

       do l = 1, N_g

          if (i_iter_scat .EQ. 1) then
             source(l,i,:) = photon_destruct(l,i,:)*r +  (1d0-photon_destruct(l,i,:))*J_star_ini(l,i,:)
          else
             r = source(l,i,:)

          end if


          do j = 1, N_mu



             ! Own boundary treatment
             f1 = mu(j)/(tau_approx_scat(l,i,1+1)-tau_approx_scat(l,i,1))

             ! own test against instability
             if (f1 > inv_del_tau_min) then
                f1 = inv_del_tau_min
             end if
             if (f1 .NE. f1) then
                f1 = inv_del_tau_min
             end if

             b(1) = 1d0 + 2d0 * f1 * (1d0 + f1)
             c(1) = -2d0*f1**2d0
             a(1) = 0d0

             ! Calculate the local approximate lambda iterator
             lambda_loc(l,i,1) = lambda_loc(l,i,1) + &
                  w_gauss_mu(j)/(1d0 + 2d0 * f1 * (1d0 + f1))

             do k = 1+1, struc_len-1

                f1 = 2d0*mu(j)/(tau_approx_scat(l,i,k+1)-tau_approx_scat(l,i,k-1))
                f2 = mu(j)/(tau_approx_scat(l,i,k+1)-tau_approx_scat(l,i,k))
                f3 = mu(j)/(tau_approx_scat(l,i,k)-tau_approx_scat(l,i,k-1))

                ! own test against instability
                if (f1 > 0.5d0*inv_del_tau_min) then
                   f1 = 0.5d0*inv_del_tau_min
                end if
                if (f1 .NE. f1) then
                   f1 = 0.5d0*inv_del_tau_min
                end if
                if (f2 > inv_del_tau_min) then
                   f2 = inv_del_tau_min
                end if
                if (f2 .NE. f2) then
                   f2 = inv_del_tau_min
                end if
                if (f3 > inv_del_tau_min) then
                   f3 = inv_del_tau_min
                end if
                if (f3 .NE. f3) then
                   f3 = inv_del_tau_min
                end if

                b(k) = 1d0 + f1*(f2+f3)
                c(k) = -f1*f2
                a(k) = -f1*f3

                ! Calculate the local approximate lambda iterator
                lambda_loc(l,i,k) = lambda_loc(l,i,k) + &
                     w_gauss_mu(j)/(1d0+f1*(f2+f3))

             end do

             ! Own boundary treatment
             f1 = mu(j)/(tau_approx_scat(l,i,struc_len)-tau_approx_scat(l,i,struc_len-1))

             ! own test against instability
             if (f1 > inv_del_tau_min) then
                f1 = inv_del_tau_min
             end if
             if (f1 .NE. f1) then
                f1 = inv_del_tau_min
             end if

  !!$              b(struc_len) = 1d0 + 2d0*f1**2d0
  !!$              c(struc_len) = 0d0
  !!$              a(struc_len) = -2d0*f1**2d0
  !!$
  !!$              ! Calculate the local approximate lambda iterator
  !!$              lambda_loc(l,i,struc_len) = lambda_loc(l,i,struc_len) + &
  !!$                   w_gauss_mu(j)/(1d0 + 2d0*f1**2d0)

             ! TEST PAUL SCAT
             b(struc_len) = 1d0
             c(struc_len) = 0d0
             a(struc_len) = 0d0

             ! r(struc_len) = I_J(struc_len) = 0.5[I_plus + I_minus]
             ! where I_plus is the light that goes downwards and
             ! I_minus is the light that goes upwards.
             !!!!!!!!!!!!!!!!!! ALWAYS NEEDED !!!!!!!!!!!!!!!!!!
             I_plus = I_plus_surface(j, l, i)

                            !!!!!!!!!!!!!!! EMISSION ONLY TERM !!!!!!!!!!!!!!!!
             I_minus = surf_emi(i)*planck(struc_len) &
                           !!!!!!!!!!!!!!! SURFACE SCATTERING !!!!!!!!!!!!!!!!
                           ! ----> of the emitted/scattered atmospheric light
                           ! + surf_refl(i) * SUM(I_plus_surface(:, l, i) * w_gauss_mu) ! OLD PRE 091220
                           + surf_refl(i) * 2d0 * SUM(I_plus_surface(:, l, i) * mu * w_gauss_mu)
                           ! ----> of the direct stellar beam (depends on geometry)
             if  (trim(adjustl(geom)) .NE. 'non-isotropic') then
               I_minus = I_minus + surf_refl(i) &
                    ! * SUM(I_star_calc(l,:, struc_len, i) * w_gauss_mu) ! OLD PRE 091220
                    * 2d0 * SUM(I_star_calc(l,:, struc_len, i) * mu * w_gauss_mu)
             else
               !I_minus = I_minus + surf_refl(i) *J_star_ini(l,i,struc_len)  !to be checked! ! OLD PRE 091220
               I_minus = I_minus + surf_refl(i) *J_star_ini(l,i,struc_len) * 4d0 * mu_star
             end if

             !sum to get I_J
             r(struc_len)=0.5*(I_plus + I_minus)

             ! Calculate the local approximate lambda iterator
             lambda_loc(l,i,struc_len) = lambda_loc(l,i,struc_len) + &
                  w_gauss_mu(j)/(1d0 + 2d0*f1**2d0)

              call tridag_own(a,b,c,r,I_J(:,j),struc_len)

             I_H(1,j) = -I_J(1,j)

             do k = 1+1, struc_len-1
                f1 = mu(j)/(tau_approx_scat(l,i,k+1)-tau_approx_scat(l,i,k))
                f2 = mu(j)/(tau_approx_scat(l,i,k)-tau_approx_scat(l,i,k-1))
                if (f1 > inv_del_tau_min) then
                   f1 = inv_del_tau_min
                end if
                if (f2 > inv_del_tau_min) then
                   f2 = inv_del_tau_min
                end if
                deriv1 = f1*(I_J(k+1,j)-I_J(k,j))
                deriv2 = f2*(I_J(k,j)-I_J(k-1,j))
                I_H(k,j) = -(deriv1+deriv2)/2d0

                ! TEST PAUL SCAT
                if (k .EQ. struc_len - 1) then
                   I_plus_surface(j, l, i) = &
                        I_J(struc_len,j)  - deriv1
                end if
                ! END TEST PAUL SCAT
             end do

             I_H(struc_len,j) = 0d0

             ! TEST PAUL SCAT
             !I_plus_surface(j, l, i) = I_J(struc_len-1,j)+I_H(struc_len-1,j)
             ! END TEST PAUL SCAT

          end do

          J_bol_g = 0d0

          do j = 1, N_mu

             J_bol_g = J_bol_g + I_J(:,j) * w_gauss_mu(j)
             flux(i) = flux(i) - I_H(1,j)*mu(j) &
                  * 4d0*pi * w_gauss_ck(l) * w_gauss_mu(j)
          end do

          ! Save angle-dependent surface flux
          if (GCM_read) then
             do j = 1, N_mu
                I_GCM(j,i) = I_GCM(j,i) - 2d0*I_H(1,j)*w_gauss_ck(l)
             end do
          end if

          J_planet_scat(l,i,:) = J_bol_g

       end do

    end do

    do k = 1, struc_len
       do i = 1, freq_len_p_1-1
          do l = 1, N_g
             if (photon_destruct(l,i,k) < 1d-10) THEN
                photon_destruct(l,i,k) = 1d-10
             end if
          end do
       end do
    end do

    do i = 1, freq_len_p_1-1
       call planck_f_lr(struc_len,temp(1:struc_len),border_freqs(i),border_freqs(i+1),r)
       do l = 1, N_g
         source(l,i,:) = (photon_destruct(l,i,:)*r+(1d0-photon_destruct(l,i,:))* &
               (J_star_ini(l,i,:)+J_planet_scat(l,i,:)-lambda_loc(l,i,:)*source(l,i,:))) / &
               (1d0-(1d0-photon_destruct(l,i,:))*lambda_loc(l,i,:))
       end do
    end do

    source_planet_scat_n3 = source_planet_scat_n2
    source_planet_scat_n2 = source_planet_scat_n1
    source_planet_scat_n1 = source_planet_scat_n
    source_planet_scat_n  = source

    if (mod(i_iter_scat,4) .EQ. 0) then
       !write(*,*) 'Ng acceleration!'
       call NG_source_approx(source_planet_scat_n,source_planet_scat_n1, &
            source_planet_scat_n2,source_planet_scat_n3,source, &
            N_g,freq_len_p_1,struc_len)
    end if

    conv_val = MAXVAL(ABS((flux-flux_old)/flux))
    if ((conv_val < 1d-2) .AND. (i_iter_scat > 9)) then
        exit
    end if
    
 end do

  ! Calculate the contribution function.
  ! Copied from flux_ck, here using "source" as the source function
  ! (before it was the Planck function).

  contr_em = 0d0
  if (contribution) then

    do i_mu = 1, N_mu

       ! Transmissions for a given incidence angle
       transm_mu = exp(-tau_approx_scat/mu(i_mu))

       do i_str = 1, struc_len
          do i_freq = 1, freq_len_p_1-1
             ! Integrate transmission over g-space
             transm_all(i_freq,i_str) = sum(transm_mu(:,i_freq,i_str)*w_gauss_ck)
          end do
       end do

       ! Do the actual radiative transport
       do i_freq = 1, freq_len_p_1-1
          ! Spatial transmissions at given wavelength
          transm_all_loc = transm_all(i_freq,:)
          ! Calc Eq. 9 of manuscript (em_deriv.pdf)
          do i_str = 1, struc_len
             r(i_str) = sum(source(:,i_freq,i_str)*w_gauss_ck)
          end do
          do i_str = 1, struc_len-1
             contr_em(i_str,i_freq) = contr_em(i_str,i_freq)+ &
                  (r(i_str)+r(i_str+1)) * &
                  (transm_all_loc(i_str)-transm_all_loc(i_str+1)) &
                  *mu(i_mu)*w_gauss_mu(i_mu)
          end do
          contr_em(struc_len,i_freq) = contr_em(struc_len,i_freq)+ &
              !2d0*r(struc_len)*transm_all_loc(struc_len)*mu(i_mu)*w_gauss_mu(i_mu) MODIFIED TO INCLUDE SURFACE
              2d0*I_minus*transm_all_loc(struc_len)*mu(i_mu)*w_gauss_mu(i_mu)
       end do

    end do

    do i_freq = 1, freq_len_p_1-1
       contr_em(:,i_freq) = contr_em(:,i_freq)/SUM(contr_em(:,i_freq))
    end do

  end if

end subroutine feautrier_rad_trans


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine NG_source_approx(source_n,source_n1,source_n2,source_n3,source, &
     N_g,freq_len_p_1,struc_len)

  implicit none
  INTEGER :: struc_len, freq_len_p_1, N_g, i, i_ng, i_freq
  DOUBLE PRECISION :: tn(struc_len), tn1(struc_len), tn2(struc_len), &
       tn3(struc_len), temp_buff(struc_len), &
       source_n(N_g,freq_len_p_1-1,struc_len), source_n1(N_g,freq_len_p_1-1,struc_len), &
       source_n2(N_g,freq_len_p_1-1,struc_len), source_n3(N_g,freq_len_p_1-1,struc_len), &
       source(N_g,freq_len_p_1-1,struc_len), source_buff(N_g,freq_len_p_1-1,struc_len)
  DOUBLE PRECISION :: Q1(struc_len), Q2(struc_len), Q3(struc_len)
  DOUBLE PRECISION :: A1, A2, B1, B2, C1, C2
  DOUBLE PRECISION :: a, b

  do i_freq = 1, freq_len_p_1-1
     do i_ng = 1, N_g

        tn = source_n(i_ng,i_freq,1:struc_len)
        tn1 = source_n1(i_ng,i_freq,1:struc_len)
        tn2 = source_n2(i_ng,i_freq,1:struc_len)
        tn3 = source_n3(i_ng,i_freq,1:struc_len)

        Q1 = tn - 2d0*tn1 + tn2
        Q2 = tn - tn1 - tn2 + tn3
        Q3 = tn - tn1

        ! test
        Q1(1) = 0d0
        Q2(1) = 0d0
        Q3(1) = 0d0

        A1 = sum(Q1*Q1)
        A2 = sum(Q2*Q1)
        B1 = sum(Q1*Q2)
        B2 = sum(Q2*Q2)
        C1 = sum(Q1*Q3)
        C2 = sum(Q2*Q3)

        if ((abs(A1) >= 1d-250) .AND. &
             (abs(A2) >=1d-250) .AND. &
             (abs(B1) >=1d-250) .AND. &
             (abs(B2) >=1d-250) .AND. &
             (abs(C1) >=1d-250) .AND. &
             (abs(C2) >=1d-250)) THEN

           a = (C1*B2-C2*B1)/(A1*B2-A2*B1)
           b = (C2*A1-C1*A2)/(A1*B2-A2*B1)

           temp_buff = (1d0-a-b)*tn + a*tn1 + b*tn2

           do i = 1,struc_len
              if (temp_buff(i) <= 0d0) then
                 temp_buff(i) = 0d0
              end if
           end do

           do i = 1,struc_len
              if (temp_buff(i) .NE. temp_buff(i)) then
                 return
              end if
           end do

           source_buff(i_ng,i_freq,1:struc_len) = temp_buff

        else

           source_buff(i_ng,i_freq,1:struc_len) = source(i_ng,i_freq,1:struc_len)

        end if

     end do
  end do

  source = source_buff

end subroutine NG_source_approx

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!**********************************************************
! RANDOM SEED GENERATOR BELOW TAKEN FROM
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!**********************************************************

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid, getpid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to OR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine wrap_quicksort_swap(length, array)

  implicit none
  integer, intent(in) :: length
  double precision, intent(inout) :: array(length, 2)
  double precision :: swapped_array(2, length)

  swapped_array(1, :) = array(:, 1)
  swapped_array(2, :) = array(:, 2)
  call quicksort_own_2d_swapped(length, swapped_array)
  array(:,1) = swapped_array(1,:)
  array(:,2) = swapped_array(2,:)

end subroutine wrap_quicksort_swap

recursive subroutine quicksort_own_2d_swapped(length, array)

  implicit none
  integer, intent(in) :: length
  double precision, intent(inout) :: array(2,length)
  integer :: partition_index
  integer :: ind_up, &
       ind_down, &
       ind_down_start
  double precision :: buffer(2), compare(2)
  logical :: found

  found = .False.

  partition_index = length
  compare = array(:, partition_index)

  ind_down_start = length-1

  do ind_up = 1, length-1

     if (array(1,ind_up) > compare(1)) then

        found = .True.

        do ind_down = ind_down_start, 1, -1

           if (ind_down == ind_up) then

              array(:,partition_index) = array(:,ind_down)
              array(:,ind_down) = compare

              if ((length-ind_down) > 1) then
                 call quicksort_own_2d_swapped(length-ind_down, array(:,ind_down+1:length))
              end if
              if ((ind_down-1) > 1) then
                 call quicksort_own_2d_swapped(ind_down-1, array(:,1:ind_down-1))
              end if
              return

           else if (array(1,ind_down) < compare(1)) then

              buffer = array(:,ind_up)
              array(:,ind_up) = array(:,ind_down)
              array(:,ind_down) = buffer
              ind_down_start = ind_down
              exit

           end if

        end do

     end if

  end do

  if (found .EQV. .FALSE.) then

     if ((length-1) > 1 ) then
        call quicksort_own_2d_swapped(length-1, array(:,1:length-1))
     end if
  end if

end subroutine quicksort_own_2d_swapped

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Tridag, own implementation, following the numerical recipes book.

subroutine tridag_own(a,b,c,res,solution,length)

  implicit none

  ! I/O
  integer, intent(in) :: length
  double precision, intent(in) :: a(length), &
       b(length), &
       c(length), &
       res(length)
  double precision, intent(out) :: solution(length)

  ! Internal variables
  integer :: ind
  double precision :: buffer_scalar, &
       buffer_vector(length)

  ! Test if b(1) == 0:
  if (b(1) .EQ. 0) then
     stop "Error in tridag routine, b(1) must not be zero!"
     end if

  ! Begin inversion
  buffer_scalar = b(1)
  solution(1) = res(1) / buffer_scalar

  do ind = 2, length
     buffer_vector(ind) = c(ind-1)/buffer_scalar
     buffer_scalar = b(ind) - a(ind) * buffer_vector(ind)
     if (buffer_scalar .EQ. 0) then
        write(*,*) "Tridag routine failed!"
        solution = 0d0
  return
     end if
     solution(ind) = (res(ind) - &
          a(ind)*solution(ind-1))/buffer_scalar
  end do

  do ind = length-1, 1, -1
     solution(ind) = solution(ind) &
          - buffer_vector(ind+1) * solution(ind + 1)
  end do

end subroutine tridag_own


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine planck_f_lr(PT_length,T,nul,nur,B_nu)

  use constants_block
  implicit none
  INTEGER                         :: PT_length
  DOUBLE PRECISION                :: T(PT_length),B_nu(PT_length)
  DOUBLE PRECISION                ::  nu1, nu2, nu3, nu4, nu5, nu_large, nu_small, &
       nul, nur, diff_nu

  !~~~~~~~~~~~~~

  B_nu = 0d0
  ! Take mean using Boole's method
  nu_large = max(nul,nur)
  nu_small = min(nul,nur)
  nu1 = nu_small
  nu2 = nu_small+DBLE(1)*(nu_large-nu_small)/4d0
  nu3 = nu_small+DBLE(2)*(nu_large-nu_small)/4d0
  nu4 = nu_small+DBLE(3)*(nu_large-nu_small)/4d0
  nu5 = nu_large
  diff_nu = nu2-nu1
  B_nu = B_nu + 1d0/90d0*( &
       7d0* 2d0*hplanck*nu1**3d0/c_l**2d0/(exp(hplanck*nu1/kB/T)-1d0) + &
       32d0*2d0*hplanck*nu2**3d0/c_l**2d0/(exp(hplanck*nu2/kB/T)-1d0) + &
       12d0*2d0*hplanck*nu3**3d0/c_l**2d0/(exp(hplanck*nu3/kB/T)-1d0) + &
       32d0*2d0*hplanck*nu4**3d0/c_l**2d0/(exp(hplanck*nu4/kB/T)-1d0) + &
       7d0* 2d0*hplanck*nu5**3d0/c_l**2d0/(exp(hplanck*nu5/kB/T)-1d0))

end subroutine planck_f_lr

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine calc_rosse_opa(HIT_kappa_tot_g_approx,HIT_border_freqs,temp,HIT_N_g,HIT_coarse_borders, &
     kappa_rosse, w_gauss)

  use constants_block
  implicit none
  INTEGER                         :: HIT_N_g,HIT_coarse_borders
  DOUBLE PRECISION                :: HIT_border_freqs(HIT_coarse_borders)
  DOUBLE PRECISION                :: HIT_kappa_tot_g_approx(HIT_N_g,HIT_coarse_borders-1)
  DOUBLE PRECISION                :: temp, kappa_rosse, w_gauss(HIT_N_g), B_nu_dT(HIT_coarse_borders-1), &
       numerator
  INTEGER                         :: i

  !~~~~~~~~~~~~~

  call star_planck_div_T(HIT_coarse_borders,temp,HIT_border_freqs,B_nu_dT)

  kappa_rosse = 0d0
  numerator = 0d0

  do i = 1, HIT_coarse_borders-1
     kappa_rosse = kappa_rosse + &
          B_nu_dT(i) * sum(w_gauss/HIT_kappa_tot_g_approx(:,i)) * &
          (HIT_border_freqs(i)-HIT_border_freqs(i+1))
     numerator = numerator + &
          B_nu_dT(i) * &
          (HIT_border_freqs(i)-HIT_border_freqs(i+1))
  end do

  kappa_rosse = numerator / kappa_rosse

end subroutine calc_rosse_opa

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine star_planck_div_T(freq_len,T,nu,B_nu_dT)

  use constants_block
  implicit none
  INTEGER                         :: freq_len
  DOUBLE PRECISION                :: T,B_nu_dT(freq_len-1),nu(freq_len)
  DOUBLE PRECISION                :: buffer(freq_len-1),nu_use(freq_len-1)
  INTEGER                         :: i

  !~~~~~~~~~~~~~

  do i = 1, freq_len-1
     nu_use(i) = (nu(i)+nu(i+1))/2d0
  end do

  buffer = 2d0*hplanck**2d0*nu_use**4d0/c_l**2d0
  B_nu_dT = buffer / ((exp(hplanck*nu_use/kB/T/2d0)-exp(-hplanck*nu_use/kB/T/2d0))**2d0)/kB/T**2d0

end subroutine star_planck_div_T

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine calc_planck_opa(HIT_kappa_tot_g_approx,HIT_border_freqs,temp,HIT_N_g,HIT_coarse_borders, &
     kappa_planck, w_gauss)

  use constants_block
  implicit none
  INTEGER                         :: HIT_N_g,HIT_coarse_borders
  DOUBLE PRECISION                :: HIT_border_freqs(HIT_coarse_borders)
  DOUBLE PRECISION                :: HIT_kappa_tot_g_approx(HIT_N_g,HIT_coarse_borders-1)
  DOUBLE PRECISION                :: temp, kappa_planck, w_gauss(HIT_N_g), B_nu(HIT_coarse_borders-1), &
       norm

  INTEGER                         :: i

  call star_planck(HIT_coarse_borders,temp,HIT_border_freqs,B_nu)

  kappa_planck = 0d0
  norm = 0d0
  do i = 1, HIT_coarse_borders-1
     kappa_planck = kappa_planck + &
          B_nu(i) * sum(HIT_kappa_tot_g_approx(:,i)*w_gauss) * &
          (HIT_border_freqs(i)-HIT_border_freqs(i+1))
     norm = norm + &
          B_nu(i) * (HIT_border_freqs(i)-HIT_border_freqs(i+1))
  end do

  kappa_planck = kappa_planck / norm

end subroutine calc_planck_opa

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine star_planck(freq_len,T,nu,B_nu)

  use constants_block
  implicit none
  INTEGER                         :: freq_len
  DOUBLE PRECISION                :: T,B_nu(freq_len-1), B_nu_l(freq_len),nu(freq_len)
  DOUBLE PRECISION                :: buffer(freq_len-1), buffer_l(freq_len),nu_use(freq_len-1)
  INTEGER                         :: i, integ
  DOUBLE PRECISION                :: diff_nu, nu1, nu2, nu3, nu4, nu5, nu_large, nu_small, diff_nu_sum
  DOUBLE PRECISION                :: norm_B1, norm_B2, norm_B3, norm_corr

  !~~~~~~~~~~~~~

  B_nu = 0d0

  ! Take mean using Boole's method
  do i = 1, freq_len-1
     nu_large = max(nu(i),nu(i+1))
     nu_small = min(nu(i),nu(i+1))
     nu1 = nu_small
     nu2 = nu_small+DBLE(1)*(nu_large-nu_small)/4d0
     nu3 = nu_small+DBLE(2)*(nu_large-nu_small)/4d0
     nu4 = nu_small+DBLE(3)*(nu_large-nu_small)/4d0
     nu5 = nu_large
     diff_nu = nu2-nu1
     B_nu(i) = B_nu(i) + 1d0/90d0*( &
             7d0* 2d0*hplanck*nu1**3d0/c_l**2d0/(exp(hplanck*nu1/kB/T)-1d0) + &
             32d0*2d0*hplanck*nu2**3d0/c_l**2d0/(exp(hplanck*nu2/kB/T)-1d0) + &
             12d0*2d0*hplanck*nu3**3d0/c_l**2d0/(exp(hplanck*nu3/kB/T)-1d0) + &
             32d0*2d0*hplanck*nu4**3d0/c_l**2d0/(exp(hplanck*nu4/kB/T)-1d0) + &
             7d0* 2d0*hplanck*nu5**3d0/c_l**2d0/(exp(hplanck*nu5/kB/T)-1d0))
  end do

end subroutine star_planck
