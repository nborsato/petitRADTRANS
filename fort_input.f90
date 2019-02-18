!!$*****************************************************************************
!!$*****************************************************************************
!!$*****************************************************************************
!!$ fort_input.f90: utility functions to read, interpolate and mix opacities
!!$                 for the petitRADTRANS radiative transfer package
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

!!$ Subroutine to get length of frequency grid in correlated-k mode

subroutine get_freq_len(path,freq_len,g_len)

  implicit none
  ! I/O
  character*250, intent(in) :: path
  integer, intent(out) :: freq_len, g_len

  g_len = 1
  open(unit=10,file=trim(adjustl(path))//'/opacities/lines/corr_k/H2O/kappa_g_info.dat')
  read(10,*) freq_len, g_len
  close(10)
  freq_len = freq_len-1
  ! g_len = g_len-2 ! for OLD 32 grid

end subroutine get_freq_len

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to read in frequency grid

subroutine get_freq(path,freq_len,freq)

  implicit none
  ! I/O
  character*150, intent(in) :: path
  integer, intent(in) :: freq_len
  double precision, intent(out) :: freq(freq_len)
  ! internal
  integer :: i_freq, freq_len_use_ck
  double precision :: buffer, freq_use_ck(freq_len+1)

  ! Because freqs fot c-k are stored as borders!
  freq_len_use_ck = freq_len + 1
  open(unit=10,file=trim(adjustl(path))//'/opacities/lines/corr_k/H2O/kappa_g_info.dat')
  read(10,*)
  do i_freq = 1, freq_len_use_ck-2
     read(10,*) buffer, freq_use_ck(i_freq)
  end do
  read(10,*) freq_use_ck(freq_len_use_ck), &
       freq_use_ck(freq_len_use_ck-1)
  close(10)
  freq = (freq_use_ck(1:freq_len_use_ck-1)+freq_use_ck(2:freq_len_use_ck))/2d0

end subroutine get_freq

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to read in the molecular opacities (c-k or line-by-line)

subroutine read_in_molecular_opacities(path,species_names_tot,freq_len,g_len,species_len,opa_TP_grid_len, &
     opa_grid_kappas, mode, arr_min, arr_max)

  implicit none
  ! I/O
  character*150, intent(in) :: path
  character*5000, intent(in) :: species_names_tot
  integer, intent(in) :: freq_len,g_len,species_len, opa_TP_grid_len, arr_min, arr_max
  character*3, intent(in) :: mode
  double precision, intent(out) :: opa_grid_kappas(g_len,freq_len,species_len,opa_TP_grid_len)
  ! Internal
  character*2 :: species_id
  character*150 :: path_names(opa_TP_grid_len)
  character*400 :: path_read_stream
  !character*150 :: species_names(species_len)
  integer :: species_name_inds(2,species_len)
  double precision :: molparam, read_val, buffer
  integer :: i_spec, i_file, i_str, curr_spec_ind, &
       i_kg, curr_N_g_int, curr_cb_int

  ! Get single species names
  curr_spec_ind = 1
  species_name_inds(1,curr_spec_ind) = 1
  do i_str = 1, 5000
     if (curr_spec_ind > species_len) then
        EXIT
     end if
     if (species_names_tot(i_str:i_str) .EQ. ':') then
        species_name_inds(2,curr_spec_ind) = i_str-1
        curr_spec_ind = curr_spec_ind+1
        if (curr_spec_ind <= species_len) then
           species_name_inds(1,curr_spec_ind) = i_str+1
        end if
     end if
  end do

  ! Get paths of opacity files
  open(unit=20,file=trim(adjustl(path))//'/opa_input_files/opa_filenames.txt')
  do i_file = 1, opa_TP_grid_len
     read(20,*) path_names(i_file)
  end do
  close(20)

  write(*,*)
  ! Read opas for every species...
  do i_spec = 1, species_len
     ! Get species file ID and molparam
     if (mode .EQ. 'c-k') then
        open(unit=20,file=trim(adjustl(path))//'/opacities/lines/corr_k/' &
             //trim(adjustl(species_names_tot(species_name_inds(1,i_spec): &
             species_name_inds(2,i_spec))))//'/molparam_id.txt')
     else if (mode .EQ. 'lbl') then
        open(unit=20,file=trim(adjustl(path))//'/opacities/lines/line_by_line/' &
             //trim(adjustl(species_names_tot(species_name_inds(1,i_spec): &
             species_name_inds(2,i_spec))))//'/molparam_id.txt')
     end if
     write(*,*) ' Read line opacities of '//trim(adjustl(species_names_tot(species_name_inds(1, &
          i_spec):species_name_inds(2,i_spec))))//'...'
     read(20,*)
     read(20,'(A2)') species_id
     read(20,*)
     read(20,*) molparam
     close(20)
     ! ...for every P-T grid point...
     do i_file = 1, opa_TP_grid_len
        ! Open opacity file
        if (mode .EQ. 'c-k') then

           open(unit=20,file=trim(adjustl(path))//'/opacities/lines/corr_k/' &
                //trim(adjustl(species_names_tot(species_name_inds(1,i_spec): &
                species_name_inds(2,i_spec))))//'/sigma_'//species_id// &
                adjustl(trim(path_names(i_file))), form='unformatted')
           
        else if (mode .EQ. 'lbl') then

           path_read_stream =  trim(adjustl(path))//'/opacities/lines/line_by_line/' &
                //trim(adjustl(species_names_tot(species_name_inds(1,i_spec): &
                species_name_inds(2,i_spec))))//'/sigma_'//species_id// &
                adjustl(trim(path_names(i_file)))
           call read_kappa(arr_min, arr_max, freq_len, &
                path_read_stream, opa_grid_kappas(1,:,i_spec,i_file))
           
        end if
        ! ...for every frequency point.
        if (mode .EQ. 'c-k') then
           do i_kg = 1, g_len*freq_len
              curr_cb_int = (i_kg-1)/g_len+1
              curr_N_g_int = i_kg - (curr_cb_int-1)*g_len
              read(20) opa_grid_kappas(curr_N_g_int,curr_cb_int,i_spec,i_file)
           end do
           close(20)
        end if
        
     end do
     opa_grid_kappas(:,:,i_spec,:) = opa_grid_kappas(:,:,i_spec,:)/molparam
  end do

  write(*,*) 'Done.'
  write(*,*)
  
end subroutine read_in_molecular_opacities

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to read in the molecular opacities (c-k or line-by-line)

subroutine read_in_cloud_opacities(path,species_names_tot,species_modes_tot,N_cloud_spec, &
     N_cloud_lambda_bins,rho_cloud_particles,cloud_specs_abs_opa,cloud_specs_scat_opa, &
     cloud_aniso,cloud_lambdas,cloud_rad_bins,cloud_radii)

  implicit none
  ! Params
  integer, parameter :: N_cloud_rad_bins = 130
  
  ! I/O
  character*150, intent(in) :: path
  character*5000, intent(in) :: species_names_tot,species_modes_tot
  integer, intent(in) :: N_cloud_spec,N_cloud_lambda_bins

  double precision, intent(out) :: rho_cloud_particles(N_cloud_spec)
  double precision, intent(out) :: cloud_specs_abs_opa(N_cloud_rad_bins,N_cloud_lambda_bins,N_cloud_spec), &
       cloud_specs_scat_opa(N_cloud_rad_bins,N_cloud_lambda_bins,N_cloud_spec), &
       cloud_aniso(N_cloud_rad_bins,N_cloud_lambda_bins,N_cloud_spec), cloud_lambdas(N_cloud_lambda_bins), &
       cloud_rad_bins(N_cloud_rad_bins+1), cloud_radii(N_cloud_rad_bins)

  ! Internal
  integer :: i_str, curr_spec_ind, i_cloud, i_cloud_read, i_cloud_lamb, i_size, i_lamb
  integer :: species_name_inds(2,N_cloud_spec), species_mode_inds(2,N_cloud_spec)
  character*80  :: cloud_opa_names(N_cloud_spec), cloud_name_buff, buff_line, path_add
  double precision :: cloud_dens_buff, buffer
  character*2 :: cloud_opa_mode(N_cloud_spec)

  ! Get single cloud species names
  curr_spec_ind = 1
  species_name_inds(1,curr_spec_ind) = 1
  do i_str = 1, 5000
     if (curr_spec_ind > N_cloud_spec) then
        EXIT
     end if
     if (species_names_tot(i_str:i_str) .EQ. ':') then
        species_name_inds(2,curr_spec_ind) = i_str-1
        curr_spec_ind = curr_spec_ind+1
        if (curr_spec_ind <= N_cloud_spec) then
           species_name_inds(1,curr_spec_ind) = i_str+1
        end if
     end if
  end do

  ! Get single cloud species modes
  curr_spec_ind = 1
  species_mode_inds(1,curr_spec_ind) = 1
  do i_str = 1, 5000
     if (curr_spec_ind > N_cloud_spec) then
        EXIT
     end if
     if (species_modes_tot(i_str:i_str) .EQ. ':') then
        species_mode_inds(2,curr_spec_ind) = i_str-1
        curr_spec_ind = curr_spec_ind+1
        if (curr_spec_ind <= N_cloud_spec) then
           species_mode_inds(1,curr_spec_ind) = i_str+1
        end if
     end if
  end do

  ! Read in cloud densities
  rho_cloud_particles = -1d0
  DO i_cloud = 1, N_cloud_spec
     
     cloud_opa_names(i_cloud) = species_names_tot(species_name_inds(1,i_cloud): &
          species_name_inds(2,i_cloud))
     
     open(unit=10,file=trim(adjustl(path))//'/opa_input_files/cloud_names.dat')
     open(unit=11,file=trim(adjustl(path))//'/opa_input_files/cloud_densities.dat')
     do i_cloud_read = 1, 1000000
        read(10,*,end=199) cloud_name_buff
        read(11,*) cloud_dens_buff
        if (trim(adjustl(cloud_name_buff)) .EQ. &
             trim(adjustl(cloud_opa_names(i_cloud)))) then
           rho_cloud_particles(i_cloud) = cloud_dens_buff
        end if
     end do
199  close(10)
     close(11)
     IF (rho_cloud_particles(i_cloud) < 0d0) THEN
        WRITE(*,*) 'ERROR! DENSITY FOR CLOUD SPECIES '//trim( &
             adjustl(cloud_opa_names(i_cloud))) &
             //'NOT FOUND!'
        STOP
     END IF
  END DO
  
  ! Read in cloud opacities
  cloud_specs_abs_opa = 0d0
  cloud_specs_scat_opa = 0d0
  cloud_aniso = 0d0

  open(unit=10,file=trim(adjustl(path))// &
       '/opacities/continuum//clouds/MgSiO3_c/amorphous/mie/bin_borders.dat')
  read(10,*)
  do i_cloud_lamb = 1, N_cloud_rad_bins
     read(10,*) cloud_rad_bins(i_cloud_lamb)
  end do
  read(10,*) cloud_rad_bins(N_cloud_rad_bins+1)
  close(10)
  
  open(unit=11,file=trim(adjustl(path))// &
       '/opacities/continuum//clouds/MgSiO3_c/amorphous/mie/particle_sizes.dat')
  read(11,*)
  do i_cloud_lamb = 1, N_cloud_rad_bins
     read(11,'(A80)') buff_line
     read(buff_line(17:len(buff_line)),*) cloud_radii(i_cloud_lamb)
  end do
  close(11)

  open(unit=10,file=trim(adjustl(path))// &
       '/opacities/continuum//clouds/MgSiO3_c/amorphous/mie/opa_0001.dat')
  do i_cloud_lamb = 1,11
     read(10,*)
  end do
  do i_cloud_lamb = 1, N_cloud_lambda_bins
     read(10,*) cloud_lambdas(i_cloud_lamb)
     cloud_lambdas(i_cloud_lamb) = cloud_lambdas(i_cloud_lamb) / 1d4
  end do
  close(10)

  DO i_cloud = 1, N_cloud_spec

     cloud_opa_mode(i_cloud) = species_modes_tot(species_mode_inds(1,i_cloud): &
          species_mode_inds(2,i_cloud))
     
     path_add = trim(adjustl( &
          cloud_opa_names(i_cloud)(1:len(trim(adjustl( &
          cloud_opa_names(i_cloud))))-3)))

     if (trim(adjustl( &
          cloud_opa_names(i_cloud)(len(trim(adjustl( &
          cloud_opa_names(i_cloud))))-2: &
          len(trim(adjustl( &
          cloud_opa_names(i_cloud))))))) .EQ. '(c)') then
        path_add = trim(adjustl(path_add))//'_c'
     else if (trim(adjustl( &
          cloud_opa_names(i_cloud)(len(trim(adjustl( &
          cloud_opa_names(i_cloud))))-2: &
          len(trim(adjustl( &
          cloud_opa_names(i_cloud))))))) .EQ. '(L)') then
        path_add = trim(adjustl(path_add))//'_L'
     end if

     write(*,*) ' Read in opacity of cloud species ' &
          //trim(adjustl(path_add(1:len(trim(adjustl(path_add)))-2)))//' ...'

     if (cloud_opa_mode(i_cloud)(1:1) .EQ. 'a') then
        path_add = trim(adjustl(path_add))//'/amorphous'
     else if (cloud_opa_mode(i_cloud)(1:1) .EQ. 'c') then
        path_add = trim(adjustl(path_add))//'/crystalline'
     end if

     if (cloud_opa_mode(i_cloud)(2:2) .EQ. 'm') then
        path_add = trim(adjustl(path_add))//'/mie'
     else if (cloud_opa_mode(i_cloud)(2:2) .EQ. 'd') then
        path_add = trim(adjustl(path_add))//'/DHS'
        ! Decrease cloud particle density due to porosity
        rho_cloud_particles(i_cloud) = rho_cloud_particles(i_cloud)*0.75d0
     end if

     open(unit=11,file=trim(adjustl(path))// &
          '/opacities/continuum//clouds/'//trim(adjustl(path_add))// &
          '/particle_sizes.dat')
     read(11,*)
     do i_size = 1, N_cloud_rad_bins
        read(11,'(A80)') buff_line
        open(unit=10,file=trim(adjustl(path))// &
             '/opacities/continuum//clouds/'//trim(adjustl(path_add))// &
             '/'//trim(adjustl(buff_line(1:17))))
        do i_lamb = 1,11
           read(10,*)
        end do
        do i_lamb = 1, N_cloud_lambda_bins
           read(10,*) buffer, cloud_specs_abs_opa(i_size,i_lamb,i_cloud), &
                cloud_specs_scat_opa(i_size,i_lamb,i_cloud), &
                cloud_aniso(i_size,i_lamb,i_cloud)
        end do
        close(10)
     end do
     close(11)
  END DO
  write(*,*) 'Done.'
  write(*,*)


end subroutine read_in_cloud_opacities

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to interpolate the total opacity at a given PT structure

subroutine interpol_opa_ck(press,temp,struc_len,opa_TP_grid, &
     opa_grid_kappas,tp2nddim,N_PT_grid,N_species,freq_len,g_len,opa_struc_kappas)

  implicit none
  ! I/O
  INTEGER, intent(in)                   :: struc_len, N_PT_grid, g_len
  INTEGER, intent(in)                   :: N_species, freq_len, tp2nddim
  DOUBLE PRECISION, intent(in)          :: press(struc_len), temp(struc_len)
  DOUBLE PRECISION, intent(in)          :: opa_TP_grid(N_PT_grid,tp2nddim)
  DOUBLE PRECISION, intent(in)          :: opa_grid_kappas(g_len,freq_len,N_species,N_PT_grid)
  DOUBLE PRECISION, intent(out)         :: opa_struc_kappas(g_len,freq_len,N_species,struc_len)
  ! internal
  INTEGER                               :: i_str, s_temp_ind
  INTEGER                               :: press_s_t_s_ind, press_s_t_l_ind, &
       s_temp_ind_own, press_ind_own
  INTEGER                               :: PT_ind_Ts_Ps, PT_ind_Ts_Pl,PT_ind_Tl_Ps, &
       PT_ind_Tl_Pl
  DOUBLE PRECISION                      :: PorT(g_len,freq_len,N_species), &
       slopes(g_len,freq_len,N_species), &
       buffer1(g_len,freq_len,N_species),buffer2(g_len,freq_len,N_species)
  DOUBLE PRECISION                      :: buffer_Ts(g_len,freq_len,N_species), &
       buffer_Tl(g_len,freq_len,N_species),temp_min, temp_max

  !~~~~~~~~~~~~~

  press_s_t_l_ind = 0
  press_s_t_s_ind = 0
  s_temp_ind = 0

  temp_min = MINVAL(opa_TP_grid(:,1))
  temp_max = MAXVAL(opa_TP_grid(:,1))

  do i_str = 1, struc_len

     s_temp_ind_own = MAX(MIN(INT(log10(temp(i_str)/81.14113604736988d0)/ &
          log10(2995d0/81.14113604736988d0)*12d0)+1,12),1)
     press_ind_own = MAX(MIN(INT(log10(press(i_str)*1d-6)+6d0)+1,9),1)

     ! Opacity N_PT_grid indice at smaller P and T than point of interest
     PT_ind_Ts_Ps = (s_temp_ind_own-1)*10+press_ind_own
     ! Opacity N_PT_grid indice at larger P and smaller T than point of interest
     PT_ind_Ts_Pl = (s_temp_ind_own-1)*10+press_ind_own+1
     ! Opacity N_PT_grid indice at smaller P and larger T than point of interest
     PT_ind_Tl_Ps = s_temp_ind_own*10+press_ind_own
     ! Opacity N_PT_grid indice at larger P and T than point of interest
     PT_ind_Tl_Pl = s_temp_ind_own*10+press_ind_own+1

     ! Interpolate...

     !**********************************************************
     ! Interpolation to correct pressure at smaller temperatures
     !**********************************************************

     ! kappas

     ! kappa_gs at smaller T and smaller P
     buffer1 = opa_grid_kappas(:,:,:,PT_ind_Ts_Ps)
     ! kappa_gs at smaller T and larger P
     buffer2 = opa_grid_kappas(:,:,:,PT_ind_Ts_Pl)

     PorT = press(i_str)-opa_TP_grid(PT_ind_Ts_Ps,2)

     slopes = (buffer2-buffer1)/(opa_TP_grid(PT_ind_Ts_Pl,2)-opa_TP_grid(PT_ind_Ts_Ps,2))

     if (press(i_str) >= opa_TP_grid(PT_ind_Ts_Pl,2)) then
        buffer_Ts = buffer2
     else if (press(i_str) <= opa_TP_grid(PT_ind_Ts_Ps,2)) then
        buffer_Ts = buffer1
     else
        buffer_Ts = buffer1 + slopes*PorT
     end if

     !*********************************************************
     ! Interpolation to correct pressure at larger temperatures
     !*********************************************************

     ! kappas

     ! opacity at larger T and smaller P
     buffer1 = opa_grid_kappas(:,:,:,PT_ind_Tl_Ps)
     ! opacity at larger T and larger P
     buffer2 = opa_grid_kappas(:,:,:,PT_ind_Tl_Pl)

     PorT = press(i_str)-opa_TP_grid(PT_ind_Tl_Ps,2)

     ! slopes to correct to correct pressure are larger T
     slopes = (buffer2-buffer1)/(opa_TP_grid(PT_ind_Tl_Pl,2)-opa_TP_grid(PT_ind_Tl_Ps,2))

     ! total opacity at larger temperature and correct pressure
     if (press(i_str) >= opa_TP_grid(PT_ind_Tl_Pl,2)) then
        buffer_Tl = buffer2
     else if (press(i_str) <= opa_TP_grid(PT_ind_Tl_Ps,2)) then
        buffer_Tl = buffer1
     else
        buffer_Tl = buffer1 + slopes*PorT
     end if

     !***********************************************************
     ! Interpolation to correct pressure and correct temperatures
     !***********************************************************

     ! kappas

     PorT = temp(i_str)-opa_TP_grid(PT_ind_Ts_Ps,1)

     slopes = (buffer_Tl-buffer_Ts)/(opa_TP_grid(PT_ind_Tl_Ps,1)-opa_TP_grid(PT_ind_Ts_Ps,1))

     if (temp(i_str) >= temp_max) then
        opa_struc_kappas(:,:,:,i_str) = buffer_Tl
     else if (temp(i_str) <= temp_min) then
        opa_struc_kappas(:,:,:,i_str) = buffer_Ts
     else
        opa_struc_kappas(:,:,:,i_str) = &
             buffer_Ts +  slopes*PorT
     end if

  end do

end subroutine interpol_opa_ck

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to get the abundance weightes opas for ck, and for adding the continuum opas. 

subroutine mix_opas_ck(abundances,opa_struc_kappas,continuum_opa, &
     N_species,freq_len,struc_len,g_len,opa_struc_kappas_out)

  implicit none
  ! I/O
  integer, intent(in) :: N_species,freq_len,struc_len,g_len
  double precision, intent(in) :: abundances(struc_len,N_species), &
       continuum_opa(freq_len,struc_len)
  double precision, intent(in) :: opa_struc_kappas(g_len,freq_len,N_species,struc_len)
  double precision, intent(out) :: opa_struc_kappas_out(g_len,freq_len,N_species,struc_len)
  ! internal
  integer :: i_spec, i_struc, i_freq

  do i_struc = 1, struc_len
     do i_spec = 1, N_species
        opa_struc_kappas_out(:,:,i_spec,i_struc) = abundances(i_struc,i_spec)* &
             opa_struc_kappas(:,:,i_spec,i_struc)
     end do
  end do

  do i_struc = 1, struc_len
     do i_freq = 1, freq_len
        opa_struc_kappas_out(:,i_freq,1,i_struc) = &
             opa_struc_kappas_out(:,i_freq,1,i_struc) + &
             continuum_opa(i_freq,i_struc)
     end do
  end do

end subroutine mix_opas_ck

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to interpolate the CIA opacities

subroutine CIA_interpol(freq,temp,CIA_cpair_lambda,CIA_cpair_temp,CIA_cpair_alpha_grid, &
     press,mmw,mfrac,mu_part,CIA_cpair_intp_out,struc_len,freq_len)

  use constants_block
  implicit none
  ! I/O
  INTEGER, intent(in)                          :: struc_len, freq_len
  DOUBLE PRECISION, intent(in)                 :: CIA_cpair_alpha_grid(5000,30)
  DOUBLE PRECISION, intent(in)                 :: CIA_cpair_lambda(5000)
  DOUBLE PRECISION, intent(in)                 :: CIA_cpair_temp(30)
  DOUBLE PRECISION, intent(in)                 :: press(struc_len), mmw(struc_len), &
       mfrac(struc_len), mu_part
  DOUBLE PRECISION, intent(in)                 :: freq(freq_len), temp(struc_len)
  DOUBLE PRECISION, intent(out)                :: CIA_cpair_intp_out(freq_len,struc_len)
  ! Internal
  DOUBLE PRECISION                             :: buff1,buff2
  DOUBLE PRECISION                             :: lambda(freq_len), factor(struc_len)
  INTEGER                                      :: temp_ind_lower, lamb_ind_lower, i_lamb, i_struc
  

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  lambda = c_l / freq
  CIA_cpair_intp_out = 0d0

  factor = (mfrac/mu_part)**2d0*mmw/amu/(L0**2d0)*press/kB/temp

  do i_struc = 1, struc_len

     do i_lamb = 1, freq_len

!!$        if (mod(i_lamb-1,10) .NE. 0) then
!!$           CIA_cpair_intp_out(i_lamb,i_struc) = CIA_cpair_intp_out(i_lamb-1,i_struc)
!!$           cycle
!!$        end if
        
        if (lambda(i_lamb) < CIA_cpair_lambda(1)) then
           CIA_cpair_intp_out(i_lamb,i_struc) = 0d0
        else if (lambda(i_lamb) > CIA_cpair_lambda(5000)) then
           CIA_cpair_intp_out(i_lamb,i_struc) = 0d0
        else
           if (lambda(i_lamb) > CIA_cpair_lambda(4999)) then
              lamb_ind_lower = 4999
           else if (lambda(i_lamb) < CIA_cpair_lambda(2)) then
              lamb_ind_lower = 1
           else
              lamb_ind_lower = IDINT((log10(lambda(i_lamb)) - &
                   log10(0.607d0*1d-4))/(log10(250d0*1d-4) - log10(0.607d0*1d-4)) * &
                   (5d3-1d0))+1
           end if


           if (temp(i_struc) < CIA_cpair_temp(2)) then
              temp_ind_lower = 1
           else if (temp(i_struc) > CIA_cpair_temp(29)) then
              temp_ind_lower = 29
           else
              temp_ind_lower = IDINT((temp(i_struc)-100d0)/(3000d0-100d0)*(30d0-1d0))+1
           end if

           buff1 = CIA_cpair_alpha_grid(lamb_ind_lower,temp_ind_lower) + &
                (CIA_cpair_alpha_grid(lamb_ind_lower,temp_ind_lower+1) - &
                CIA_cpair_alpha_grid(lamb_ind_lower,temp_ind_lower)) / &
                (CIA_cpair_temp(temp_ind_lower+1)-CIA_cpair_temp(temp_ind_lower)) * &
                (temp(i_struc) - CIA_cpair_temp(temp_ind_lower))

           buff2 = CIA_cpair_alpha_grid(lamb_ind_lower+1,temp_ind_lower) + &
                (CIA_cpair_alpha_grid(lamb_ind_lower+1,temp_ind_lower+1) - &
                CIA_cpair_alpha_grid(lamb_ind_lower+1,temp_ind_lower)) / &
                (CIA_cpair_temp(temp_ind_lower+1)-CIA_cpair_temp(temp_ind_lower)) * &
                (temp(i_struc) - CIA_cpair_temp(temp_ind_lower))

           CIA_cpair_intp_out(i_lamb,i_struc) = buff1 + &
                (buff2-buff1)/(CIA_cpair_lambda(lamb_ind_lower+1) - &
                CIA_cpair_lambda(lamb_ind_lower)) * &
                (lambda(i_lamb)-CIA_cpair_lambda(lamb_ind_lower))

        end if

        if (CIA_cpair_intp_out(i_lamb,i_struc) < 0d0) then
           CIA_cpair_intp_out(i_lamb,i_struc) = 0d0
        end if

     end do

     CIA_cpair_intp_out(:,i_struc) = CIA_cpair_intp_out(:,i_struc)*factor(i_struc)

  end do

  !--------------
  !-- CHANGE to rather giving
  !-- the opacity at the largest / smallest
  !-- temperature grid point if temperature
  !-- is smaller or larger than the min / max
  !-- grid temperature!
  !--------------

end subroutine CIA_interpol

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to read the CIA opacities

subroutine CIA_read(cpair,opacity_path_str,CIA_cpair_lambda, &
     CIA_cpair_temp,CIA_cpair_alpha_grid)

  implicit none
  ! I/O
  CHARACTER*20, intent(in)              :: cpair
  DOUBLE PRECISION, intent(out)         :: CIA_cpair_alpha_grid(5000,30)
  DOUBLE PRECISION, intent(out)         :: CIA_cpair_lambda(5000)
  DOUBLE PRECISION, intent(out)         :: CIA_cpair_temp(30)
  CHARACTER*150, intent(in)             :: opacity_path_str
  ! internal
  INTEGER                               :: i,j

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  open(unit=10,file=trim(adjustl(opacity_path_str)) &
       //'/opacities/continuum/CIA/'//trim(adjustl(cpair))//'/temps.dat')
  do i = 1,30
     read(10,*) CIA_cpair_temp(i)
  end do
  close(10)

  open(unit=11,file=trim(adjustl(opacity_path_str)) &
       //'/opacities/continuum/CIA/'//trim(adjustl(cpair)) &
       //'/CIA_'//trim(adjustl(cpair))//'_final.dat')
  read(11,*)
  read(11,*)
  do i = 1, 5000
     read(11,'(G22.12)',advance='no') CIA_cpair_lambda(i)
     do j = 1, 29
        read(11,'(G22.12)',advance='no') CIA_cpair_alpha_grid(i,j)
     end do
     read(11,'(G22.12)') CIA_cpair_alpha_grid(i,30)
  end do
  close(11)

end subroutine CIA_read

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to get the length of the opacity arrays in the high-res case

subroutine get_arr_len_array_bords(wlen_min_read, wlen_max_read, &
     file_path, arr_len, arr_min, arr_max)

  implicit none

  ! I/O
  double precision, intent(in) :: wlen_min_read, wlen_max_read
  character*400, intent(in)    :: file_path
  integer, intent(out)         :: arr_len, arr_min, arr_max
  ! Internal
  double precision :: curr_wlen, last_wlen
  integer          :: curr_int

  ! open wavelength file
  open(file=trim(adjustl(file_path)), unit=10, form = 'unformatted', &
       ACCESS='stream')

  ! to contain the current wavelength index
  curr_int = 1

  ! to contain the the minimum and the maximum wavelength index
  ! to be used for reading in the opacities and wavelengths later.
  arr_min = -1
  arr_max = -1

  ! to contain the wavelength of the previous line reading
  last_wlen = 0d0

  do while (1>0)

     read(10,end=123) curr_wlen

     if ((curr_int .EQ. 1) .AND. (curr_wlen > wlen_min_read)) then
        write(*,*) 'ERROR! Desired minimum wavelength is too small!'
        STOP
     end if

     ! look for minimum index, bracketing the desired range
     if (arr_min .EQ. -1) then
        if ((curr_wlen > wlen_min_read) .AND. &
             (last_wlen < wlen_min_read)) then
           arr_min = curr_int - 1
        end if
     end if

     ! look for maximum index, bracketing the desired range
     if (arr_min .NE. -1) then
        if ((curr_wlen > wlen_max_read) .AND. &
             (last_wlen < wlen_max_read)) then
           arr_max = curr_int
           EXIT
        end if
     end if

     last_wlen = curr_wlen

     curr_int = curr_int + 1
  end do

123 close(10)

  if ((arr_min .EQ. -1) .OR. (arr_max .EQ. -1)) then

     write(*,*) 'ERROR! Desired wavelength range is too large,'
     write(*,*) 'or not contained within the tabulated opacity' &
          // ' wavelength range.'
     STOP

  end if

  arr_len = arr_max - arr_min + 1

end subroutine get_arr_len_array_bords

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to read the wavelength array in the high-res case

subroutine read_wlen(arr_min, arr_max, arr_len, file_path, wlen)

  implicit none
  ! I/O
  integer, intent(in)           :: arr_min, arr_max
  integer, intent(in)           :: arr_len
  character*400, intent(in)     :: file_path
  double precision, intent(out) :: wlen(arr_len)

  integer          :: i_lamb

  open(unit=49, file=trim(adjustl(file_path)), &
       form = 'unformatted', ACCESS='stream')

  read(49, pos = (arr_min-1)*8+1) wlen(1)
  do i_lamb = 2, arr_len
     read(49) wlen(i_lamb)
  end do

  close(49)

end subroutine read_wlen

!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################
!!$ #########################################################################

!!$ Subroutine to read the kappa array in the high-res case

subroutine read_kappa(arr_min, arr_max, arr_len, file_path, kappa)

  implicit none
  ! I/O
  integer, intent(in)           :: arr_min, arr_max
  integer, intent(in)           :: arr_len
  character*400, intent(in)     :: file_path
  double precision, intent(out) :: kappa(arr_len)

  integer          :: i_lamb

  open(unit=49, file=trim(adjustl(file_path)), &
       form = 'unformatted', ACCESS='stream')

  read(49, pos = (arr_min-1)*8+1) kappa(1)
  do i_lamb = 2, arr_len
     read(49) kappa(i_lamb)
  end do

  close(49)

end subroutine read_kappa
