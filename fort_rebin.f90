!!$*****************************************************************************
!!$*****************************************************************************
!!$*****************************************************************************
!!$ fort_rebin.f90: utility function to rebin synthetic spectrum to data spectrum
!!$ for the petitRADTRANS radiative transfer package.
!!$ 
!!$ Copyright 2016-2018, Paul Molliere
!!$ Maintained by Paul Molliere, molliere@strw.leidenunivl.nl
!!$ Status: under development
!!$*****************************************************************************
!!$*****************************************************************************
!!$*****************************************************************************

subroutine rebin_spectrum(nu_synth,flux_synth,nu_obs, &
     reb_synth_flux,len_obs,len_synth)

  implicit none
  ! I/O
  integer, intent(in) :: len_obs, len_synth
  double precision, intent(in) :: nu_obs(len_obs)
  double precision, intent(in) :: nu_synth(len_synth),flux_synth(len_synth)
  double precision, intent(out) :: reb_synth_flux(len_obs)
  ! internal
  double precision :: nu_obs_bin_bords(len_obs+1)
  integer :: i_nu, i_nu_synth
  double precision :: nu1, nu2, y_ax, del_nu, ddel_nu, slope, add_fl, min_flux, max_flux

  reb_synth_flux = 0d0

  ! Check stuff...
  if (len_obs < 2) then
     STOP 'rebin.f90 ERROR: Coarse wavelength grid needs at least two wavelength points!'
  end if

  ! Get intp bin bords
  nu_obs_bin_bords(1) = nu_obs(1)-(nu_obs(2)-nu_obs(1))/2d0
  do i_nu = 2, len_obs
     nu_obs_bin_bords(i_nu) = nu_obs(i_nu)-(nu_obs(i_nu)-nu_obs(i_nu-1))/2d0
  end do
  nu_obs_bin_bords(len_obs+1) = nu_obs(len_obs)+(nu_obs(len_obs)-nu_obs(len_obs-1))/2d0

  if ((nu_synth(1) >= nu_obs_bin_bords(1)) .OR. (nu_synth(len_synth) <= nu_obs_bin_bords(len_obs+1))) then
     STOP 'rebin.f90 ERROR: Fine wavelength grid needs to extend at least half'// &
          ' a bin width further than coarse grid (on both sides)!'
  end if

  ! Start interpolation
  i_nu_synth = 1
  
  do i_nu = 1, len_obs

     do while (nu_synth(i_nu_synth) < nu_obs_bin_bords(i_nu))
        i_nu_synth = i_nu_synth+1
     end do

     del_nu = 0d0
     
     do while (nu_synth(i_nu_synth) < nu_obs_bin_bords(i_nu+1))

        nu1 = max(nu_synth(i_nu_synth-1),nu_obs_bin_bords(i_nu))
        nu2 = min(nu_synth(i_nu_synth),nu_obs_bin_bords(i_nu+1))
        ddel_nu = nu2-nu1
        slope = (flux_synth(i_nu_synth)-flux_synth(i_nu_synth-1))/ &
             (nu_synth(i_nu_synth)-nu_synth(i_nu_synth-1))
        y_ax = flux_synth(i_nu_synth-1)

        min_flux = min(flux_synth(i_nu_synth),flux_synth(i_nu_synth-1))
        max_flux = max(flux_synth(i_nu_synth),flux_synth(i_nu_synth-1))

        ! y+b*(x-a) , int from x1 .. x2
        ! y*(x2-x1)-b*a*(x2-x1)+b/2d0*(x2**2d0-x1**2d0)
        add_fl = (y_ax-slope*nu_synth(i_nu_synth-1))*(nu2-nu1) + &
             slope*(nu2-nu1)*(nu2+nu1)/2d0

        reb_synth_flux(i_nu) = reb_synth_flux(i_nu) + add_fl
        
        del_nu = del_nu + ddel_nu

        i_nu_synth = i_nu_synth + 1

     end do

     nu1 = max(nu_synth(i_nu_synth-1),nu_obs_bin_bords(i_nu))
     nu2 = min(nu_synth(i_nu_synth),nu_obs_bin_bords(i_nu+1))
     ddel_nu = nu2-nu1
     slope = (flux_synth(i_nu_synth)-flux_synth(i_nu_synth-1))/ &
          (nu_synth(i_nu_synth)-nu_synth(i_nu_synth-1))
     y_ax = flux_synth(i_nu_synth-1)

     min_flux = min(flux_synth(i_nu_synth),flux_synth(i_nu_synth-1))
     max_flux = max(flux_synth(i_nu_synth),flux_synth(i_nu_synth-1))

     ! y+b*(x-a) , int from x1 .. x2
     ! y*(x2-x1)-b*a*(x2-x1)+b/2d0*(x2**2d0-x1**2d0)
     add_fl = (y_ax-slope*nu_synth(i_nu_synth-1))*(nu2-nu1) + &
          slope*(nu2**2d0-nu1**2d0)/2d0

     reb_synth_flux(i_nu) = reb_synth_flux(i_nu) + add_fl

     del_nu = del_nu + ddel_nu

     reb_synth_flux(i_nu) = reb_synth_flux(i_nu)/del_nu

  end do
  
end subroutine rebin_spectrum
