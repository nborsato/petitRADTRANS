Adding opacities
================

petitRADTRANS has an extensive database of line opacities. However, it is very
likely that we are missing the one atom / molecule that you want.
Here we give an example on how to calculate opacities and reformatting
them for use in petitRADTRANS. In this way you will be able to add
additional species, that we do not have added yet. If you are
particularly nice, you can share these opacities with us, we will then
make them available to the other petitRADTRANS users via this website, while properly
attributing your contribution.

From line lists to opacities (using ExoCross)
_____________________________________________

First, a line list needs to be converted into actual opacities.
In this example we will show you how to do this using ExoCross, the
open-source opacity calculator of the `Exomol`_ database.
ExoCross can be downloaded `here <https://github.com/Trovemaster/exocross>`_, is described in
`Yurchenko et al. (2018)`_ and documented `here
<https://exocross.readthedocs.io>`_. 

.. _Exomol: http://www.exomol.com
.. _Yurchenko et al. (2018): https://arxiv.org/abs/1801.09803

First, download the ExoCross source, go into the folder containing the
source and the makefile called "makefile". Adapt that to your liking.
For example, if you have the gfortran compiler, but not ifort, make
sure that the flag using ifort is commented out, and that it uses
gfortran. The relevant lines in "makefile" should look like this:

.. code-block:: bash

    #FOR  = ifort
    #FFLAGS =  -O3 -qopenmp -traceback  -ip                                                                                        
    FOR = gfortran
    FFLAGS = -O2 -fopenmp -std=f2008

Then, build ExoCross by typing ``make`` in the terminal. Sometimes the compiler will
complain that lines within the ExoCross source are too long. Just open
the source and introduce a line break there manually, like this:

.. code-block:: fortran

    ! This is an example for a line that is too long
    DOUBLE PRECISION :: very_long_variable_name_number_one, very_long_variable_name_number_two, very_long_variable_name_number_three

    ! This is how you introduce line breaks
    DOUBLE PRECISION :: very_long_variable_name_number_one, &
       very_long_variable_name_number_two, &
       very_long_variable_name_number_three

So the ``&`` is th line break operator. After fixing this, recompile
using ``make``.

In this example we will calculate the opacities of the NaH molecule.
All necessary files for calculating opacities can be found on the Exomol
website, just `click here`_.

.. _click here: http://www.exomol.com/data/molecules/NaH/23Na-1H/Rivlin/

The following files need to be downloaded:

- 23Na-1H__Rivlin.states.bz2
- 23Na-1H__Rivlin.trans.bz2
- 23Na-1H__Rivlin.pf

Please unzip the .bz2 files before use.

Next, make an input file for carrying out the calculations, in this
example we call it NaH_input.inp. This is what it looks like:

.. code-block:: bash

    absorption
    voigt
    verbose 3
    offset 60.
    mass 24
    temperature 1000.000000
    pressure 0.00001
    range 39. 91000.
    R 1000000
    pffile 23Na-1H__Rivlin.pf
    output NaH_1000K_1em5bar.out
    states 23Na-1H__Rivlin.states
    transitions
      "23Na-1H__Rivlin.trans"
    end
    species
      0 gamma 0.06 n 0.5 t0 296 ratio 1.
    end

This calculates the opacity of NaH with the following settings

- ``offset `` results in a line cutoff of 60 :math:`{\rm
  cm}^{-1}`. While being an important effect that also speeds up
  calculations, the choice of a cutoff is often arbitrary because the
  physics behind it remain difficult to model, see, for example the
  discussion in `Grimm & Heng
  (2015)`_. Here we use the equivalent width of the line decrease
  function given by `Hartmann et al. (2002)`_, for :math:`\rm CH_4`
  broadened by :math:`\rm H_2`.
- NaH has a mass of 24 (in amu)
- The opacity is calculated at a temperature of 1000 K
- The opacity is calculated at a pressure of :math:`10^{-5}` bar
- The opacity is calculated in the range from 39 to 91000 :math:`{\rm
  cm}^{-1}`. This corresponds to a wavelength range from 0.1099 to
  256.4103 micron, therefore bracketing the full petitRADTRANS
  wavelength range (0.11 to 250 micron at low resolution). This large
  a range is needed. Therefore, do not change this. Note that the opacities in
  the high-resolution mode of petitRADTRANS ultimately only go from
  0.3 to 28 microns.
- The resolution of the calculations carried out here is for a
  wavelength spacing of :math:`\lambda/\Delta\lambda=10^6`.
- The ``pfile`` line gives the relative path to the partition function
  file, that you have already downloaded from Exomol.
- The ``states`` line gives the relative path to the states
  file, that you have already downloaded from Exomol.
- The lines below ``transitions`` line give the relative paths to the transition
  files, that you have already downloaded from Exomol. For NaH this is
  only one file. For molecules with a lot more lines this can be
  multiple files.
- The lines below ``species`` define the pressure broadening to be
  used. This pressure boradening (width of the Lorentz profile) is of
  the form :math:`\gamma \cdot (T_{0}/T)^n ({\rm ratio}\cdot
  P/{\rm 1 \ bar})`, in units of :math:`\rm cm^{-1}`.  The choice here is a compromise between the
  various values reported for the broadening by :math:`\rm H_2/He` of
  various absorbers, e.g. in `Amundsen et al. (2014)`_, `Gharib-Nezhad &
  Line (2018)`_. Also see the text around Equation 12 in `Sharp &
  Burrows (2007)`_ for more information. Sometimes more detailed
  broadening information is available on Exomol, `see here`_.
  
.. _Hartmann et al. (2002): http://adsabs.harvard.edu/abs/2002JQSRT..72..117H
.. _Grimm & Heng (2015): https://arxiv.org/abs/1503.03806
.. _Amundsen et al. (2014): https://arxiv.org/abs/1402.0814
.. _Gharib-Nezhad & Line (2018): https://arxiv.org/abs/1809.02548v2
.. _Sharp & Burrows (2007): https://arxiv.org/abs/astro-ph/0607211
.. _see here: http://www.exomol.com/data/data-types/broadening_coefficients/

If more detailed broadening information is avaiable (not for NaH) you can replace
the lines below ``species`` with something like

.. code-block:: bash
		
    species
      0 gamma 0.06 n 0.5 t0 296 file path_toH2_broadening_information_file model J ratio 0.860000
      1 gamma 0.06 n 0.5 t0 296 file path_toHe_broadening_information_file model J ratio 0.140000
    end

The above setting is for a primordial composition atmosphere, where
:math:`\rm H_2` and He roughly make up 86 % and 14 % of the
atmosphere, respectively (i.e. these are volume mixing ratios, not
mass fractions). The :math:`\gamma` and :math:`n` values given before
the path to the boradening files are what is used for rotational
quantum numbers (:math:`J`) not covered by the broadening files.

Finally, the opacities are calculated by running ExoCross from the
terminal command line via

.. code-block:: bash

     ./xcross.exe < NaH_input.inp > test_run.out

The resulting wavelength-dependent opacity will be in the "NaH_1000K_1em5bar.out.xsec" file, in our
example here.
In the end quite a few opacity points need to be calculated for
petitRADTRANS (at 130 or 200 different pressure-temperature
conbinations, see below). This is doable on a local machine for smaller
linelists such as NaH, but may require the use of a cluster for much
larger linelists. There also exsists the so-called superline
treatment `(see Yurchenko et al. 2018)`_
, where multiple lines are combined into one, this can speed
up calculations a lot, but is not recommended if you want to calculate
high-resolution spectra with petitRADTRANS (because line positions
will shift if multiple lines are combined into one on a fixed
wavelength grid).

.. _(see Yurchenko et al. 2018): https://arxiv.org/abs/1801.09803

	
Preparing ExoCross opacities for petitRADTRANS
______________________________________________


For creating opacities for use in petitRADTRANS, calculate the
molecular opacities from Exomol with ExoCross using the settings
outlined above. Change parameters where applicable (temperature,
pressure, molecule mass, broadening information...).

The opacities need to be calculated at the 130 pressure temperature points
of petitRADTRANS which you can find in the file
`PTgrid.dat <https://gitlab.com/mauricemolli/petitRADTRANS/blob/b4e305de65f298c5c0b09568756aa005477489b2/docs/content/files/PTgrid.dat>`_. Temeratures go from 80 up to 3000 K,
in a log-uniform way. If you want to be ready for the future, please calculate opacities
using `PTgrid_new.dat <https://gitlab.com/mauricemolli/petitRADTRANS/blob/b4e305de65f298c5c0b09568756aa005477489b2/docs/content/files/PTgrid_new.dat>`_, where we have added a
few more points at high temperatures (increasing the temperature resolution there) and extend
the temperature range to 4000 K (note that currently petitRADTRANS sets
:math:`\kappa(T>3000 K)` to :math:`\kappa(T=3000 K)` for the opacity
:math:`\kappa`, if tempertatures get too high). The new grid has 200 points in total. The ability of
petitRADTRANS to use the high-temperature grid (`PTgrid_new.dat <https://gitlab.com/mauricemolli/petitRADTRANS/blob/b4e305de65f298c5c0b09568756aa005477489b2/docs/content/files/PTgrid_new.dat>`_) will
be added ASAP. Shoot us a `email`_ to find out when / pressuring us to do
this even more quickly.

.. important::
   
    For later use in petitRADTRANS it is important that the opacity
    files have the correct names. Every molecule needs to be assigned
    a random two digit integer (e.g. "06"). It does not matter if two
    different molecules have the same number. Then, the opacity file
    at (for example) T = 200 K, P = 1 bar must be called
    "sigma_06_200.K_1.000000bar.dat". The exact names of the files
    used in petitRADTRANS can be found `in this file here`_. petitRADTRANS will be
    looking for these files, and throw an error message and crash if
    they are not named properly.

.. _email: molliere@mpia.de
.. _in this file here: https://gitlab.com/mauricemolli/petitRADTRANS/blob/b4e305de65f298c5c0b09568756aa005477489b2/docs/content/files/PTnames.dat

Now, let's turn towards preparing the ExoCross results for
petitRADTRANS. We will assume that you have calculated the opacites at
all 130 (or 200) pressure-temperature points. The high-resolution
wavelength setup between ExoCross and our
classical petitCODE/petitRADTRANS opacity calculator is slightly
different. ExoCross' wavelength spacing varies a bit around the
user-defined resolution, whereas our routines preparing the opacity
files for petitRADTRANS assume that the wavelength spacing is exactly
:math:`\lambda/\Delta\lambda=10^6`, from 0.11 to 250 microns.
Hence we will first have to rebin the ExoCross results to the
petitCODE/petitRADTRANS grid. To this end, please download the
petitRADTRANS high resolution grid (`wlen_petitRADTRANS.dat`_).

.. _`wlen_petitRADTRANS.dat`: https://www.dropbox.com/s/2lyo8ot3nq4rx43/wlen_petitRADTRANS.dat?dl=0

Next, rebin all ExoCross opacity files to that wavelength file, like
shown below, using Python, here for simplicity we use the NaH opacity file
calculated above.

.. code-block:: bash

    import numpy as np
    from scipy.interpolate import interp1d
    
    # Read the opacity file from ExoCross
    dat = np.genfromtxt('NaH_1000K_1em5bar.out.xsec')
    wavelength = 1./dat[:,0]
    sigma = dat[:,1]

    # Invert them to go from a accending wavenumber ordering
    # to an accending wavelength ordering.
    wavelength = wavelength[::-1]
    sigma = sigma[::-1]

    # Read the fiducial petitRADTRANS wavelength grid
    wavelength_petit = np.genfromtxt('wlen_petitRADTRANS.dat')

    # Interpolate the ExoCross calculation to that grid
    sig_interp = interp1d(wavelength, sigma)
    sig_interpolated_petit = sig_interp(wavelength_petit)

    # Save rebinned calculation
    np.savetxt('NaH_1000K_1em5bar_petit_grid.dat', \
       np.column_stack((wavelength_petit, \
                                    sig_interpolated_petit)))

Now we can create the correlated-k tables (or just "k-tables") and high-resolution opacity files from
these formatted files. Please `email`_ us to get the relevant Fortran
source to do this, we will send you four files called

- calc_k_g_r1000_ptrad.f90: this converts the opacity data to
  petitRADTRANS k-tables (these are the opacities for the
  low-resolution mode of petitRADTRANS, at :math:`\lambda/\Delta\lambda=1000`.
- retrieval_NP_16_ggrid.dat: this is the 16-point Gaussian quadrature
  grid that petitRADTRANS uses as the g-coordinate for the k-tables.
- make_short.f90: this cuts the opacities to the right 0.3 to 28
  micron range for the high-resolution calculations
  :math:`\lambda/\Delta\lambda=10^6`.
- short_stream_lambs_mass.dat: input file for make_short.f90.

.. _email: molliere@mpia.de

You do not need to understand anything about k-tables to do this step
here, we just wanted to explain what the routines are for.

To start, put the names of all opacity files you want to convert into a file called
"sigma_list.ls". Do not include the paths to these files, just the
file names. Hence will have to run the Fortran conversion routines in the
folder where the opacity files are. In our simple example (just one
NaH file at 1000 K and :math:`10^{-5}` bar, its content looks like this:

.. code-block:: bash
		
    NaH_1000K_1em5bar_petit_grid.dat

Let's start with the k-table calculation, for the low-resolution
opacity mode of petitRADTRANS. Open calc_k_g_r1000_ptrad.f90 and
modify it to have the correct mass for the molecular species that you
are interested in (NaH has 24 amu, so just put 24, like below):

.. code-block:: fortran

    ! (c) Paul Molliere 2014
    
     program calc_k_g
    
      implicit none

      !-----------------------------------------------------------
      !            |||               |||                |||      !
      !           \|||/             \|||/              \|||/     !
      !             v                 v                  v       !    
      !----------------------------------------------------------!
      !----------------------------------------------------------!
      ! DO NOT FORGET TO CHANGE THE MASS OF THE MOLECULE
      ! EVERY TIME!!!
      DOUBLE PRECISION, parameter   :: mol_mass_amu = 24d0  !<---!
      !----------------------------------------------------------!
      !----------------------------------------------------------!
      !             ^                 ^                  ^       !
      !           /|||\             /|||\              /|||\     !
      !            |||               |||                |||      !
      !----------------------------------------------------------!
 
Next, compile the Fortran source:

.. code-block:: bash
		
    gfortran -o calc_k_g_r1000_ptrad calc_k_g_r1000_ptrad.f90

Lastly, create a folder called kappa_gs_r1000. Now, take care that the opacity files, the compiled Fortran routine,
sigma_list.ls, retrieval_NP_16_ggrid.dat and the kappa_gs_r1000 folder
are all in the same folder. And that you are in this folder. Type

.. code-block:: bash
		
    ./calc_k_g_r1000_ptrad

and all k-tables will be generated and placed into the kappa_gs_r1000
folder.

For the high resolution mode, generate a folder called "short_stream".
Next open the short_stream_lambs_mass.dat file and adapt its content
to have the correct molecule mass. Do not change the wavelength
boundary values in this file. For NaH, with mass 24, it should look
like this:

.. code-block:: bash
		
    # Minimum wavelength in cm
    0.3d-4
    # Maximum wavelength in cm
    28d-4
    # Molecular mass in amu
    24d0

Next, compile the high-resolution opacity conversion routine:

.. code-block:: bash
		
    gfortran -o make_short make_short.f90

Now, again take care that the opacity files, the compiled Fortran routine,
sigma_list.ls, short_stream_lambs_mass.dat and the short_stream folder
are all in the same folder. And that you are in this folder. Type

.. code-block:: bash
		
    ./make_short

and all high resolution opacity tables will be generated and placed into the short_stream
folder.

Installing the new opacity files in petitRADTRANS
_________________________________________________

The new opacity files are now ready to be installed. Before that
create a file called "molparam_id.txt" with the following content

.. code-block:: bash
		
    #### Species ID (A2) format
    06
    #### molparam value
    1.0

Simply exchange the "06" two-digit integer with the one that you have
chosen for your molecule (or leave it at 06 if you chose to keep it).
Copy the "molparam_id.txt" file to the short_stream and kappa_gs_r1000
folders. Now we are ready for installation. In the folder where
petitRADTRANS is installed, there also is a input_data folder. To
install a new species (e.g. NaH), create a folder called NaH in the
input_data/opacities/lines/corr_k/ and
input_data/opacities/lines/line_by_line folders. Copy the contents of
the kappa_gs_r1000 and short_stream folders to the NaH folders in the
corr_k and line_by_line folders, respectively. The opacities are now
installed and ready for use!
