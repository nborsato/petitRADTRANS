.. petitRADTRANS documentation master file, created by
   sphinx-quickstart on Tue Jan 15 15:07:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

petitRADTRANS documentation
=========================================

Welcome to the **petitRADTRANS** (pRT) documentation. pRT is a
Python package for calculating transmission and emission spectra
of exoplanets, at low (:math:`\lambda/\Delta\lambda=1000`) and high
(:math:`\lambda/\Delta\lambda=10^6`) resolution, for clear and cloudy
atmospheres.

.. important::

   In addition to transmission spectra, **pRT now includes scattering also for emission spectra**, if specifically turned on (note that scattering increases the runtime), see `Scattering for Emission Spectra <content/notebooks/emis_scat.html>`_.

petitRADTRANS is available under the MIT License, and documented in
`Mollière et al. (2019) <https://arxiv.org/abs/1904.11504>`_, for the general code, and `Mollière et al. (2020) <https://arxiv.org/abs/2006.09394>`_, Alei et al. (in prep.), for the scattering implementation. Please cite these papers if you make use of petitRADTRANS in your work.

.. _contact: molliere@mpia.de

This documentation webpage currently contains an installation guide, a
tutorial, a first code documentation, and an implemented retrieval
example for mock JWST emission and transmission spectra.
Also, we give a tutorial on how to include opacities that may be
missing from our database.

News
____

**December 2020: stellar and planetary surface scattering added**
    pRT now includes the scattering of the incoming stellar flux for irradiated planets. Also a scattering surface for terrestrial planets has been added, see `Scattering for Emission Spectra <content/notebooks/emis_scat.html>`_. The surface albedo and emissivity are both freely tunable parameters (as a function of wavelength). Surface scattering is treated to be isotropic (that is, assuming that the surface is Lambertian).

**September 2020: self-scattering for emission spectra and chemical equilibrium interpolation now available**
    pRT now includes scattering also for emission spectra, if specifically turned on (note that scattering increases the runtime), see `Scattering for Emission Spectra <content/notebooks/emis_scat.html>`_. Currently the self-scattering by the planetary atmosphere is included, which is appropriate for, for example, brown dwarf and directly imaged atmospheres. In addition, you can now download our chemical equilibrium interpolation package, which is documented in `Interpolating chemical equilibrium abundances <content/notebooks/poor_man.html>`_
    
**September 2020: petitRADTRANS opacities available on the Exomol website**
    Opacity tables created specifically in the petitRADTRANS format
    are now available on the `Exomol website <http://www.exomol.com/data/data-types/opacity/>`_, also see `Chubb et al. (2020) <https://arxiv.org/abs/2009.00687>`_ for the accompanying paper.
    The opacities can be installed in petitRADTRANS in an easy plug-and-play fashion.
    Please see Section `Adding opacities <content/opa_add.html>`_ for more information.

**September 2020: More high-temperature atom and ion opacities available**
    We have added more atom and ion opacities, bringing the total list to
    Al, AlII, AlIII, AlIV, AlV, AlVI, B, BII, BIII, Be, BeII, C, CII, CIII, CIV,
    Ca, CaII, Cr, Fe, FeII, K, KII, KIII, KIV, KV, KVI, Li, Mg, MgII, MgIII, MgIV,
    MgV, MgVI, N, NII, NIII, NIV, NV, Na, NaII, NaIII, NaIV, NaV, NaVI, Si, SiII, Ti,
    TiII, V, VII, Y

**May 2019: high-temperature atom and ion opacities now available**
    We have added the opacities of Fe, Fe+, Mg, Mg+, Li, Ca, Ca+,
    Si, Si+, O, O+, Al, Al+, Ti, Ti+, V and V+, up to temperatures of
    4000 K. As usual, if the atmospheric temperatures increase
    above 4000 K, petitRADTRANS will use the absorbers respective
    opacities at 4000 K. **Please make sure to install the latest
    petitRADTRANS version to make use of the high-temperature
    points of the new opacity tables!**

Developers
___________

- Paul Mollière
- Eleonora Alei
- Evert Nasedkin

Contributors
________________
- Karan Molaverdikhani
- Mantas Zilinskas



.. toctree::
   :maxdepth: 2
   :caption: Guide:

   content/installation
   content/tutorial
   content/available_opacities
   content/retrieval_examples
   content/opa_add

.. toctree::
   :maxdepth: 2
   :caption: Code documentation:

   content/code
