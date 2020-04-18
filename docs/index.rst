.. petitRADTRANS documentation master file, created by
   sphinx-quickstart on Tue Jan 15 15:07:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

petitRADTRANS documentation
=========================================

Welcome to the **petitRADTRANS** documentation. petitRADTRANS is a
Python package for calculating transmission and emission spectra
of exoplanets, at low (:math:`\lambda/\Delta\lambda=1000`) and high
(:math:`\lambda/\Delta\lambda=10^6`) resolution, for clear and cloudy
atmospheres.

.. important::
   Scattering is currently not included for
   emission spectra in petitRADTRANS (but it is for the transmission
   spectra). We intend to migrate the scattering implementation of
   `petitCODE <http://www.mpia.de/homes/molliere/#petitcode>`_
   to petitRADTRANS soon.

petitRADTRANS is available under the MIT License, and documented in
`Mollière et al. (2019) <https://arxiv.org/abs/1904.11504>`_.
Please cite this paper if you make use of petitRADTRANS in your work.

.. _contact: molliere@mpia.de

This documentation webpage currently contains an installation guide, a
tutorial, a first code documentation, and an implemented retrieval
example for mock JWST emission and transmission spectra.
Also, we give a tutorial on how to calculate opacities that may be
missing from our database.

News
____

**27 May 2019: high-temperature atom and ion opacities now available**
    We have added the opacities of Fe, Fe+, Mg, Mg+, Li, Ca, Ca+,
    Si, Si+, O, O+, Al, Al+, Ti, Ti+, V and V+, up to temperatures of
    4000 K. As usual, if the atmospheric temperatures increase
    above 4000 K, petitRADTRANS will use the absorbers respective
    opacities at 4000 K. **Please make sure to install the latest
    petitRADTRANS version to make use of the high-temperature
    points of the new opacity tables!**

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
