.. _avail_opas:

Available opacity species
=========================

Line absorbers
______________

**Line absorbers, low resolution mode** (``"c-k"``, with :math:`\lambda/\Delta\lambda=1000`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. important::
   In low resolution mode (``"c-k"``), most of the molecular opacitites are calculated considering only the main isotopologue. This is different only for CO and TiO, where the contribution of all isotopologues is considered. For CO because the secondary isotopes of carbon, for example :math:`\rm ^{13}C`, are quite abundant when compared to the main isotope, that is :math:`\rm ^{12}C/^{13}C\sim 100`, and because CO has very strong and sparse lines. Not including these lines therefore has a noticeable effect already at low resolution. For TiO all isotopologues are included because the relative ratios between the Ti isotopes are quite large. Apart from these two species, the main isotopologue treatment compared very well to codes including all isotopologues, at this low resolution, see `Baudino et al. (2017) <http://adsabs.harvard.edu/abs/2017ApJ...850..150B>`_.

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Species name
     - Required in abuncance dictionary
     - Description
   * - C2H2
     - C2H2
     - Main isotopologue, HITRAN
   * - CH4
     - CH4
     - Main isotopologue, Exomol
   * - CO
     - CO
     - Main isotopologue, HITEMP/Kurucz
   * - CO2
     - CO2
     - Main isotopologue, HITEMP
   * - CO_all_iso
     - CO_all_iso
     - All isotopologues, HITEMP/Kurucz
   * - H2
     - H2
     - Main isotopologue, HITRAN
   * - H2O
     - H2O
     - Main isotopologue, HITEMP
   * - H2S
     - H2S
     - Main isotopologue, HITRAN
   * - HCN
     - HCN
     - Main isotopologue, Exomol
   * - HDO
     - HDO
     - Main isotopologue, HITRAN
   * - K
     - K
     - Main isotopologue, VALD, Allard wings
   * - K_lor_cut
     - K_lor_cut
     - Main isotopologue, VALD, Lorentzian wings
   * - NH3
     - NH3
     - Main isotopologue, Exomol
   * - NH3_HITRAN
     - NH3_HITRAN
     - Main isotopologue, HITRAN
   * - Na
     - Na
     - Main isotopologue, VALD, Allard wings
   * - Na_lor_cut
     - Na_lor_cut
     - Main isotopologue, VALD, Lorentzian wings
   * - O3
     - O3
     - Main isotopologue, HITRAN
   * - OH
     - OH
     - Main isotopologue, HITEMP
   * - PH3
     - PH3
     - Main isotopologue, Exomol
   * - PH3_HITRAN
     - PH3_HITRAN
     - Main isotopologue, HITRAN
   * - SiO_main_iso
     - SiO_main_iso
     - Main isotopologue, Exomol
   * - TiO
     - TiO
     - All isotopologues, B. Plez
   * - VO
     - VO
     - Main isotopologue, B. Plez

**Line absorbers, high resolution mode** (``"lbl"``, with :math:`\lambda/\Delta\lambda=10^6`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 10 80
   :header-rows: 1

   * - Species name
     - Required in abundance dictionary
     - Description
   * - C2H2_main_iso
     - C2H2_main_iso
     - Main isotopologue, HITRAN
   * - CH4_212
     - CH4_212
     - :math:`\rm CH_3D`, HITRAN
   * - CH4_main_iso
     - CH4_main_iso
     - Main isotopologue, Exomol
   * - CO2_main_iso
     - CO2_main_iso
     - Main isotopologue, HITEMP
   * - CO_27
     - CO_27
     - :math:`\rm ^{12}C^{17}O`, HITRAN
   * - CO_28
     - CO_28
     - :math:`\rm ^{12}C^{18}O`, HITRAN
   * - CO_36
     - CO_36
     - :math:`\rm ^{13}C^{16}O`, HITRAN
   * - CO_37
     - CO_37
     - :math:`\rm ^{13}C^{17}O`, HITRAN
   * - CO_38
     - CO_38
     - :math:`\rm ^{13}C^{18}O`, HITRAN
   * - CO_all_iso
     - CO_all_iso
     - All isotopologues
   * - CO_main_iso
     - CO_main_iso
     - Main isotopologue, HITEMP
   * - H2O_162
     - H2O_162
     - :math:`\rm HDO`, HITRAN
   * - H2O_171
     - H2O_171
     - :math:`\rm H_2 \ ^{17}O`, HITRAN
   * - H2O_172
     - H2O_172
     - :math:`\rm HD^{17}O`, HITRAN
   * - H2O_181
     - H2O_181
     - :math:`\rm H_2 \ ^{18}O`, HITRAN
   * - H2O_182
     - H2O_182
     - :math:`\rm HD^{18}O`, HITRAN
   * - H2O_main_iso
     - H2O_main_iso
     - Main isotopologue, HITEMP
   * - H2S_main_iso
     - H2S_main_iso
     - Main isotopologue, HITRAN
   * - H2_12
     - H2_12
     - :math:`\rm HD`, HITRAN
   * - H2_main_iso
     - H2_main_iso
     - Main isotopologue, HITRAN
   * - HCN_main_iso
     - HCN_main_iso
     - Main isotopologue, Exomol
   * - K
     - K
     - Main isotopologue, VALD, Allard wings
   * - NH3_main_iso
     - NH3_main_iso
     - Main isotopologue, Exomol
   * - Na
     - Na
     - Main isotopologue, VALD, Allard wings
   * - O3_main_iso
     - O3_main_iso
     - Main isotopologue, HITRAN
   * - PH3_main_iso
     - PH3_main_iso
     - Main isotopologue, Exomol
   * - SiO_main_iso
     - SiO_main_iso
     - Main isotopologue, Exomol
   * - TiO_all_iso
     - TiO_all_iso
     - All isotopologues, B. Plez
   * - VO
     - VO
     - Main isotopologue, B. Plez

Cloud opacities
_______________

.. list-table::
   :widths: 10 10 80
   :header-rows: 1
		 
   * - Species name
     - Required in abundance dictionary
     - Description
   * - Al2O3(c)_cm
     - Al2O3(c)
     - Crystalline, Mie scattering (spherical)
   * - Al2O3(c)_cd
     - Al2O3(c)
     - Crystalline, DHS (irregular shape)
   * - Fe(c)_am
     - Fe(c)
     - Amorphous, Mie scattering (spherical)
   * - Fe(c)_ad
     - Fe(c)
     - Amorphous, DHS (irregular shape)
   * - Fe(c)_cm
     - Fe(c)
     - Crystalline, Mie scattering (spherical)
   * - Fe(c)_cd
     - Fe(c)
     - Crystalline, DHS (irregular shape)
   * - H2O(c)_cm
     - H2O(c)
     - Crystalline, Mie scattering (spherical)
   * - H2O(c)_cd
     - H2O(c)
     - Crystalline, DHS (irregular shape)
   * - KCL(c)_cm
     - KCL(c)
     - Crystalline, Mie scattering (spherical)
   * - KCL(c)_cd
     - KCL(c)
     - Crystalline, DHS (irregular shape)
   * - Mg05Fe05SiO3(c)_am
     - Mg05Fe05SiO3(c)
     - Amorphous, Mie scattering (spherical)
   * - Mg05Fe05SiO3(c)_ad
     - Mg05Fe05SiO3(c)
     - Amorphous, DHS (irregular shape)
   * - Mg2SiO4(c)_am
     - Mg2SiO4(c)
     - Amorphous, Mie scattering (spherical)
   * - Mg2SiO4(c)_ad
     - Mg2SiO4(c)
     - Amorphous, DHS (irregular shape)
   * - Mg2SiO4(c)_cm
     - Mg2SiO4(c)
     - Crystalline, Mie scattering (spherical)
   * - Mg2SiO4(c)_cd
     - Mg2SiO4(c)
     - Crystalline, DHS (irregular shape)
   * - MgAl2O4(c)_cm
     - MgAl2O4(c)
     - Crystalline, Mie scattering (spherical)
   * - MgAl2O4(c)_cd
     - MgAl2O4(c)
     - Crystalline, DHS (irregular shape)
   * - MgFeSiO4(c)_am
     - MgFeSiO4(c)
     - Amorphous, Mie scattering (spherical)
   * - MgFeSiO4(c)_ad
     - MgFeSiO4(c)
     - Amorphous, DHS (irregular shape)
   * - MgSiO3(c)_am
     - MgSiO3(c)
     - Amorphous, Mie scattering (spherical)
   * - MgSiO3(c)_ad
     - MgSiO3(c)
     - Amorphous, DHS (irregular shape)
   * - MgSiO3(c)_cm
     - MgSiO3(c)
     - Crystalline, Mie scattering (spherical)
   * - MgSiO3(c)_cd
     - MgSiO3(c)
     - Crystalline, DHS (irregular shape)
   * - Na2S(c)_cm
     - Na2S(c)
     - Crystalline, Mie scattering (spherical)
   * - Na2S(c)_cd
     - Na2S(c)
     - Crystalline, DHS (irregular shape)
   * - SiC(c)_cm
     - SiC(c)
     - Crystalline, Mie scattering (spherical)
   * - SiC(c)_cd
     - SiC(c)
     - Crystalline, DHS (irregular shape)
   
		 
Rayleigh scatterers
___________________

.. list-table::
   :widths: 10 10
   :header-rows: 1
		 
   * - Species name
     - Required in abundance dictionary
   * - H2
     - H2
   * - He
     - He
   * - H2O
     - H2O
   * - CO2
     - CO2
   * - O2
     - O2
   * - N2
     - N2

Continuum opacity sources
_________________________

.. list-table::
   :widths: 10 10 80
   :header-rows: 1
		 
   * - Species name
     - Required in abundance dictionary
     - Descripton
   * - H2-H2
     - H2
     - Collision induced absorption (CIA)
   * - H2-He
     - H2, He
     - Collision induced absorption (CIA)
   * - H-
     - H, H-, e-
     - H- bound-free and free-free opacity
