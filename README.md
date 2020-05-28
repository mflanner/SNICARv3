# SNICARv3

This code calculates snow spectral albedo and radiative fluxes at the
interfaces of non-homogenous snow layers using a two-stream approximation.

TO RUN THE CODE:
-----------------------------------------------------------------------
(1) Download the snicar_v3.m and snicar_v3_drv.m source files.
(2) Download one of the following optics library packages:
   - snicar_v3_optics.tar.gz (471 MB), which includes all optics files
     EXCEPT those of algae
   - snicar_v3_optics_with_algae.tar.gz (7 GB), which includes
     everything in the above package AND algae. Please note that the
     algae library contains > 125,000 files.
(3) Unpack your optics library into the directory of your choice,
    using: "tar xvfz <file>"
(4) In snicar_v3.m, set variable "dir_op_root" to the directory where
    you unpacked the optics library.
(5) Run the template driver routine "snicar_v3_drv.m" and see basic
    output fields and a plot of spectral albedo!

Please see the header of snicar_v3.m for more information about input
fields.


NEW FEATURES IN THIS RELEASE:
--------------------------------------------------------------------
- Snow algae (Cook et al.,2017)
- Non-spherical ice particles (He et al.,2017)
- CO2 ice (Hansen et al.,2005; Singh et al.,2016)
- Additional dust optical properties:
   - Saharan dust (Balkanski et al.,2007)
   - San Juan Mountains, Colorado (Skiles et al.,2017)
   - Greenland (Polashenski et al.,2015)
   - Martian dust (Wolff et al.2009,2010; Singh et al.,2016)
- Volcanic ash (Flanner et al.,2014)
- A larger size bin (5-50um radius) for dust and ash particles
- Multiple options for H2O ice refractive indices
   - Warren et al.(1984)
   - Warren and Brandt (2008)
   - Picard et al. (2016)
- Several new options for surface spectral irradiance weighting
- Extension of the simulated spectrum down to 0.2um
