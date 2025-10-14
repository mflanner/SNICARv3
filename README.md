# ❄️ SNICAR-ADv3 + Graphical User Interface

This code uses a two-stream approximation to calculate snow spectral albedo and radiative fluxes at the interfaces of any number of snow layers with distinct properties. Now with a GUI for offline modeling within the MATLAB environment.

## 🛠️ How to Use
- (1) Download the `snicarAD_v3.m` and `snicarAD_v3_drv.m` source files. **Optional**: Download the `SNICARv3.mlapp` file, which is a MATLAB graphical user interface (GUI).
- (2) Download one of the following optics library packages from "Links to optics libraries": (a) `snicar_v3_optics.tar.gz` (637 MB), which includes all optics files EXCEPT those of algae; or (b) `snicar_v3_optics_with_algae.tar.gz` (7.1 GB), which includes everything in the above package AND algae. Please note that the algae library contains > 125,000 files. 
- (3) Unpack the optics library into the directory of your choice, using: "`tar xvfz <file>`". If you are using Windows10, right-click the optics folder and use `7-Zip` to extract the contents of the `.tar` file.
- (4) In `snicarAD_v3.m`, set variable "`dir_op_root`" to the directory where you unpacked the optics library. Use the full name of the search path to ensure that the model runs properly. For example:

```matlab
<Ln184> dir_op_root = 'C:\Users\JohnDoe\Desktop\SNICAR\OpticsLibrary\snicar_480band\'; 
```

- (5) Run the template driver routine "`snicarAD_v3_drv.m`" and see basic output fields and a plot of spectral albedo.

Please see the header of `snicarAD_v3.m` for more information about input fields.

## ✨ New Features
### October 2025
- Introducing `SNICARv3.mlapp`: A graphical user interface that can be run offline in the MATLAB environment. Simply follow the instructions in the ["How to Use"](#how-to-Use) section and then type `SNICARv3` in the MATLAB command window to open the GUI.

### January 2021 
- Replacement of the tri-diagonal matrix two-stream solver with the Adding-Doubling solver (Dang et al.,2019)
- Solar zenith angle-dependent surface spectral irradiances
- Updated optical properties of light-absorbing constituents

### June 2020
- Snow algae (Cook et al.,2017)
- Non-spherical ice particles (He et al.,2017)
- CO2 ice (Hansen et al.,2005; Singh et al.,2016)
- Additional dust optical properties: (a) Saharan dust (Balkanski et al.,2007); (b) San Juan Mountains, Colorado (Skiles et al.,2017); (c) Greenland (Polashenski et al.,2015); (d) Martian dust (Wolff et al.2009,2010; Singh et al.,2016).
- Volcanic ash (Flanner et al.,2014)
- A larger size bin (5-50um radius) for dust and ash particles
- Multiple options for H2O ice refractive indices: (a) Warren et al.(1984); (b) Warren and Brandt (2008); (c) Merged Picard et al. (2016) / Warren and Brandt (2008).
- Several new options for surface spectral irradiance weighting
- Extension of the simulated spectrum down to 0.2um

### 🎓 How to Cite
Please cite both of the following in peer-review publications:

* Flanner, M. G., Arnheim, J. B., Cook, J. M., Dang, C., He, C., Huang, X., Singh, D., Skiles, S. M., Whicker, C. A., & Zender, C. S. (2021). SNICAR-ADv3: A community tool for modeling spectral snow albedo. *Geoscientific Model Development*, *14*(12), 7673–7704. https://doi.org/10.5194/gmd-14-7673-2021
* Flanner, M. G. (2023). SNICAR-ADv3 [MATLAB]. GitHub. https://github.com/mflanner/SNICARv3

For `BibLaTeX`:
```tex
@article{flanner2021,
  title = {{{SNICAR-ADv3}}: {{A}} Community Tool for Modeling Spectral Snow Albedo},
  shorttitle = {{{SNICAR-ADv3}}},
  author = {Flanner, Mark G. and Arnheim, Julian B. and Cook, Joseph M. and Dang, Cheng and He, Cenlin and Huang, Xianglei and Singh, Deepak and Skiles, S. McKenzie and Whicker, Chloe A. and Zender, Charles S.},
  date = {2021-12-21},
  journaltitle = {Geoscientific Model Development},
  volume = {14},
  number = {12},
  pages = {7673--7704},
  doi = {10.5194/gmd-14-7673-2021}
}

@software{flanner2023,
  title = {{{SNICAR-ADv3}}},
  author = {Flanner, Mark G.},
  date = {2023-08-10},
  url = {https://github.com/mflanner/SNICARv3}
}
```

