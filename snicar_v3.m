% The Snow, Ice, and Aerosol Radiative (SNICAR) Model, version 3.0
%
% Author: Mark Flanner (flanner@umich.edu)
% Non-spherical ice representation written by Cenlin He (cenlinhe@ucar.edu)
% References cited and open-source license are at end of file.
%
% This routine calculates snow spectral albedo and radiative fluxes at
% the interfaces of non-homogeneous snow layers using a two-stream
% approximation (Toon et al.,1989). The original SNICAR model is
% described by Flanner et al. (2005,2007).  Incorporation of SNICAR
% into the Community Land Model is described by Lawrence et
% al. (2018). A new version of SNICAR that employs the adding-doubling
% solver (SNICAR-AD) is described by Dang et al. (2019), and we plan
% to ultimately merge features of this release (v3) with SNICAR-AD.
%
% This function reads particle optical properties from an external
% library of NetCDF files, based on user-defined properties of each
% snow layer.  The root directory containing this collection of files
% must be specified in variable "dir_op_root".
%
% Input parameters are described below. The user can run this code as
% a script instead of a function by commenting out the function line
% and setting '1==1' below the function definition.

%%%%%%%%%    NEW FEATURES IN THIS RELEASE   %%%%%%%%%%%%
%  - Snow algae (Cook et al.,2017)
%  - Non-spherical ice particles (He et al.,2017)
%  - CO2 ice (Hansen et al.,2005; Singh et al.,2016)
%  - Additional dust optical properties:
%     - Saharan dust (Balkanski et al.,2007)
%     - San Juan Mountains, Colorado (Skiles et al.,2017)
%     - Greenland (Polashenski et al.,2015)
%     - Martian dust (Wolff et al.2009,2010; Singh et al.,2016)
%  - Volcanic ash (Flanner et al.,2014)
%  - A larger size bin (5-50um radius) for dust and ash particles
%  - Multiple options for H2O ice refractive indices
%    - Warren et al.(1984)
%    - Warren and Brandt (2008)
%    - Picard et al. (2016)
%  - Several new options for surface spectral irradiance weighting
%  - Extension of the simulated spectrum down to 0.2um



%%%%%%%%%%  Input parameters: %%%%%%%%%%%
%
% direct_beam:   Direct or diffuse incident radiation
%                  1 = direct-beam (clear-sky)
%                  0 = diffuse (cloudy)
%
% atm:           Atmospheric profile used to obtain surface-incident spectral flux distribution
%                and subsequent broadband albedo:
%                  1 = mid-latitude winter
%                  2 = mid-latitude summer
%                  3 = sub-Arctic winter
%                  4 = sub-Arctic summer
%                  5 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
%                  6 = High Mountain (summer, surface pressure of 556 hPa)
%                  7 = Top-of-atmosphere
%
% flx_dwn_bb:    Broadband surface-incident solar flux [W/m2]. 
%                  Spectral flux fractions, determined by the
%                  atmospheric profile and direct_beam/cloudy flag, 
%                  are scaled by this constant

% coszen:        Cosine of solar zenith angle (only applies when direct_beam=1)
%                  Range: 0-1
%
% ice_ri:        Flag for ice refractive index data to use
%                  1 = Warren et al (1984)
%                  2 = Warren and Brandt (2008)
%                  3 = Picard et al (2016)
%                  4 = CO2 ice, Hansen et al (2005)
%
% R_sfc_all_wvl  Broadband albedo of underlying surface 
%                  Range: 0-1
%                  (user can also define a spectrally-dependent albedo below)
%
% dz:            Array of snow layer thicknesses [meters]. Top layer has index 1. 
%                  *** THE LENGTH OF THIS ARRAY DEFINES THE NUMBER OF SNOW LAYERS ***
%
% rho_snw:       Array of snow layer densities [kg m-3].
%                  Must have same length as dz
%
% rds_snw:       Array of snow layer effective grain radii [microns]
%                  Must have same length as dz
%
% options for asperical ice particles added by Cenlin He:
% sno_shp:      Snow grain shape
%               1=sphere (original SNICAR scheme)
%               2=spheroid; 3=hexagonal plate; 4=koch snowflake
% sno_fs:       Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
%               0=use recommended default value (He et al., 2017);
%               others(0<fs<1)= use user-specified value
%               only activated when sno_shp > 1 (i.e. nonspherical)
% sno_ar:       Aspect ratio: ratio of grain width to length
%               0=use recommended default value (He et al., 2017);
%               others(0.1<fs<20)= use user-specified value
%               only activated when sno_shp > 1 (i.e. nonspherical)
%
% mss_cnc_sot1:  Layer-specific array of mass mixing ratio of black carbon species 1 (uncoated BC)
%                  Units of parts per billion, ng g-1
%                  NOTE: This is mass of impurity per mass of snow
%                  (i.e., mass of impurity / mass of ice+impurity)
%
% mss_cnc_sot2:  Layer-specific array of mass mixing ratio of black carbon species 2 (sulfate-coated)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_brc1:  Layer-specific array of mass mixing ratio of brown carbon species 1
%                  Units of parts per billion, ng g-1
%
% mss_cnc_brc2:  Layer-specific array of mass mixing ratio of brown carbon species 2 (sulfate-coated)
%                  Units of parts per billion, ng g-1
%
% dust_type:     Type of dust to use:
%                  1 = Saharan dust (Balkanski et al., 2007, central hematite)
%                  2 = San Juan Mountains, CO (Skiles et al, 2017)
%                  3 = Greenland (Polashenski et al., 2015, central absorptivity)
%                  4 = Martian dust (Wolff et al., 2009;2010, Singh et al., 2016)
%
% mss_cnc_dst1:  Layer-specific array of mass mixing ratio of dust species 1 (radii of 0.05-0.5um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_dst2:  Layer-specific array of mass mixing ratio of dust species 2 (radii of 0.5-1.25um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_dst3:  Layer-specific array of mass mixing ratio of dust species 3 (radii of 1.25-2.5um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_dst4:  Layer-specific array of mass mixing ratio of dust species 4 (radii of 2.5-5.0um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_dst5:  Layer-specific array of mass mixing ratio of dust species 5 (radii of 5.0-50um)
%                  Units of parts per billion, ng g-1
%
% ash_type:      Type of volcanic ash to use
%                  1 = Eyjafjallajokull ash (Flanner et al.,2014)
%
% mss_cnc_ash1:  Layer-specific array of mass mixing ratio of volcanic ash species 1 (radii of 0.05-0.5um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_ash2:  Layer-specific array of mass mixing ratio of volcanic ash species 2 (radii of 0.5-1.25um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_ash3:  Layer-specific array of mass mixing ratio of volcanic ash species 3 (radii of 1.25-2.5um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_ash4:  Layer-specific array of mass mixing ratio of volcanic ash species 4 (radii of 2.5-5.0um)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_ash5:  Layer-specific array of mass mixing ratio of volcanic ash species 5 (radii of 5.0-50um)
%                  Units of parts per billion, ng g-1
%
% cell_nbr_conc: Layer-specific array of snow algal cell abundance
%                  Units of cells/mL or cells/gH20
%
% alg_rds:       Layer-specific mean cell radius of Gaussian distribution [um]
%                  Valid values are: [1 2 5 10 15 20 25 30 40 50] um
%
% dcmf_pig_chla: Dry cell mass fraction of chlorophyll-a
%                  Valid values are: [0.000 0.005 0.010 0.015 0.020 0.025 0.030]
%
% dcmf_pig_chlb: Dry cell mass fraction of chlorophyll-b
%                  Valid values are: [0.000 0.005 0.010 0.015 0.020 0.025 0.030]
%
% dcmf_pig_cara: Dry cell mass fraction of photoprotective carotenoids
%                  Valid values are: 0.00-0.15 by 0.01 increments
%
% dcmf_pig_carb: Dry cell mass fraction of photosynthetic carotenoids  
%                  Valid values are: 0.00-0.15 by 0.01 increments
%
%
%
%%%%%  Output data: %%%%%
%
% All output data is contained in the structure "data_out",
% described at end of this file


%%%%%%%%%%%%% BEGIN CODE: %%%%%%%%%%%%%%%

function data_out = snicar_v3(input_args)
    

% ROOT DIRECTORY FOR ALL OPTICAL AND SPECTRAL IRRADIANCE FILES
dir_op_root = '/data/flanner/mie/snicar_480band/';

% IF USER PREFERS TO RUN THIS AS A SCRIPT INSTEAD OF A FUNCTION:
%   (1) COMMENT OUT THE FUNCTION CALL ABOVE
%   (2) SET 1==1 BELOW AND DEFINE ALL INPUT IN THIS BLOCK

if (1==0)
    
    clear;
    
    % RADIATIVE TRANSFER CONFIGURATION:
    direct_beam   = 1;   % 1= Direct-beam incident flux, 0= Diffuse incident flux
                         % NOTE that cloudy-sky spectral fluxes are loaded when direct_beam=0
    
    % ATMOSPHERIC PROFILE for surface-incident flux:
    %    1 = mid-latitude winter
    %    2 = mid-latitude summer
    %    3 = sub-Arctic winter
    %    4 = sub-Arctic summer
    %    5 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
    %    6 = High Mountain (summer, surface pressure of 556 hPa)
    %    7 = Top-of-atmosphere
    % NOTE that clear-sky spectral fluxes are loaded when direct_beam=1,
    % and cloudy-sky spectral fluxes are loaded when direct_beam=0
    atm = 1;

    % Broadband surface-incident solar flux [W/m2]:
    %  (used to scale spectral flux fractions)
    flx_dwn_bb = 1.0;

    % COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM RT
    coszen = 0.5;
  
    % ICE REFRACTIVE INDEX DATASET TO USE:
    ice_ri = 3;
    
    % REFLECTANCE OF SURFACE UNDERLYING SNOW:
    % (value applied to all wavelengths.  user can also specify
    % spectrally-dependent ground albedo below)
    R_sfc_all_wvl = 0.25;

    % SNOW LAYER THICKNESSES [m]:
    %dz = [0.02 0.02 0.05];
    dz = [1000];
    nbr_lyr = length(dz);  % number of snow layers
  
    % SNOW DENSITY FOR EACH LAYER (units: kg/m3)
    rho_snw(1:nbr_lyr) = 150;

    % SNOW GRAIN SIZE FOR EACH LAYER (units: microns):
    rds_snw(1:nbr_lyr) = 1000;
  
    % Options added by Cenlin He for nonspherical ice particles based on
    % the parameterizations described by He et al. (2017,
    % doi:10.1175/JCLI-D-17-0300.1)
    
    sno_shp(1:nbr_lyr)  = 2;    % Snow grain shape option
                                % 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake

    sno_fs(1:nbr_lyr)   = 0;    % Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                % 0=use recommended default value (He et al. 2017);
                                % others(0<fs<1)= use user-specified value
                                % only activated when sno_shp > 1 (i.e. nonspherical)
                                
    sno_ar(1:nbr_lyr)   = 0;    % Aspect ratio: ratio of grain width to length
                                % 0=use recommended default value (He et al. 2017);
                                % others(0.1<fs<20)= use user-specified value
                                % only activated when sno_shp > 1 (i.e. nonspherical)
    
    % type of dust:
    dust_type = 1;              % 1=Saharan, 2=Colorado, 3=Greenland, 4=Mars
    
    % type of volcanic ash:
    ash_type = 1;               % 1 = Eyjafjallajokull
    
    % PARTICLE MASS MIXING RATIOS (units: ng g-1)
    % NOTE: This is mass of impurity per mass of snow
    %  (i.e., mass of impurity / mass of ice+impurity)
    mss_cnc_sot1(1:nbr_lyr)  = 0.0;    % uncoated black carbon [ng/g]
    mss_cnc_sot2(1:nbr_lyr)  = 0.0;    % coated black carbon [ng/g]
    mss_cnc_brc1(1:nbr_lyr)  = 0.0;    % uncoated brown carbon [ng/g]
    mss_cnc_brc2(1:nbr_lyr)  = 0.0;    % coated brown carbon [ng/g]
    mss_cnc_dst1(1:nbr_lyr)  = 0.0;    % dust species 1 [ng/g]
    mss_cnc_dst2(1:nbr_lyr)  = 0.0;    % dust species 2 [ng/g]
    mss_cnc_dst3(1:nbr_lyr)  = 0.0;    % dust species 3 [ng/g]
    mss_cnc_dst4(1:nbr_lyr)  = 0.0;    % dust species 4 [ng/g]
    mss_cnc_dst5(1:nbr_lyr)  = 0.0;    % dust species 5 [ng/g]
    mss_cnc_ash1(1:nbr_lyr)  = 0.0;    % volcanic ash species 1 [ng/g]
    mss_cnc_ash2(1:nbr_lyr)  = 0.0;    % volcanic ash species 2 [ng/g]
    mss_cnc_ash3(1:nbr_lyr)  = 0.0;    % volcanic ash species 3 [ng/g]
    mss_cnc_ash4(1:nbr_lyr)  = 0.0;    % volcanic ash species 4 [ng/g]
    mss_cnc_ash5(1:nbr_lyr)  = 0.0;    % volcanic ash species 5 [ng/g]

    cell_nbr_conc(1:nbr_lyr) = 0.0;    % algae [cells/mL]
    alg_rds(1:nbr_lyr)       = 10;     % mean cell radius (um)
    dcmf_pig_chla(1:nbr_lyr) = 0.02;   % dry cell mass fraction of chlorophyll-a
    dcmf_pig_chlb(1:nbr_lyr) = 0.02;   % dry cell mass fraction of chlorophyll-b
    dcmf_pig_cara(1:nbr_lyr) = 0.05;   % dry cell mass fraction of photoprotective_carotenoids
    dcmf_pig_carb(1:nbr_lyr) = 0.00;   % dry cell mass fraction of photosynthetic_carotenoids    
    
else
    % assign all input variables from function input structure:
    
    direct_beam   = input_args.direct_beam;     % direct beam or diffuse
    atm           = input_args.atm;             % atmospheric profile
    flx_dwn_bb    = input_args.flx_dwn_bb;      % broadband surface insolation
    coszen        = input_args.coszen;          % cosine of solar zenith angle
    ice_ri        = input_args.ice_ri;          % ice refractive index data
    R_sfc_all_wvl = input_args.R_sfc_all_wvl;   % albedo of underlying surface
    dz            = input_args.dz;              % snow layer thicknesses
    rho_snw       = input_args.rho_snw;         % snow layer densities
    rds_snw       = input_args.rds_snw;         % snow layer grain radii
    sno_shp       = input_args.sno_shp;         % snow layer grain shape
    sno_fs        = input_args.sno_fs;          % snow layer grain shape factor 
    sno_ar        = input_args.sno_ar;          % snow layer grain aspect ratio
    dust_type     = input_args.dust_type;       % type of dust
    ash_type      = input_args.ash_type;        % type of ash
    mss_cnc_sot1  = input_args.mss_cnc_sot1;    % uncoated black carbon [ng/g]
    mss_cnc_sot2  = input_args.mss_cnc_sot2;    % coated black carbon [ng/g]
    mss_cnc_brc1  = input_args.mss_cnc_brc1;    % uncoated brown carbon [ng/g]
    mss_cnc_brc2  = input_args.mss_cnc_brc2;    % coated brown carbon [ng/g]
    mss_cnc_dst1  = input_args.mss_cnc_dst1;    % dust species 1 [ng/g]
    mss_cnc_dst2  = input_args.mss_cnc_dst2;    % dust species 2 [ng/g]
    mss_cnc_dst3  = input_args.mss_cnc_dst3;    % dust species 3 [ng/g]
    mss_cnc_dst4  = input_args.mss_cnc_dst4;    % dust species 4 [ng/g]
    mss_cnc_dst5  = input_args.mss_cnc_dst5;    % dust species 5 [ng/g]
    mss_cnc_ash1  = input_args.mss_cnc_ash1;    % volcanic ash species 1 [ng/g]
    mss_cnc_ash2  = input_args.mss_cnc_ash2;    % volcanic ash species 2 [ng/g]
    mss_cnc_ash3  = input_args.mss_cnc_ash3;    % volcanic ash species 3 [ng/g]
    mss_cnc_ash4  = input_args.mss_cnc_ash4;    % volcanic ash species 4 [ng/g]
    mss_cnc_ash5  = input_args.mss_cnc_ash5;    % volcanic ash species 5 [ng/g]
    cell_nbr_conc = input_args.cell_nbr_conc;   % algae [cells/mL]
    alg_rds       = input_args.alg_rds;         % mean cell radius (um)
    dcmf_pig_chla = input_args.dcmf_pig_chla;   % dry cell mass fraction of chlorophyll-a
    dcmf_pig_chlb = input_args.dcmf_pig_chlb;   % dry cell mass fraction of chlorophyll-b
    dcmf_pig_cara = input_args.dcmf_pig_cara;   % dry cell mass fraction of photoprotective carotenoids
    dcmf_pig_carb = input_args.dcmf_pig_carb;   % dry cell mass fraction of photosynthetic carotenoids
end;


% Sub-directories of NetCDF files for (1) optical properties of
% light-absorbing impurities, (2) optical properties of snow algae,
% and (3) surface spectral irradiance profiles:
dir_lai     = strcat(dir_op_root,'lai/');
dir_alg     = strcat(dir_op_root,'alg_pig/');
dir_spc     = strcat(dir_op_root,'fsds/');


% Flag to apply the Delta approximation (Joseph,1976) 
%   1 = Use DELTA approximation (strongly recommended)
%   0 = don't use DELTA approximation
delta_aprx    = 1;


% Two-Stream Approximation Type:
%   1 = Eddington
%   2 = Quadrature
%   3 = Hemispheric Mean  
% NOTE: Delta-Eddington Approximation is likely best for the visible
% spectrum, but can provide negative albedo in the near-IR under
% diffuse light conditions. Hence, the Hemispheric Mean approximation
% is recommended for general use.
aprx_typ      = 3;


% Set wavelength grid (um) and wavelength number:
wvl     = [0.205:0.01:4.995];
nbr_wvl = length(wvl);


% file substrings for ice Mie parameters:
if (ice_ri == 1)
    stb1 = 'ice_Wrn84';
elseif (ice_ri == 2)
    stb1 = 'ice_Wrn08';
elseif (ice_ri == 3)
    stb1 = 'ice_Pic16';
elseif (ice_ri == 4)
    stb1 = 'co2ice';
end;

% subdirectory for ice optical properties: 
dir_ice  = strcat(dir_op_root,stb1,'/');


% filenames for impurity optical properties:
fl_sot1  = 'mie_sot_ChC90_dns_1317.nc';
fl_sot2  = 'miecot_slfsot_ChC90_dns_1317.nc';

fl_brc1  = 'brC_Kirch_BCsd.nc';
fl_brc2  = 'brC_Kirch_BCsd_slfcot.nc';

if (dust_type==1)
    fl_dst1  = 'dust_balkanski_central_size1.nc';
    fl_dst2  = 'dust_balkanski_central_size2.nc';
    fl_dst3  = 'dust_balkanski_central_size3.nc';
    fl_dst4  = 'dust_balkanski_central_size4.nc';
    fl_dst5  = 'dust_balkanski_central_size5.nc';
elseif (dust_type==2)
    fl_dst1  = 'dust_skiles_size1.nc';
    fl_dst2  = 'dust_skiles_size2.nc';
    fl_dst3  = 'dust_skiles_size3.nc';
    fl_dst4  = 'dust_skiles_size4.nc';
    fl_dst5  = 'dust_skiles_size5.nc';
elseif (dust_type==3)
    fl_dst1  = 'dust_greenland_central_size1.nc';
    fl_dst2  = 'dust_greenland_central_size2.nc';
    fl_dst3  = 'dust_greenland_central_size3.nc';
    fl_dst4  = 'dust_greenland_central_size4.nc';
    fl_dst5  = 'dust_greenland_central_size5.nc';
elseif (dust_type==4)
    fl_dst1  = 'dust_mars_size1.nc';
    fl_dst2  = 'dust_mars_size2.nc';
    fl_dst3  = 'dust_mars_size3.nc';
    fl_dst4  = 'dust_mars_size4.nc';
    fl_dst5  = 'dust_mars_size5.nc';
end;

if (ash_type==1)
    fl_ash1  = 'volc_ash_eyja_central_size1.nc';
    fl_ash2  = 'volc_ash_eyja_central_size2.nc';
    fl_ash3  = 'volc_ash_eyja_central_size3.nc';
    fl_ash4  = 'volc_ash_eyja_central_size4.nc';
    fl_ash5  = 'volc_ash_eyja_central_size5.nc';
end;

% create cell structure of impurity file names
f1  = strcat(dir_lai,fl_sot1);
f2  = strcat(dir_lai,fl_sot2);
f3  = strcat(dir_lai,fl_brc1);
f4  = strcat(dir_lai,fl_brc2);
f5  = strcat(dir_lai,fl_dst1);
f6  = strcat(dir_lai,fl_dst2);
f7  = strcat(dir_lai,fl_dst3);
f8  = strcat(dir_lai,fl_dst4);
f9  = strcat(dir_lai,fl_dst5);
f10 = strcat(dir_lai,fl_ash1);
f11 = strcat(dir_lai,fl_ash2);
f12 = strcat(dir_lai,fl_ash3);
f13 = strcat(dir_lai,fl_ash4);
f14 = strcat(dir_lai,fl_ash5);

tmp_char = strvcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14);
files    = cellstr(tmp_char);


% NUMBER OF PARTICLE SPECIES IN SNOW (ICE AND ALGAE EXCLUDED)
tmp_sz  = size(files);
nbr_aer = tmp_sz(1);


% REFLECTANCE OF UNDERLYING SURFACE
% (OPTIONAL: User can set spectrally-dependent surface albedo here):
R_sfc(1:nbr_wvl,1) = R_sfc_all_wvl;

% files containing surface-incident spectral flux fractions:
stb_spc1 = 'swnb_480bnd_';
var_spc  = 'flx_frc_sfc';

if (atm==1)
    stb_atm = 'mlw';
elseif (atm==2)
    stb_atm = 'mls';
elseif (atm==3)
    stb_atm = 'saw';
elseif (atm==4)
    stb_atm = 'sas';
elseif (atm==5)
    stb_atm = 'smm';
elseif (atm==6)
    stb_atm = 'hmn';
elseif (atm==7)
    stb_atm = 'toa';
    var_spc = 'flx_frc_toa';
end;

%mu_not=cos((slr_znt/360)*2*pi);
mu_not = coszen;

% Set incident solar flux spectral distribution:
if (direct_beam == 1)
    fi_spc                    = strcat(dir_spc,stb_spc1,stb_atm,'_clr.nc');
    flx_slr                   = ncread(fi_spc,var_spc);    
    flx_slr(find(flx_slr==0)) = 1E-30;
    
    % direct-beam incident spectral flux [W/m2/band]
    Fs(1:nbr_wvl,1) = flx_dwn_bb.*flx_slr./(mu_not*pi);

    % diffuse incident spectral flux [W/m2/band]:
    Fd(1:nbr_wvl,1) = 0;                    
    
elseif (direct_beam == 0)
    fi_spc                    = strcat(dir_spc,stb_spc1,stb_atm,'_cld.nc');
    flx_slr                   = ncread(fi_spc,var_spc);    
    flx_slr(find(flx_slr==0)) = 1E-30;
    
    % direct-beam incident spectral flux [W/m2/band]
    Fd(1:nbr_wvl,1) = flx_dwn_bb.*flx_slr;

    % diffuse incident spectral flux [W/m2/band]
    Fs(1:nbr_wvl,1) = 0;
end;

% Set indices for visible (0.2-0.7um) and near-IR (0.7-5.0um) bands
vis_max_idx = 50;
nir_max_idx = length(wvl);

  
%%% Constants for aspherical ice particles %%%
%
% g_snw asymmetry factor parameterization coefficients (6 bands) from
% Table 3 & Eqs. 6-7 in He et al. (2017)
% assume same values for 4-5 um band, which leads to very small biases (<3%)
g_wvl = [0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00]; % wavelength (um) division point
g_wvl_center = g_wvl(2:8)/2 + g_wvl(1:7)/2 ; % center point for wavelength band
g_b0 = [9.76029E-01,9.67798E-01,1.00111E+00,1.00224E+00,9.64295E-01,9.97475E-01,9.97475E-01];
g_b1 = [5.21042E-01,4.96181E-01,1.83711E-01,1.37082E-01,5.50598E-02,8.48743E-02,8.48743E-02];
g_b2 = [-2.66792E-04,1.14088E-03,2.37011E-04,-2.35905E-04,8.40449E-04,-4.71484E-04,-4.71484E-04];
% Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007 
g_F07_c2 = [1.349959e-1,1.115697e-1,9.853958e-2,5.557793e-2,-1.233493e-1,0.0,0.0];
g_F07_c1 = [-3.987320e-1,-3.723287e-1,-3.924784e-1,-3.259404e-1,4.429054e-2,-1.726586e-1,-1.726586e-1];
g_F07_c0 = [7.938904e-1,8.030084e-1,8.513932e-1,8.692241e-1,7.085850e-1,6.412701e-1,6.412701e-1];
g_F07_p2 = [3.165543e-3,2.014810e-3,1.780838e-3,6.987734e-4,-1.882932e-2,-2.277872e-2,-2.277872e-2];
g_F07_p1 = [1.140557e-1,1.143152e-1,1.143814e-1,1.071238e-1,1.353873e-1,1.914431e-1,1.914431e-1];
g_F07_p0 = [5.292852e-1,5.425909e-1,5.601598e-1,6.023407e-1,6.473899e-1,4.634944e-1,4.634944e-1];
		

%%% Constants for snow algae %%%
sigma_alg = alg_rds.*0.1;  % standard deviation of Gaussian size distribution [um] (consistent with Mie)
rho_alg   = 1080;          % algae density [kg/m3] (consistent with Mie calculations)


nbr_lyr   = length(dz);      % number of snow layers

% Establish layer-specific filenames of Mie ice and algae optical properties:
for n=1:nbr_lyr

    % 1. ice:
    fl_ice = strcat(dir_ice,stb1,'_',sprintf('%04d',rds_snw(n)),'.nc');
        
    % 2. algae, based on cell size and pigment concentrations:
    if (cell_nbr_conc(n) > 0.0)
        fl_alg = strcat(dir_alg,...
                        'alg_sph_r',sprintf('%03d',round(alg_rds(n))),'um_',...
                        'chla',sprintf('%03d',round(dcmf_pig_chla(n)*1000)),'_',...
                        'chlb',sprintf('%03d',round(dcmf_pig_chlb(n)*1000)),'_',...
                        'cara',sprintf('%03d',round(dcmf_pig_cara(n)*1000)),'_',...
                        'carb',sprintf('%03d',round(dcmf_pig_carb(n)*1000)),...
                        '.nc');
    
        % check that derived alg-impurity file exists:
        if (exist(fl_alg,'file')~=2)
            error(strcat('Snow algae impurity file: ',fl_alg,' does not exist.'));
        end;
    end; 
    
    % single-scatter albedo, mass extinction coefficient, and
    % asymmetry paramater of ice:
    omega_ice(:,n)       = ncread(fl_ice,'ss_alb');
    ext_cff_mss_ice(:,n) = ncread(fl_ice,'ext_cff_mss');
    
    %%% shape-dependent asymmetry factors (Cenlin He) %%%%
    if (sno_shp(n) == 1)  % spheres from Mie properties (original SNICAR scheme)
        g_ice(:,n)       = ncread(fl_ice,'asm_prm'); % asymmetry paramater
        
    elseif (sno_shp(n) == 2) % 2=spheroid, He et al. (2017) parameterization
        diam_ice = 2.0 .* rds_snw(n); % effective snow grain diameter        
        if (sno_fs(n) == 0)
            fs_sphd = 0.929; % default shape factor for spheroid; He et al. (2017), Table 1
        else
            fs_sphd = sno_fs(n);
        end
        fs_hex = 0.788; % shape factor for hexagonal plate (reference);
        if (sno_ar(n) == 0)
            AR_tmp = 0.5; % default aspect ratio for spheroid; He et al. (2017), Table 1
        else
            AR_tmp = sno_ar(n);
        end
        g_ice_Cg_tmp = g_b0 .* (fs_sphd/fs_hex).^g_b1 .* diam_ice.^g_b2; % Eq.7, He et al. (2017)
        gg_ice_F07_tmp = g_F07_c0 + g_F07_c1 .* AR_tmp + g_F07_c2 .* AR_tmp^2;% Eqn. 3.1 in Fu (2007)
        
    elseif (sno_shp(n) == 3) % 3=hexagonal plate, He et al. 2017 parameterization
        diam_ice = 2.0 .* rds_snw(n); % effective snow grain diameter        
        if (sno_fs(n) == 0)
            fs_hex0 = 0.788; % default shape factor for hexagonal plates; He et al. (2017), Table 1
        else
            fs_hex0 = sno_fs(n);
        end
        fs_hex = 0.788; % shape factor for hexagonal plate (reference);
        if (sno_ar(n) == 0)
            AR_tmp = 2.5; % default aspect ratio for hexagonal plate; He et al. (2017), Table 1
        else
            AR_tmp = sno_ar(n);
        end
        g_ice_Cg_tmp = g_b0 .* (fs_hex0/fs_hex).^g_b1 .* diam_ice.^g_b2; % Eq.7, He et al. (2017)
        gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 .* log(AR_tmp) + g_F07_p2 .* (log(AR_tmp))^2;   % Eqn. 3.3 in Fu (2007) 
        
    elseif (sno_shp(n) == 4) % 4=koch snowflake, He et al. (2017) parameterization
        diam_ice = 2.0 .* rds_snw(n) ./ 0.544; % effective snow grain diameter
        if (sno_fs(n) == 0)
            fs_koch = 0.712; % default shape factor for koch snowflake; He et al. (2017), Table 1
        else
            fs_koch = sno_fs(n);
        end
        fs_hex = 0.788; % shape factor for hexagonal plate (reference);
        if (sno_ar(n) == 0)
            AR_tmp = 2.5; % default aspect ratio for koch snowflake; He et al. (2017), Table 1
        else
            AR_tmp = sno_ar(n);
        end
        g_ice_Cg_tmp = g_b0 .* (fs_koch/fs_hex).^g_b1 .* diam_ice.^g_b2; % Eq.7, He et al. (2017)
        gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 .* log(AR_tmp) + g_F07_p2 .* (log(AR_tmp))^2;   % Eqn. 3.3 in Fu (2007)
         
    end
    
    if (sno_shp(n) > 1)
        % 6 wavelength bands for g_ice to be interpolated into 480-bands of SNICAR
        % shape-preserving piecewise interpolation into 480-bands 
        g_Cg_intp = pchip(g_wvl_center,g_ice_Cg_tmp,wvl) ;
        gg_F07_intp = pchip(g_wvl_center,gg_ice_F07_tmp,wvl) ;
        g_ice_F07 = gg_F07_intp + (1.0 - gg_F07_intp) ./ omega_ice(:,n)' ./ 2; % Eq.2.2 in Fu (2007)
        g_ice(:,n) = g_ice_F07 .* g_Cg_intp; % Eq.6, He et al. (2017)
        g_ice(381:480,n) = g_ice(380); % assume same values for 4-5 um band, with very small biases (<3%)
    end;
    
    g_ice(g_ice > 0.99) = 0.99; % avoid unreasonable values (so far only occur in large-size spheroid cases)
    %%%%%%%%%%%%%%%%%
    
    
    % single-scatter albedo, mass extinction coefficient, and
    % asymmetry paramater of algae:
    if (cell_nbr_conc(n) > 0)
        omega_alg(:,n)       = ncread(fl_alg,'ss_alb');
        ext_cff_mss_alg(:,n) = ncread(fl_alg,'ext_cff_mss');
        g_alg(:,n)           = ncread(fl_alg,'asm_prm');
    end;
end


% Read Mie LAI parameters (layer-independent)

% read NetCDF properties
for j=1:nbr_aer
    fl_in                = char(files(j));
    omega_aer(:,j)       = ncread(fl_in,'ss_alb');
    g_aer(:,j)           = ncread(fl_in,'asm_prm');
    if ((j==2) | (j==4))
        % unique variable name for mass extinction coefficient of
        % coated aerosols (i.e., coated BC and coated BrC)
        ext_cff_mss_aer(:,j) = ncread(fl_in,'ext_cff_mss_ncl');
    else
        ext_cff_mss_aer(:,j) = ncread(fl_in,'ext_cff_mss');
    end;
end

% Set aerosol concentration matrix:
mss_cnc_aer(1:nbr_lyr,1)  = mss_cnc_sot1;
mss_cnc_aer(1:nbr_lyr,2)  = mss_cnc_sot2;
mss_cnc_aer(1:nbr_lyr,3)  = mss_cnc_brc1;
mss_cnc_aer(1:nbr_lyr,4)  = mss_cnc_brc2;

mss_cnc_aer(1:nbr_lyr,5)  = mss_cnc_dst1;
mss_cnc_aer(1:nbr_lyr,6)  = mss_cnc_dst2;
mss_cnc_aer(1:nbr_lyr,7)  = mss_cnc_dst3;
mss_cnc_aer(1:nbr_lyr,8)  = mss_cnc_dst4;
mss_cnc_aer(1:nbr_lyr,9)  = mss_cnc_dst5;

mss_cnc_aer(1:nbr_lyr,10)  = mss_cnc_ash1;
mss_cnc_aer(1:nbr_lyr,11)  = mss_cnc_ash2;
mss_cnc_aer(1:nbr_lyr,12)  = mss_cnc_ash3;
mss_cnc_aer(1:nbr_lyr,13)  = mss_cnc_ash4;
mss_cnc_aer(1:nbr_lyr,14)  = mss_cnc_ash4;


% convert to units of kg/kg:
mss_cnc_aer  = mss_cnc_aer.*10^-9;

if (cell_nbr_conc(n) > 0)
    % mean algal cell volume (3rd moment of Gaussian distribution) by layer:
    mean_vol_cell = 4/3*pi .* (alg_rds.^3 + 3.*alg_rds.*sigma_alg.^2); %[um^3/cell]

    % mean mass per cell by layer:
    mass_per_cell = mean_vol_cell.*1E-18.*rho_alg; % [kg/cell]
    
    % mass concentration of algae by layer
    mss_cnc_alg   = cell_nbr_conc.*1000.*mass_per_cell; % [kg/kg]
end;


% BEGIN RT SOLVER:

% Calculate effective tau, omega, g for the (ice+algae+impurity) system
for n=1:nbr_lyr
    
    % Snow column mass [kg/m^2] (array)
    % Mass of snow is ice+impurities
    L_snw(n)     = rho_snw(n)*dz(n);
    
    % burden and optical thickness of algae
    if (cell_nbr_conc(n) > 0)
        L_alg(n)     = L_snw(n)*mss_cnc_alg(n);
        tau_alg(:,n) = L_alg(n).*ext_cff_mss_alg(:,n);
    else
        L_alg(n)     = 0.0;
        tau_alg(:,n) = 0.0;
    end;
    
    % burdens and optical thicknesses of LAI
    for j=1:nbr_aer
        L_aer(n,j)     = L_snw(n)*mss_cnc_aer(n,j);
        tau_aer(:,n,j) = L_aer(n,j).*ext_cff_mss_aer(:,j);
    end
 
    % ice mass = snow mass - impurity mass (generally a tiny correction)
    L_ice(n)     = L_snw(n) - L_alg(n) - sum(L_aer(n,:));
    
    if (L_ice(n) < 0)
        error(['Impurity load cannot exceed snow load. Snow mass ' ...
               'is assumed to be that of ice+impurities, so the sum ' ...
               'of impurity mixing ratios cannot exceed 1']);
    end;

    % optical thickness due to ice:
    tau_ice(:,n) = L_ice(n).*ext_cff_mss_ice(:,n);
    
    tau_sum(1:nbr_wvl,1)   = 0.0;
    omega_sum(1:nbr_wvl,1) = 0.0;
    g_sum(1:nbr_wvl,1)     = 0.0;
  
    for j=1:nbr_aer
        tau_sum   = tau_sum + tau_aer(:,n,j);
        omega_sum = omega_sum + (tau_aer(:,n,j).*omega_aer(:,j));
        g_sum     = g_sum + (tau_aer(:,n,j).*omega_aer(:,j).*g_aer(:,j));
    end
  
    % add algae contribution to weighted LAI sums:
    if (cell_nbr_conc(n) > 0)
        tau_sum   = tau_sum + tau_alg(:,n);
        omega_sum = omega_sum + (tau_alg(:,n).*omega_alg(:,n));
        g_sum     = g_sum + (tau_alg(:,n).*omega_alg(:,n).*g_alg(:,n));
    end;
    
    % cleaner: add ice contribution to weighted LAI+alg sums:
    %tau_sum   = tau_sum + tau_ice(:,n);
    %omega_sum = omega_sum + (tau_ice(:,n).*omega_ice(:,n));
    %g_sum     = g_sum + (tau_ice(:,n).*omega_ice(:,n).*g_ice(:,n));

    % cleaner:
    %tau(:,n)   = tau_sum;
    %omega(:,n) = (1./tau(:,n)).*omega_sum;
    %g(:,n)     = (1./(tau(:,n).*omega(:,n))).*g_sum;

    % original: weighted sums, including ice:
    tau(:,n)   = tau_sum + tau_ice(:,n);
    omega(:,n) = (1./tau(:,n)).*(omega_sum+ (omega_ice(:,n).*tau_ice(:,n)));
    g(:,n)     = (1./(tau(:,n).*omega(:,n))) .* (g_sum+ (g_ice(:,n).*omega_ice(:,n).*tau_ice(:,n)));
end
  
  
% Perform Delta-transformations, if called for
if (delta_aprx == 1)
    g_star     = g./(1+g);
    omega_star = ((1-(g.^2)).*omega) ./ (1-(omega.*(g.^2)));
    tau_star   = (1-(omega.*(g.^2))).*tau;
else
    g_star     = g;
    omega_star = omega;
    tau_star   = tau;
end;

% Calculate total column optical depth
% tau_clm(:,n) = total optical depth of layers above n, or optical
% depth from upper model boundary to top of layer n
tau_clm(1:nbr_wvl,1:1) = 0.0;
for n=2:nbr_lyr
    tau_clm(:,n) = tau_clm(:,n-1)+tau_star(:,n-1);
end


% Boundary Condition: Upward (reflected) direct beam flux at lower model boundary
% (Eq. 37)
S_sfc = R_sfc.*mu_not.*exp(-(tau_clm(:,nbr_lyr)+tau_star(:,nbr_lyr))./mu_not).*pi.*Fs;


% Apply 2-stream approximation technique (Toon et al., Table 1.)

if (aprx_typ == 1)
    % Eddington:
    gamma1 = (7-(omega_star.*(4+(3*g_star))))./4;
    gamma2 = -(1-(omega_star.*(4-(3*g_star))))./4;
    gamma3 = (2-(3*g_star.*mu_not))./4;
    gamma4 = 1-gamma3;
    mu_one = 0.5;
  
elseif (aprx_typ == 2)
    % Quadrature: 
    gamma1 = sqrt(3).*(2-(omega_star.*(1+g_star)))./2;
    gamma2 = omega_star.*sqrt(3).*(1-g_star)./2;
    gamma3 = (1-(sqrt(3).*g_star.*mu_not))./2;
    gamma4 = 1-gamma3;
    mu_one = 1/sqrt(3);
  
elseif (aprx_typ == 3)
    % Hemispheric mean:
    gamma1 = 2 - (omega_star.*(1+g_star));
    gamma2 = omega_star.*(1-g_star);
    gamma3 = (1-(sqrt(3).*g_star.*mu_not))./2;
    gamma4 = 1-gamma3;
    mu_one = 0.5;
end;

% Eq. 21 and 22
lambda = sqrt(abs((gamma1.^2) - (gamma2.^2)));
GAMMA  = gamma2./(gamma1+lambda);

% Eq. 44
e1 = 1+(GAMMA.*exp(-lambda.*tau_star));
e2 = 1-(GAMMA.*exp(-lambda.*tau_star));
e3 = GAMMA + exp(-lambda.*tau_star);
e4 = GAMMA - exp(-lambda.*tau_star);


% Before calculating the C functions:
%
% C functions are indeterminate if [lambda^2 = 1/(mu_not^2)].
% Upper bound of lambda^2 = 4 for delta-hemispheric mean, and
% lambda^2 = 3 for delta-Eddington.
%
% So, problem can only arise when 0.5 < mu_not < 1.0.  
%
% Assuming that abs(lambda^2 - 1/(mu_not^2)) < 0.01 is dangerous:
%  Let x= lambda^2 - 1/(mu_not^2),
%  dx/dmu_not = 2*mu_not^-3
%  dx = 0.01 = 2*(mu_not^-3)*dmu_not
% For the range of possible mu_not:
%  0.04 > dmu_not > 0.005
% So, changing mu_not by 0.04 would be sufficient when this condition arises.
%
% This adjustment is NOT implemented here because of vectorized
% code. This adjustment MUST be implemented in ESMs like CLM


% Calculate C function for each layer
% (evaluated at the top and bottom of each layer n)
% (Eq. 23 and 24)

for n=1:nbr_lyr

    if (sum(Fs) > 0.0)
      
        % Stability check, see explanation above:
        %for p=1:nbr_wvl
        %temp1 = ( lambda(p,n)^2) - 1/mu_not^2);
        %if abs(temp) < 0.01
        %  mu_not_sv = mu_not;
        %  if (temp1 > 0)
        %    mu_not = mu_not + 0.04;
        %  else
        %    mu_not = mu_not - 0.04
        %  end;
        %end;
        %end
        %if (n==1)
        %    mu_not=mu_not+0.00;
        %else
        %    mu_not = mu_not-0.00;
        %end;
	
        C_pls_btm(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not).*...
                          (((gamma1(:,n)-(1/mu_not)).*gamma3(:,n))+ ...
                           (gamma4(:,n).*gamma2(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
  
            
        C_mns_btm(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not).*...
                          (((gamma1(:,n)+(1/mu_not)).*gamma4(:,n))+ ...
                           (gamma2(:,n).*gamma3(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));


        C_pls_top(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-tau_clm(:,n)./mu_not).*...
                          (((gamma1(:,n)-(1/mu_not)).*gamma3(:,n))+(gamma4(:,n).*gamma2(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
        
            
        C_mns_top(:,n) = (omega_star(:,n).*pi.*Fs.*...
                          exp(-tau_clm(:,n)./mu_not).*...
                          (((gamma1(:,n)+(1/mu_not)).*gamma4(:,n))+(gamma2(:,n).*gamma3(:,n))))./...
            ((lambda(:,n).^2)-(1/(mu_not^2)));
        
    else
        % no direct-beam flux:
        C_pls_btm(1:nbr_wvl,n) = 0.0;
        C_mns_btm(1:nbr_wvl,n) = 0.0;
        C_pls_top(1:nbr_wvl,n) = 0.0;
        C_mns_top(1:nbr_wvl,n) = 0.0;
    end;
      
end
  

% Eq. 41-43
for i=1:2*nbr_lyr
    % Boundary values for i=1 and i=2nbr_lyr, specifics for i=odd and i=even    
    if (i==1)
        A(1:nbr_wvl,1) = 0.0;
        B(1:nbr_wvl,1) = e1(:,1);
        D(1:nbr_wvl,1) = -e2(:,1);
        E(1:nbr_wvl,1) = Fd(:)-C_mns_top(:,1);
					  
    elseif(i==2*nbr_lyr)
        A(:,i) = e1(:,nbr_lyr)-(R_sfc.*e3(:,nbr_lyr));
        B(:,i) = e2(:,nbr_lyr)-(R_sfc.*e4(:,nbr_lyr));
        D(:,i) = 0.0;
        E(:,i) = S_sfc(:) - C_pls_btm(:,nbr_lyr) + (R_sfc.*C_mns_btm(:,nbr_lyr));
    
    elseif(mod(i,2)==1)  % If odd and i>=3 (n=1 for i=3. this is confusing)
        n      = floor(i/2);
        A(:,i) = (e2(:,n).*e3(:,n))-(e4(:,n).*e1(:,n));
        B(:,i) = (e1(:,n).*e1(:,n+1))-(e3(:,n).*e3(:,n+1));
        D(:,i) = (e3(:,n).*e4(:,n+1))-(e1(:,n).*e2(:,n+1));
        E(:,i) = (e3(:,n).*(C_pls_top(:,n+1)-C_pls_btm(:,n))) + ...
                 (e1(:,n).*(C_mns_btm(:,n)-C_mns_top(:,n+1)));
    
    elseif(mod(i,2)==0)  % If even and i<=2nbr_lyr
        n      = (i/2);
        A(:,i) = (e2(:,n+1).*e1(:,n))-(e3(:,n).*e4(:,n+1));
        B(:,i) = (e2(:,n).*e2(:,n+1))-(e4(:,n).*e4(:,n+1));
        D(:,i) = (e1(:,n+1).*e4(:,n+1))-(e2(:,n+1).*e3(:,n+1));
        E(:,i) = (e2(:,n+1).*(C_pls_top(:,n+1)-C_pls_btm(:,n))) + ...
                 (e4(:,n+1).*(C_mns_top(:,n+1)-C_mns_btm(:,n))); 
    end;
end


% Eq. 45
AS(1:nbr_wvl,2*nbr_lyr) = A(:,2*nbr_lyr)./B(:,2*nbr_lyr);
DS(1:nbr_wvl,2*nbr_lyr) = E(:,2*nbr_lyr)./B(:,2*nbr_lyr);
  
% Eq. 46
for i=(2*nbr_lyr-1):-1:1
    X(1:nbr_wvl,i) = 1./(B(:,i)-(D(:,i).*AS(:,i+1)));
    AS(:,i)        = A(:,i).*X(:,i);
    DS(:,i)        = (E(:,i)-(D(:,i).*DS(:,i+1))).*X(:,i);
end

% Eq. 47
Y(1:nbr_wvl,1) = DS(:,1);
for i=2:2*nbr_lyr
    Y(:,i) = DS(:,i)-(AS(:,i).*Y(:,i-1));
end;


for n=1:nbr_lyr
    % Direct beam flux at the base of each layer (Eq. 50)
    direct(1:nbr_wvl,n) = mu_not*pi*Fs.*exp(-(tau_clm(:,n)+tau_star(:,n))./mu_not);

    % Net flux (positive upward = F_up-F_down) at the base of each
    % layer (Eq. 48)
    F_net(1:nbr_wvl,n) = (Y(:,(2*n-1)).*(e1(:,n)-e3(:,n))) +...
        (Y(:,(2*n)).*(e2(:,n)-e4(:,n))) + ...
        C_pls_btm(:,n) - C_mns_btm(:,n) - direct(:,n);

    % Mean intensity at the base of each layer (Eq. 49):
    intensity(1:nbr_wvl,n) = (1/mu_one).*...
        ( Y(:,(2*n-1)).*(e1(:,n)+e3(:,n)) + ...
          Y(:,(2*n)).*(e2(:,n)+e4(:,n)) + C_pls_btm(:,n) + C_mns_btm(:,n)) +...
        (direct(:,n)./mu_not);
  
    intensity(1:nbr_wvl,n) = intensity(1:nbr_wvl,n)./(4*pi);
end
  

% Upward flux at upper model boundary (Eq. 31):
F_top_pls = (Y(:,1).*(exp(-lambda(:,1).*tau_star(:,1))+GAMMA(:,1))) + ...
    (Y(:,2).*(exp(-lambda(:,1).*tau_star(:,1))-GAMMA(:,1))) + ...
    C_pls_top(:,1);


for n=1:nbr_lyr
    % Upward flux at the bottom of each layer interface (Eq. 31)
    F_up(1:nbr_wvl,n) = ...
        Y(:,2*n-1).*(exp(0) + GAMMA(:,n).*exp(-lambda(:,n).*tau_star(:,n))) +...
        Y(:,2*n)  .*(exp(0) - GAMMA(:,n).*exp(-lambda(:,n).*tau_star(:,n))) +...
        C_pls_btm(:,n);
  
    % Downward flux at the bottom of each layer interface (Eq. 32,
    % plus direct-beam component):
    F_down(1:nbr_wvl,n) = ...
        Y(:,2*n-1).*(GAMMA(:,n).*exp(0) + exp(-lambda(:,n).*tau_star(:,n))) + ...
        Y(:,2*n)  .*(GAMMA(:,n).*exp(0) - exp(-lambda(:,n).*tau_star(:,n))) + ...
        C_mns_btm(:,n) + direct(:,n);

    % Derived net (upward-downward) flux
    % (should equal F_net, Eq. 48)
    F_net2(1:nbr_wvl,n) = F_up(1:nbr_wvl,n) - F_down(1:nbr_wvl,n);
  
    % planar intensity:
    intensity2(1:nbr_wvl,n) = F_up(1:nbr_wvl,n) + F_down(1:nbr_wvl,n);
end

% surface planar intensity:
intensity2_top(1:nbr_wvl) = F_top_pls + ((mu_not*pi*Fs)+Fd);


% Net spectral flux at lower model boundary = bulk transmission through entire
% media = absorbed radiation by underlying surface [W/m2/band]:
F_btm_net = -F_net(:,nbr_lyr);


% Hemispheric spectral albedo:
albedo = F_top_pls./((mu_not*pi*Fs)+Fd);


% Net spectral flux at upper model boundary [W/m2/band]:
F_top_net(1:nbr_wvl,1) = F_top_pls - ((mu_not*pi*Fs)+Fd);

  
% Absorbed flux in each layer [W/m2/band]:
for n=1:nbr_lyr
    if(n==1)
        F_abs(1:nbr_wvl,1) = F_net(:,1)-F_top_net;
    else
        F_abs(:,n) = F_net(:,n) - F_net(:,n-1);
    end;
end


% Spectrally-integrated absorption in each layer [W/m2/layer]:
for n=1:nbr_lyr
    F_abs_slr(n) = sum(F_abs(:,n));
    F_abs_vis(n) = sum(F_abs(1:vis_max_idx,n));
    F_abs_nir(n) = sum(F_abs(vis_max_idx+1:nir_max_idx,n));
end

% Spectrally-integrated absorption by underlying surface [W/m2]:
F_abs_btm     = sum(F_btm_net);
F_abs_vis_btm = sum(F_btm_net(1:vis_max_idx));
F_abs_nir_btm = sum(F_btm_net(vis_max_idx+1:nir_max_idx));


% Radiative heating rate in each layer:
heating_rate = F_abs_slr./(L_ice.*2117);   %[K/s], 2117 = specific heat ice (J kg-1 K-1)
heating_rate = heating_rate.*3600;         %[K/hr]
				      
				      
% Energy conservation check:
% Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
energy_sum = (mu_not*pi*Fs)+Fd - (sum(F_abs,2) + F_btm_net + F_top_pls);

if (sum(abs(energy_sum)) > 1e-10)
    energy_conservation_error = sum(abs(energy_sum))
    %error(strcat('Energy conservation error of: ',num2str(sum(abs(energy_sum)))));
end;

% spectrally-integrated terms (remove semi-colons to print these):
sum(energy_sum);        % energy conservation total error
sum((mu_not*pi*Fs)+Fd); % total incident insolation (W m-2)
sum(sum(F_abs));        % total energy absorbed by all snow layers
sum(F_btm_net);         % total energy absorbed by underlying substrate


% Spectrally-integrated solar, visible, and near-IR albedos:
alb_bb  = sum(flx_slr.*albedo)./sum(flx_slr);

alb_vis = sum(flx_slr(1:vis_max_idx).*albedo(1:vis_max_idx))/...
          sum(flx_slr(1:vis_max_idx));

alb_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*albedo(vis_max_idx+1:nir_max_idx))/...
          sum(flx_slr(vis_max_idx+1:nir_max_idx));


% Spectrally-integrated visible and near-IR total snowpack+ground
% absorption [W/m2]:
abs_vis = sum(flx_slr(1:vis_max_idx).*(1-albedo(1:vis_max_idx)));
abs_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*(1-albedo(vis_max_idx+1:nir_max_idx)));


% Spectrally-integrated downwelling solar fluxes at top of snowpack [W/m2]:
flx_dwn_top_slr  = sum((mu_not*pi*Fs))+sum(Fd);
flx_dwn_top_vis  = sum((mu_not*pi*Fs(1:vis_max_idx)))+sum(Fd(1:vis_max_idx));
flx_dwn_top_nir  = sum((mu_not*pi*Fs(vis_max_idx+1:nir_max_idx)))+sum(Fd(vis_max_idx+1:nir_max_idx));



%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_out.wvl             = wvl;                % spectral wavelength bands (um)
data_out.flx_dwn_spc     = (mu_not*pi*Fs+Fd)'; % spectral downwelling flux at model top [W/m2/band]
data_out.albedo          = albedo';            % spectral hemispheric albedo
data_out.alb_slr         = alb_bb;             % solar broadband albedo
data_out.alb_vis         = alb_vis;            % visible (0.2-0.7um) albedo
data_out.alb_nir         = alb_nir;            % near-IR (0.7-5.0um) albedo

data_out.abs_snw_slr     = sum(F_abs_slr);     % total solar absorption by entire snow column (not including underlying substrate) [W/m2]
data_out.abs_snw_vis     = sum(F_abs_vis);     % visible solar absorption by entire snow column (not including underlying substrate) [W/m2]
data_out.abs_snw_nir     = sum(F_abs_nir);     % near-IR solar absorption by entire snow column (not including underlying substrate) [W/m2]
data_out.abs_spc         = sum(F_abs,2)';      % spectral absorption by entire snow column [W/m2/band] 

data_out.abs_snw_top_slr = F_abs_slr(1);       % top snow layer solar absorption [W/m2]
data_out.abs_snw_top_vis = F_abs_vis(1);       % top snow layer VIS absorption [W/m2]
data_out.abs_snw_top_nir = F_abs_nir(1);       % top snow layer NIR absorption [W/m2]

data_out.abs_ground_slr  = F_abs_btm;          % total solar absorption by underlying substrate [W/m2]
data_out.abs_ground_vis  = F_abs_vis_btm;      % visible absorption by underlying substrate [W/m2]
data_out.abs_ground_nir  = F_abs_nir_btm;      % near-IR absorption by underlying substrate [W/m2]

data_out.flx_dwn_top_slr = flx_dwn_top_slr;    % downwelling solar broadband flux on upper boundary [W/m2]
data_out.flx_dwn_top_vis = flx_dwn_top_vis;    % downwelling visible broadband flux on upper boundary [W/m2]
data_out.flx_dwn_top_nir = flx_dwn_top_nir;    % downwelling visible broadband flux on upper boundary [W/m2]


if (1==0)
    hold on;
    plot(wvl,albedo,'linewidth',3);
    axis([0.2 2.5 0 1]);
    grid on;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%  REFERENCES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Balkanski, Y., M. Schulz, T. Claquin, and S. Guibert (2007),
% Reevaluation of mineral aerosol radiative forcings suggests a better
% agreement with satellite and AERONET data, Atmos. Chem. Phys., 7
% (1), 81-95, doi:10.5194/acp-7-81-2007.
%
% Cook, J. M., Hodson, A. J., Gardner, A. S., Flanner, M., Tedstone,
% A. J., Williamson, C., Irvine-Fynn, T. D. L., Nilsson, J., Bryant,
% R., and Tranter, M. (2017), Quantifying bioalbedo: a new physically
% based model and discussion of empirical methods for characterising
% biological influence on ice and snow albedo, The Cryosphere, 11,
% 2611-2632, doi: 10.5194/tc-11-2611-2017.
%
%
% Dang, C., Zender, C. S., and Flanner, M. G. (2019), Intercomparison
% and improvement of two-stream shortwave radiative transfer schemes
% in Earth system models for a unified treatment of cryospheric
% surfaces, The Cryosphere, 13, 2325-2343,
% doi:10.5194/tc-13-2325-2019.
%
% Flanner, M. G., A. S. Gardner, S. Eckhardt, A. Stohl, J. Perket
% (2014), Aerosol radiative forcing from the 2010 Eyjafjallajokull
% volcanic eruptions, J. Geophys. Res. Atmos., 119, 9481-9491,
% doi:10.1002/2014JD021977.
%
% Flanner, M. G., C. S. Zender, J. T. Randerson, and P. J. Rasch
% (2007), Present day climate forcing and response from black carbon
% in snow, J. Geophys. Res., 112, D11202, doi:10.1029/2006JD008003.
%
% Flanner, M. G., and C. S. Zender (2005), Snowpack radiative heating:
% Influence on Tibetan Plateau climate, Geophys. Res. Lett., 32,
% L06501, doi:10.1029/2004GL022076.
%
% Fu, Q. (2007), A new parameterization of an asymmetry factor of
% cirrus clouds for climate models, J. Atmos. Sci., 64, 4140-4150,
% https://doi.org/10.1175/2007JAS2289.1.
%
% Hansen, G. B. (2005), Ultraviolet to near-infrared absorption
% spectrum of carbon dioxide ice from 0.174 to 1.8 um,
% J. Geophys. Res., 110, E11003, doi:10.1029/2005JE002531.
%
% He, C., Y. Takano, K.-N. Liou, et al. (2017), Impact of snow grain
% shape and black carbon-snow internal mixing on snow optical
% properties: Parameterizations for climate models, J. Climate,
% 30(24), 10019-10036, doi:10.1175/JCLI-D-17-0300.1.
%
% Lawrence, D. and others (2018), Technical Description of version 5.0
% of the Community Land Model (CLM),
% http://www.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
%
% Picard, G., Libois, Q., and Arnaud, L. (2016), Refinement of the ice
% absorption spectrum in the visible using radiance profile
% measurements in Antarctic snow, The Cryosphere, 10, 2655-2672,
% doi:10.5194/tc-10-2655-2016.
%
% Polashenski, C. M., J. E. Dibb, M. G. Flanner, J. Y. Chen,
% Z. R. Courville, A. M. Lai, J. J. Schauer, M. M. Shafer, and
% M. Bergin (2015), Neither dust nor black carbon causing apparent
% albedo decline in Greenland's dry snow zone: Implications for MODIS
% C5 surface reflectance, Geophys. Res. Lett., 42, 9319-9327,
% doi:10.1002/2015GL065912.
%
% Singh, D., and M. G. Flanner (2016), An improved carbon dioxide snow
% spectral albedo model: Application to Martian conditions,
% J. Geophys. Res. Planets, 121, 2037-2054, doi:10.1002/2016JE005040.
%
% Skiles, S. McKenzie, Painter, Thomas, and Okin, Greg (2017), A
% method to retrieve the spectral complex refractive index and single
% scattering optical properties of dust deposited in mountain snow,
% Journal of Glaciology. Vol. 63, 133-147, doi:10.1017/jog.2016.126.
%
% Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam (1989),
% Rapid calculation of radiative heating rates and photodissociation
% rates in inhomogeneous multiple scattering atmospheres,
% J. Geophys. Res., 94(D13), 16,287-16,301.
%
% Warren, S. G. (1984), Optical constants of ice from the ultraviolet
% to the microwave, Appl. Opt., 23, 1206-1225.
%
% Warren, S. G., and R. E. Brandt (2008), Optical constants of ice
% from the ultraviolet to the microwave: A revised compilation,
% J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744.
%
% Wolff, M. J., M. D. Smith, R. T. Clancy, R. Arvidson, M. Kahre,
% F. Seelos IV, S. Murchie, and H. Savijarvi (2009), Wavelength
% dependence of dust aerosol single scattering albedo as observed by
% the Compact Reconnaissance Imaging Spectrometer, J. Geophys. Res.,
% 114, E00D04, doi:10.1029/2009JE003350.
%
% Wolff, M. J., R. T. Clancy, J. D. Goguen, M. C. Malin, and
% B. A. Cantor (2010), Ultraviolet dust aerosol properties as observed
% by MARCI, Icarus, 208(1), 143-155.



% Copyright (c) 2020 Mark Flanner and SNICAR contributors
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
