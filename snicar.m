function results = snicar(varargin)
%SNICARv3 (Flanner et al., 2021)
%   A user-friendly version of the snicarAD_v3 function
%
%Syntax
% results = snicar(Name=Value)
%
% where
%  Name=Value :: Input arguments given as name-value pairs
%  results    :: Struct object containing the model results. See
%                results.README for details on each field in the struct.
%
%===============================================================================
%Example
% >> results = snicar(SZA=30,...
%                     AtmosphereType=6,...
%                     SnowpackDensity=400,...
%                     SnowGrainRadius=500,...
%                     AlbedoGround=0.55,...
%                     BlackCarbon=5.0,...
%                     Dust=[2.1 6.9 3.4 0.9 0.1],...
%                     DustType=3);
%                     
%
%===============================================================================
%Name-Value Pairs (default values shown first)
% NAME                    VALUE
% 'IncidentRadiation'     1=direct beam (clear sky)
%                         0=diffuse (cloudy sky)
%
% 'SolarZenithAngle'      60 (degrees, can be any value 0-89 degrees; using
%                             'SZA' for the name argument is also valid)
%
% 'AtmosphereType'        1=mid-latitude winter
%                         2=mid-latitude summer
%                         3=sub-arctic winter
%                         4=sub-arctic summer
%                         5=Summit, Greenland (sub-Arctic summer, 796 hPa)
%                         6=High Mountain (summer, 556 hPa)
%                         7=Top of Atmosphere
%                         
% 'SnowpackThickness'     100 (meters)
%
% 'SnowpackDensity'       200 (kg/m^3)
%
% 'IceRefractiveIndex'    3=Picard et al. (2016) / Warren & Brandt (2008)
%                         1=Warren (1984) / Perovich & Govoni (1991)
%                         2=Warren & Brandt (2008)
%                         4=CO2 ice (Hansen, 2005; Singh & Flanner, 2016)
%
% 'SnowGrainRadius'       100 (microns, can be any value 30-1500 microns)
%
% 'SnowGrainShape'        3=Hexagonal plates (aspect ratio 2.5)
%                         1=Spheres
%                         2=Spheroids (aspect ratio 0.5)
%                         4=Koch snowflakes (aspect ratio 2.5)
%
% 'AlbedoGround'          0.25 (Albedo of underlying ground, 0-1)
%
% 'BlackCarbon'           0 (ppb)
%
% 'BlackCaronCoated'      0 (ppb, Sulfate-coated black carbon)
%
% 'BrownCarbon'           0 (ppb)
%
% 'BrownCarbonCoated'     0 (ppb, Sulfate-coated brown carbon)
%
% 'Dust'                  [0 0 0 0 0] (ppm, given as vector for sizes:
%                                      Size 1 = 0.1-1.0um
%                                      Size 2 = 1.0-2.5um
%                                      Size 3 = 2.5-5.0um
%                                      Size 4 = 5.0-10.0um
%                                      Size 5 = 10.0-50.0um)
%
% 'DustType'              1=Sahara (Balkanski et al., 2007)
%                         2=San Juan Mtns, CO (Skiles, et al., 2017)
%                         3=Greenland (Polashenski et al., 2015)
%                         4=Mars (Wolff et al.,2009,2010; Singh & Flanner, 2016)
%
% 'VolcanicAsh'           [0 0 0 0 0] (ppm, given as vector for sizes:
%                                      Size 1 = 0.1-1.0um
%                                      Size 2 = 1.0-2.5um
%                                      Size 3 = 2.5-5.0um
%                                      Size 4 = 5.0-10.0um
%                                      Size 5 = 10.0-50.0um)
%
% 'SnowAlgae'             0 (cells/mL)
%
% 'AlgaeRadius'           10 (microns, valid values: 1,2,5,10,15,20,25,30,40,50)
%
% 'DryCellMass'           [0.015 0.005 0.05 0] (vector, valid values are:
%                                Elem1: 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03
%                                Elem2: 0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03
%                                Elem3: 0 to 0.15 discretized by 0.01
%                                Elem4: 0 to 0.15 discretized by 0.01
%                                where,
%                                Elem1 = Chlorophyll-a
%                                Elem2 = Chlorophyll-b
%                                Elem3 = Photoprotective carotenoids
%                                Elem4 = Photosynthetic carotenoids
%
%See also
% snicarAD_v3

%Note on copyrights:
%  The original work for SNICAR-ADv3 was published by Flanner et al.
%  (2021a). The source code for the original model is on Github (Flanner et
%  al., 2021b) and Zenodo (Flanner et al., 2021c). The code in this
%  function was written by Austin M. Weber, but the original authors retain
%  the copyrights for the open source software. A copy of the software
%  license is provided at the end of the script. Please continue to cite
%  the appropriate references when using SNICAR-ADv3 in your research.
%
%References
% Flanner, M. G., Arnheim, J. B., Cook, J. M., Dang, C., He, C., Huang, X.,
%  Singh, D., Skiles, S. M., Whicker, C. A., & Zender, C. S. (2021a).
%  SNICAR-ADv3: A community tool for modeling spectral snow albedo.
%  Geoscientific Model Development, 14(12), 7673–7704.
%  https://doi.org/10.5194/gmd-14-7673-2021
%
% Flanner et al. (2021b). SNICAR-ADv3 [Software]. Github.
%  https://github.com/mflanner/SNICARv3. Last accessed 28 Oct. 2025.
%
% Flanner et al. (2021c). SNICAR-ADv3 [Software]. Zenodo.
%  https://doi.org/105281/zenodo.5176213
%


% Input parsing
inP = inputParser();
addParameter(inP,'IncidentRadiation',1,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'SolarZenithAngle',60,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'SZA',60,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'AtmosphereType',1,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'SnowpackThickness',100,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'SnowpackDensity',200,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'IceRefractiveIndex',3,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'SnowGrainRadius',100,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'SnowGrainShape',3,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'AlbedoGround',0.25,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'BlackCarbon',0,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'BlackCarbonCoated',0,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'BrownCarbon',0,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'BrownCarbonCoated',0,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'Dust',[0 0 0 0 0],@(x) isnumeric(x));
addParameter(inP,'DustType',1,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'VolcanicAsh',[0 0 0 0 0],@(x) isnumeric(x));
addParameter(inP,'SnowAlgae',0,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'AlgaeRadius',10,@(x) isnumeric(x) & isscalar(x));
addParameter(inP,'DryCellMass',[0 0 0 0],@(x) isnumeric(x));

parse(inP,varargin{:});
IR = inP.Results.IncidentRadiation;
SZA = inP.Results.SolarZenithAngle;
SZA2 = inP.Results.SZA;
ATM = inP.Results.AtmosphereType;
SPT = inP.Results.SnowpackThickness;
SPD = inP.Results.SnowpackDensity;
IRI = inP.Results.IceRefractiveIndex;
SGR = inP.Results.SnowGrainRadius;
SGS = inP.Results.SnowGrainShape;
AG  = inP.Results.AlbedoGround;
BC  = inP.Results.BlackCarbon;
BCC = inP.Results.BlackCarbonCoated;
BrC = inP.Results.BrownCarbon;
BrCC= inP.Results.BrownCarbonCoated;
Dust= inP.Results.Dust;
DustType = inP.Results.DustType;
VA  = inP.Results.VolcanicAsh;
SnowAlgae = inP.Results.SnowAlgae;
AlgaeRadius=inP.Results.AlgaeRadius;
DSM = inP.Results.DryCellMass;

% Error messages
if IR < 0 || IR > 1
  error('''IncidentRadiation'' can only be set to 1 (direct beam) or 0 (diffuse beam).')
end
if ~isequal(SZA,SZA2)
  if isequal(SZA,60)
    SZA=SZA2;
  end
end
if SZA < 0 || SZA > 89
  error('''SolarZenithAngle'' must be between 0 and 89 degrees.')
end
if ~any(ismember(1:7,ATM))
  error('''AtmosphereType'' must be an integer between 1 and 7.')
end
if SPT < 0
  error('''SnowpackThickness'' must be a positve number.')
end
if SPD <= 0 || SPD >= 1000
  error('''SnowpackDensity'' must be between 0 and 1000 kg/m3, not including 0.')
end
if ~any(ismember(1:4,IRI))
  error('''IceRefractiveIndex'' must be an integer between 1 and 4.')
end
if SGR <= 0
  error('''SnowGrainRadius'' must be a positive number.')
end
if ~any(ismember(1:4,SGS))
  error('''SnowGrainShape'' must be an integer between 1 and 4.')
end
if AG < 0 || AG > 1
  error('''AlbedoGround'' must be between 0 and 1.')
end
if BC < 0 || BCC < 0 || BrC < 0 || BrCC < 0
  error('Black/brown carbon concentrations cannot be negative.')
end
if length(Dust)>5 || length(Dust)<5
  error('''Dust'' must be a 1x5 or 5x1 numeric vector.')
end
if numel(Dust) > 5
  error('''Dust'' must be a 1x5 or 5x1 numeric vector.')
end
if any(Dust<0)
  error('''Dust'' cannot contain negative values.')
end
if ~any(ismember(1:4,DustType))
  error('''DustType'' must be an integer between 1 and 4.')
end
if length(VA)>5 || length(VA)<5
  error('''VolcanicAsh'' must be a 1x5 or 5x1 numeric vector.')
end
if numel(VA) > 5
  error('''VolcanicAsh'' must be a 1x5 or 5x1 numeric vector.')
end
if any(VA<0)
  error('''VolcanicAsh'' cannot contain negative values.')
end
if SnowAlgae < 0
  error('''SnowAlgae'' concentration cannot be negative.')
end
if ~any(ismember([1,2,5,10,15,20,25,30,40,50],AlgaeRadius))
  error('''Algae'' radius can only be set to 1, 2, 5, 10, 15, 20, 25, 30, 40, or 50 microns.')
end
if length(DSM)>4 || length(DSM)<4
  error('''DryCellMass'' must be a 1x4 or 4x1 numeric vector.')
end
if numel(DSM) > 4
  error('''DryCellMass'' must be a 1x4 or 4x1 numeric vector.')
end
if any(DSM<0)
  error('''DryCellMass'' cannot contain negative values.')
end
if ~any(ismember(DSM(1),[0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]))
  error('The first element in ''DryCellMass'' can only be set to 0, 0.005, 0.01, 0.015, 0.02, 0.025, or 0.03.')
end
if ~any(ismember(DSM(2),[0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]))
  error('The second element in ''DryCellMass'' can only be set to 0, 0.005, 0.01, 0.015, 0.02, 0.025, or 0.03.')
end
if ~any(ismember(DSM(3),[0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15]))
  disp(DSM(3))
  error('The third element in ''DryCellMass'' can only be set to a value between 0 and 0.15, discretized by 0.01.')
end
if ~any(ismember(DSM(4),[0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15]))
  error('The fourth element in ''DryCellMass'' can only be set to a value between 0 and 0.15, discretized by 0.01.')
end

% Define model inputs
input_args.direct_beam = IR;
input_args.coszen = cos(deg2rad(SZA));
input_args.atm = ATM;
input_args.dz = SPT;
nbr_lyr = 1:length(input_args.dz);
input_args.rho_snw(1:nbr_lyr) = SPD;
input_args.ice_ri = IRI;
input_args.rds_snw(1:nbr_lyr) = SGR;
input_args.sno_shp(1:nbr_lyr) = SGS;
input_args.R_sfc_all_wvl = AG;
input_args.mss_cnc_sot1(1:nbr_lyr) = BC;
input_args.mss_cnc_sot2(1:nbr_lyr) = BCC;
input_args.mss_cnc_brc1(1:nbr_lyr) = BrC;
input_args.mss_cnc_brc2(1:nbr_lyr) = BrCC;
input_args.mss_cnc_dst1(1:nbr_lyr) = Dust(1);
input_args.mss_cnc_dst2(1:nbr_lyr) = Dust(2);
input_args.mss_cnc_dst3(1:nbr_lyr) = Dust(3);
input_args.mss_cnc_dst4(1:nbr_lyr) = Dust(4);
input_args.mss_cnc_dst5(1:nbr_lyr) = Dust(5);
input_args.dust_type = DustType;
input_args.ash_type = 1; % Eyjafjallajokull
input_args.mss_cnc_ash1(1:nbr_lyr) = VA(1);
input_args.mss_cnc_ash2(1:nbr_lyr) = VA(2);
input_args.mss_cnc_ash3(1:nbr_lyr) = VA(3);
input_args.mss_cnc_ash4(1:nbr_lyr) = VA(4);
input_args.mss_cnc_ash5(1:nbr_lyr) = VA(5);
input_args.cell_nbr_conc(1:nbr_lyr) = SnowAlgae;
input_args.alg_rds(1:nbr_lyr) = AlgaeRadius;
input_args.dcmf_pig_chla(1:nbr_lyr) = DSM(1);
input_args.dcmf_pig_chlb(1:nbr_lyr) = DSM(2);
input_args.dcmf_pig_cara(1:nbr_lyr) = DSM(3);
input_args.dcmf_pig_carb(1:nbr_lyr) = DSM(4);
input_args.flx_dwn_bb = 1;

% Set shape properties
if input_args.sno_shp == 1
  input_args.sno_ar(1:nbr_lyr)   = 0;
elseif input_args.sno_shp == 2
  input_args.sno_ar(1:nbr_lyr)   = 0.5;
else
  input_args.sno_ar(1:nbr_lyr)   = 2.5;
end
if input_args.sno_shp == 1
  input_args.sno_fs(1:nbr_lyr) = 0;
elseif input_args.sno_shp == 2
  input_args.sno_fs(1:nbr_lyr) = 0.929;
elseif input_args.sno_shp == 3
  input_args.sno_fs(1:nbr_lyr) = 0.788;
elseif input_args.sno_shp == 4
  input_args.sno_fs(1:nbr_lyr)   = 0.712;
end

% Run model
results = snicarAD_v3(input_args);

% End function
end

%% License
% MIT License
% 
% Copyright (c) 2020 Mark Flanner and SNICAR contributors
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
