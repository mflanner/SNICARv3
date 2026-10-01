function sza = szaLat(latitude,varargin)
%Average solar zenith angle for a given latitude
% By default, the function evaluates the average solar zenith angle at a
% given latitude for a full calendar year. However, optional arguments make
% it possible to average over a specified range of days or months. Two
% things to note: (1) averaging over the full year is only useful for
% tropical latitudes because the average solar zenith angle for a latitude
% outside of the tropics will always be equal to the latitude itself. For
% example, szaLat(45) will return 45. (2) This function is specifically
% designed to work with the Latitude=N argument in the snicar() function.
% If you specify Latitude=N using the snicar() function, then szaLat(N)
% will be called and overwrite the solar zenith angle specified by the
% user. In other words, calling
%
%  >> snicar(SZA=30, Latitude=10)
%
% will completely replace the SZA argument with the output of szaLat(10).
%
%
%SYNTAX
% sza = szaLat(latitude)
% sza = szaLat(latitude,Day=N)
% sza = szaLat(latitude,Month=N)
%
%
%REQUIRED INPUTS
% latitude :: Latitude in degrees (must be > -90 and < 90).
%
%NAME-VALUE ARGUMENTS
% Day=N   :: Day of the year (must be >0 and <366). Can specify a scalar
%            value or a vector. For 01 January, set Day=1. For 13 March,
%            set Day=73. For 31 December, set Day=365. For 01 January
%            through 25 January, set Day=1:25 or Day=[1 2 3 {...} 25]. 
%
% Month=N :: Month of the year (must be >0 and <13). Can specify a scalar
%            value or a vector. For January, set Month=1. For March, set
%            Month=3. For a variety of months, let's say Dec-Jan-Feb, set
%            Month=[12 1 2]; 
%
%
%EXAMPLES
% sza = szaLat(30); % Will return 30 because 30° latitude is outside of
%                   % the tropics
%
% sza = szaLat(-9); % Will return 16.0424, which is the average solar
%                   % zenith angle for the full year at -9° latitude.
%
% sza = szaLat(-9, Day=31); % Will return 8.7823, the solar zenith angle
%                           % at -9° latitude on 31 January.
%
% sza = szaLat(30, Day=31); % Will return 47.7823, which is the solar
%                           % zenith angle at 30° latitude on 31 Jan. Note
%                           % that the SZA for extra-tropical latitudes
%                           % varies throughout the year. It is only when
%                           % you average the SZAs for each day of the year
%                           % that the szaLat() function would return the
%                           % extra-tropical latitude itself.
%
% sza = szaLat(-9, Day=[1,15,33]; % Will return 11.5089, which is the mean
%                                 % solar zenith angle at -9° latitude for
%                                 % the days 01 Jan, 15 Jan, and 02 Feb.
%
% sza = szaLat(-9, Month=7); % Will return 30.1015, which is the average
%                            % solar zenith angle at -9° latitude during
%                            % seventh month (July).
%
% sza - szaLat(-9, Month=[12 1 2]); % Will return 10.2888, which is the
%                                   % average solar zenith angle at -9°
%                                   % latitude during the months of Dec,
%                                   % Jan, and Feb.
%
%
%NOTES
% Solar zenith angle (θ_z) is defined by Duffie & Beckman (1980, p.13) as 
% "the angle between the vertical and the line to the sun, that is, the
% angle of incidence of beam radiation on a horizontal surface. Solar
% declination (δ) is the "angular position of the sun at solar noon" 
% (p.12). The declination can be approximated by Cooper's equation 
% (Cooper, 1969):
%
%    δ = 23.45sin(360*((284+n)/365))
%
% where `n` is the day of the year (01 Jan = 1; 31 Dec = 365). A slightly
% better estimate is defined using Spencer's equation (Spencer, 1971):
%
%    δ = (180/π)(0.006918-0.399912*cos(B)+0.070257*sin(B)-0.006758*cos(2B)+
%         0.000907*sin(2B)-0.002697*cos(3B)+0.00148*sin(3B))
%
% where B is equal to (n-1)*(360/365). Solar zenith angle may be estimated
% (in the special case of solar noon) by the formula:
%
%    θ_z = |φ − δ|
%
% where `φ` is latitude in degrees and `δ` is the declination (see Duffie
% & Beckman, 1980, p. 17, Eq. 1.6.9).
%
%REFERENCES
% Cooper, P. I. (1969). The Absorption of Solar Radiation in Solar Stills.
% Solar Energy, 12(3).
%
% Duffie, John A., and William A. Beckman (1980). Solar engineering of
% thermal processes. New York: Wiley.
%
% Spencer, J. W. (1971). Fourier Series Representation of the Position of
% the Sun. Search, 2(5), 172.
%
%
%See also
% snicar

latitude_validation_function = @(x) isnumeric(x) & isscalar(x) ...
                                    & (x>-90) & (x<90);

day_validation_function =     @(x) isnumeric(x) & isvector(x) ...
                                   & all(x>=1 & x<=365) ...
                                   & (numel(unique(x))==numel(x));

month_validation_function =   @(x) isnumeric(x) & isvector(x) ...
                                   & all(x>=1 & x<=365) ...
                                   & (numel(unique(x))==numel(x));

% Input parsing
inP = inputParser();
addRequired(inP,'latitude',latitude_validation_function);
addParameter(inP,'day',[],day_validation_function);
addParameter(inP,'month',[],month_validation_function);
parse(inP,latitude,varargin{:});
d = inP.Results.day;
m = inP.Results.month;

% Main function
coopers_fun = @(n) 23.45*sind((360/365)*(284+n)); % Solar declination

if (isempty(d)) & (isempty(m))
%Calculate average solar zenith angle for the full year
  sza365 = zeros([365 1]);
  for day_of_year = 1:length(sza365)
    sza365(day_of_year) = abs(latitude - coopers_fun(day_of_year));
  end
  sza = mean(sza365);
  return
end

if ~isempty(d)
%Evaluate solar zenith angle for the given day of the year
  if ~all(d==floor(d))
    error(['Invalid input for the ''Day'' parameter. ' ...
           'Days may only be expressed as integers!'])
  end
  if ~isempty(m)
    warning(['Cannot specify a day of year AND a month! '...
             'Evaluating solar zenith angle using the day only.'])
  end
  num_days = numel(d);
  szaD = zeros([num_days 1]);
  for day_of_year = 1:num_days
    szaD(day_of_year) = abs(latitude - coopers_fun(d(day_of_year)));
  end
  sza = mean(szaD);
  return
end


% Reaching this point means 'month' was specified but not 'day'. Creating
% dictionaries defining the days of the year in each month.
month_index = dictionary(1,'Jan',2,'Feb',3,'Mar',4,'Apr',5,'May',...
                         6,'Jun',7,'Jul',8,'Aug',9,'Sep',10,'Oct',...
                         11,'Nov',12,'Dec');
month_start_day = dictionary('Jan',1,'Feb',32,'Mar',60,'Apr',91,...
                             'May',121,'Jun',152,'Jul',182,'Aug',213,...
                             'Sep',244,'Oct',274,'Nov',305,'Dec',335);
month_end_day = dictionary('Jan',31,'Feb',59,'Mar',90,'Apr',120,...
                           'May',151,'Jun',181,'Jul',212,'Aug',243,...
                           'Sep',273,'Oct',304,'Nov',334,'Dec',365);

if ~isempty(m)
  if ~all(m==floor(m))
    error(['Invalid input for the ''Month'' parameter. ' ...
      'Months may only be expressed as integers!'])
  end
%Evaluate solar zenith angle for the given month of the year
  num_months = numel(m);
  days_of_year = [];
  for month_of_year = 1:num_months
    selected_month = month_index(m(month_of_year));
    days_of_month = month_start_day(selected_month) : month_end_day(selected_month);
    days_of_year = [days_of_year days_of_month];
  end
  num_days = numel(days_of_year);
  szaMon = zeros([num_days 1]);
  for day_of_year = 1:num_days
    szaMon(day_of_year) = abs(latitude - coopers_fun(days_of_year(day_of_year)));
  end
  sza = mean(szaMon);
  return

else
  error('Something went wrong with parsing day/month arguments.')
end
