% Filter.m
% a Class for the intensity-ratio -> temperature look up table
% intensity ratio is between corresponding pixels at two seperate
% wave lengths.

%% constants:

REF_TEMPERATURE = (293:0.5:5000)';
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)


%% imports:

% ren_cal(.csv):
path_ren_cal = "calibration_data/ren_cal.csv";
optical_data = csvread(path_ren_cal,2,0);

% spectral_response repsonse of spectral_response:
path_spectral_response = "calibration_data/SA5_Spectral_Response_Curve.csv";
spectral_response = csvread(path_spectral_response,2,0);

% band pass filter efficiency curves
BPFilter700 = csvread("calibration_data/eff700.csv",2,0);
BPFilter950 = csvread("calibration_data/eff950.csv",2,0);

%% calculations - done for each filter frequency: - currently 700nm

bp_filter = BPFilter700;
%calc filter parameters from filter data
fon=bp_filter(find(bp_filter(:,2)>(0.5*max(bp_filter(:,2))),1,'first'),1);
foff=bp_filter(find(bp_filter(:,2)>(0.5*max(bp_filter(:,2))),1,'last'),1);
fwidth=round(foff-fon);
fcentre=round((foff+fon)/2);

%range of wave lengths to do the calc (just do around the filter wave length).
wl=[(fcentre-fwidth):(fcentre+fwidth)]';

% Create grid of wl and T so we can calculte I as both vary
[WL, T] = meshgrid(wl,REF_TEMPERATURE);

% Calculate intesity emmitted vs wavelength and tempature
I = bb_spectrum(WL,T);

% Calculate amount of light collected via lens
% Solid_angle = Area_of_lens/(dist_to_object^2)
SA = pi()*10^2/420^2; % Focal length is 420mm, limiting appature is r=10mm
Asource = (ORIGINAL_PIXEL_SIZE)^2; % Area of source (20um pixel size with 1:1 mag with high speed spectral_response)
TP = Asource*SA; % Throughput in m^2.strad

I = I * TP; % Watts/m - energy emitted/second per m wave length

% Load optics transmission effeciency vs wavelength data
optical_data(:,2) = optical_data(:,2)/max(optical_data(:,2)); % normalize data
optics_eff = interp1(optical_data(:,1),optical_data(:,2),WL); % get efficiency at desired wavelengths
 
I = I .* optics_eff; % correct for optics eff

% apply colour filter
filter_eff=interp1(bp_filter(:,1),bp_filter(:,2),WL);
I=I.*filter_eff;

% camera responsivity data vs wave length
camera_resp = interp1(spectral_response(:,1),spectral_response(:,2),WL); %Amps/Watt

I = I .* camera_resp; % Camera signal in amps/m wavelength

Itotal = trapz(wl * 10^-9,I,2); % Integrate over wl to get total sensor ampage

%%

% 1 amps(A) is 1 coulomb/s  which is 6.241*10^18 e-/s as a function of elementary charge e-
% Signal output is amps*6.241e18*exposuretime/(well_capcity)*bitdepth
coloumb = 6.241e18;
ex_time = 10e-5; % Exposure time (Can vary to change the temperature range)
wc = 45000; % Well capacity of sensor, units e-
bd = 2^(imagedata.ColorBit) - 1; % 12 bit sensor 
            
% Calculate intensity value 
% Note that if temperature range is to be determined, need to
% multiply by an emissivity value ***
Itotal = (Itotal * coloumb * ex_time / wc) * bd;
