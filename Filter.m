% Filter.m
% a Class for the ratio to temperature look up

%% constants:

REF_TEMPERATURE = (293:0.5:5000)';



%% imports:

% ren_cal(.csv):
path_ren_cal = "calibration_data/ren_cal.csv";
optical_data = csvread(path_ren_cal,2,0);

% spectral repsonse of camera:
path_spectral = "calibration_data/SA5_Spectral_Response_Curve.csv";
spectral = csvread(path_spectral,2,0);

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