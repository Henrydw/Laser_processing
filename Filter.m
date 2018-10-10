

%% imports:

% ren_cal(.csv):
path_ren_cal = "calibration_data/ren_cal.csv";
optical_data = csvread(path_ren_cal,2,0);

% spectral repsonse of camera:
path_spectral = "calibration_data/SA5_Spectral_Response_Curve.csv";
spectral = csvread(path_spectral,2,0);

% band pass filter efficiency curves
BPFilter700 = csvread("calibration_data/eff700",2,0);
BPFilter950 = csvread("calibration_data/eff950",2,0);
