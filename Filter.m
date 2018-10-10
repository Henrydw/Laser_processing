

%% imports:

% ren_cal(.csv) data:
path_ren_cal = "calibration_data/ren_cal.csv";
optical_data = csvread(path_ren_cal,2,0);

% load spectral repsonse of camera:
path_spectral = "calibration_data/SA5_Spectral_Response_Curve.csv";
spectral = csvread(path_spectral,2,0);