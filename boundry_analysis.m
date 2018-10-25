% Boundry analysis
% a script to iterate though the single exposure data and caluclate:
% G - temp gradient across liquid-solid interface
% R - interface velocity


file = 
% filenames
file_name = strcat('\','100W_6400us_100000fps');

% folder location
folder_path = 'A:\Imperial College London\Hooper, Paul A - spots_v3';


[t_frame, spotradius, meanmidspottemp]=process_laser_data(folder_path,file_name);


plot(t_frame,spotradius)
title