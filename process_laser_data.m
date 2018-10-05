%% Processsing laser data
% this is a trial script to try and incorporate the 

% import laser data
[t_daq,Diode,~,~,x,y,~] = importfile(laserdatafn,laserdata.startrow,laserdata.endrow);
