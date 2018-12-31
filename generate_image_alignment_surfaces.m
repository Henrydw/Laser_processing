%% Generate image aligment surfaces
% a program to generate the image offset surafaces required to properly
% centre images on the laser.

% REASON: the f-theta lens induces an offset since the imageing wavelengths
% (700nm & 950nm) do not match the wavelength the lens is designed for 
% (1070nm) . Therefore the point of laser contact (meltpool) appears to 
% drift around the image when it should stay centre, this occurs by 
% different amounts for each image due to their seperate wavelenghts:

% ISSUE: if the melt pools cannot be properly aligned then the ratio 
% (and consequnetly the temperature) cannot be ascertained

% Inputs:

% folder_path - path of folder where the data folder is contained (it is
% assumed that the exact folder with the data will be at the end of this 
% path - eg: 'A:\Imperial College London\Hooper, Paul A - spots_v3'

% folder_name - name of the exact folder with the data in it and also the
% name of the .mat file which should contain the DAQ data - eg: 
% '100W_6400us_100000fps'

% NB this folder should contain the following:
% 100W_6400us_100000fps.mat - DAQ data
% C001H001S0001.cih - camera 1 meta data
% C002H001S0001.cih - camera 2 meta data
% C001H001S0001.mraw - camera 1 recording
% C002H001S0001.mraw - camera 2 recording


function [fits]=generate_image_alignment_surfaces(folder_path,folder_name)


% load cih data (camera information header) - creates class: imagedata
imagedata = readcih(strcat(folder_path,'\',folder_name));

% FROM asset.m NOT findspotfit.m:
% Specify the precision for reading images
imagedata.Precision = strcat('*ubit',num2str(imagedata.ColorBit));
% Specify the bits per pixel according to image data
bits_per_pixel = imagedata.ColorBit;
% Always a constant
bits_per_byte = 8;
% Specify the total bytes packed in a frame
imagedata.NBytesPerFrame = imagedata.Pixels * bits_per_pixel / bits_per_byte;

%% --------------------------- image processing

% ----------------------- load mraw settings

%number of frames being processed in total:
end_frame = imagedata.TotalFrames - 1;
end_frame = 645;

% number to images to be process per block (all at once crashes computer?)
blockSize = 100;
% tnmwfi!!!!

% ----------------------- prep

% Calculate time vector for images up to speficied end frame
start_frame = 1;
vdt = 1/imagedata.FrameRate;
t_frame = (0:(end_frame-1)) * vdt; %video frame time
t_frame = round(t_frame,9); %round to ns to avoid rounding errors

% Filter for peak (spot finding) function
filt=(fspecial('gaussian', 7,1));

% Initialise arrays
% %matrix to store x,y pos of spot in each frame
IC1_p = zeros(2,end_frame); 
IC2_p = zeros(2,end_frame);


% NOT WORKING (WHERE IS THE ESTIMATED NOISE LEVEL STORED?CALCULATED)
% Remove background noise (for images when laser is off (no spot))
%IC1(IC1 < estimated_noise_cam1) = 0;
%IC2(IC2 < estimated_noise_cam2) = 0;
% IN TANDEM WITH: *

%% ---------------------------------- diode processing

% Read laser data file
[t_daq,Diode,~,~,x,y,~] = importfile(strcat(folder_path,'\',folder_name,'\',folder_name,'.mat'));

% Create a kaiser filter to smooth DAQ data

% Hard coded settings
F_PASS = 0; % Pass frequency
F_STOP = 2500; % Stop frequency
RIPPLE = 0.0001; % Max bandpass ripple
% Create filter
dt = t_daq(3) - t_daq(2); % DAQ sample interval
Fs = 1/dt; % Sampling freq
[n, w, beta, ftype] = kaiserord([F_PASS,F_STOP], [1,0], [RIPPLE,RIPPLE], Fs);
b = fir1(n,w,ftype,kaiser(n+1,beta),'noscale');

% Filter the data and then get xy mirror pos values at video frame times
y_pos = interp1(t_daq,filtfilt(b,1,y(1:end)),t_frame);
x_pos = interp1(t_daq,filtfilt(b,1,x(1:end)),t_frame);
% Load in laser pulse data and get value at video frame times
lpulse = interp1(t_daq,Diode(1:end),t_frame);

% Find indexes when laser is on (to reduce noise from times when it is off)
lon = lpulse > 0.1; %0.1 seems a good value for 200W pulses

%----------------------- finding peaks

odds = [1:2:1000];

%estimated noise
estimated_noise_cam = 100;


% frame number
frame_number=1;
for j=1:blockSize:end_frame
    % load block
    %first and last frame to be loaded
    if frame_number+1+blockSize > end_frame
        frange = [frame_number+1,end_frame];
        blockSize=range(frange);
    else
        frange = [frame_number+1,frame_number+1+blockSize];
    end
    
    %read in frame from mraw files for each camera
    IC1=readmraw(imagedata.folderName,imagedata,frange,1);
    IC2=readmraw(imagedata.folderName,imagedata,frange,2);
    
    % Remove background noise (for images when laser is off (no spot))
    IC1(IC1 < estimated_noise_cam) = 0;
    IC2(IC2 < estimated_noise_cam) = 0;
    
    for i=1:blockSize

        frame_number = frame_number+1;

        % Calculate threshold params for peak finder function
        thresh=(max([min(max(IC1(:,:,i),[],1))  min(max(IC1(:,:,i),[],2))]));
        %find some peaks for cam 1
        p = FastPeakFind(IC1(:,:,i),thresh,filt,3,2);
        % Check to see if peaks were found
        % Obtain peak in terms of pixel value (not real size)
        % Save only coords of first found peak
        
        intensities = [];
                

    %     % *
    %     % if zero, save as NaN - allows frames which do not have laser in them
    %     % to be detected - surely this is already done by the Diode data???
        %if zero, save as nan
        if (size(p) > 0)
            %search for max peak and save
            for q = 1:2:size(p)
                intensities = [intensities,IC1(round(p(q)),round(p(q+1)),i)];
            end
            true_peak = find(intensities==max(intensities));

            IC1_p(:,frame_number) = p(odds(true_peak):odds(true_peak)+1);
        else
            IC1_p(:,frame_number) = [NaN NaN];
        end

        % Repeat above but for other camera
        thresh = (max([min(max(IC2(:,:,i),[],1))  min(max(IC2(:,:,i),[],2))]));
        p = FastPeakFind(IC2(:,:,i),thresh,filt,3,2);

        intensities = [];
        
        %     % *
        if (size(p) > 0)
            %search for max peak and save
            for q = 1:2:size(p)
                intensities = [intensities,IC2(round(p(q)),round(p(q+1)),i)];
            end
            true_peak = find(intensities==max(intensities));

            IC2_p(:,frame_number) = p(odds(true_peak):odds(true_peak)+1);
        else
            IC2_p(:,frame_number) = [NaN NaN];
        end
    end
end

% -------------------------- find fits

%only use row where a spot was found and laser on (fixes warnings about NaN
%in fit data)
idxfit1=isfinite(IC1_p(1,:))&lon;
idxfit2=isfinite(IC2_p(1,:))&lon;

datafilename = "IC_fits";
% Call fiting function to find fit values
% do this only when laser is firing
fits.IC1_x=findoffsetfit(x_pos(idxfit1),y_pos(idxfit1),IC1_p(1,idxfit1),datafilename,'IC1_x_fit',end_frame);
fits.IC1_y=findoffsetfit(x_pos(idxfit1),y_pos(idxfit1),IC1_p(2,idxfit1),datafilename,'IC1_y_fit',end_frame);
fits.IC2_x=findoffsetfit(x_pos(idxfit2),y_pos(idxfit2),IC2_p(1,idxfit2),datafilename,'IC2_x_fit',end_frame);
fits.IC2_y=findoffsetfit(x_pos(idxfit2),y_pos(idxfit2),IC2_p(2,idxfit2),datafilename,'IC2_y_fit',end_frame);

end