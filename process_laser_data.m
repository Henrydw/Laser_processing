%% Processsing laser data
% this is a list of function primalrily concerned with processing image
% data and returing the paramters from each image analysed.

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


function [t_frame, spotradius, meanmidspottemp]=process_laser_data(folder_path,folder_name)
%% Hard Coded stuff

DEFAULT_FILTER_THRESHOLD = 100;
REF_TEMPERATURE = (293:0.5:5000)';

% from IMAGE DATA class
% Temperature image relative to camera image scale
RESIZED_SCALE = 4;
% In micrometer
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
RESIZED_PIXEL_SIZE = ORIGINAL_PIXEL_SIZE / RESIZED_SCALE;
%from SPOT class
SIZE_SCALE = (RESIZED_PIXEL_SIZE * 1e3).^2;

%% Load calibration

% denote folder with data (should contain .mat *1, .cih *2, .mraw *2)
%data_folder = 'example_data';
%data_folder = "A:\Imperial College London\Hooper, Paul A - spots_v3C";
data_folder = strcat(folder_path,'\',folder_name);
%mat_folder_name = '\100W_6400us_100000fps.mat';
mat_folder_name = strcat('\',folder_name,'.mat');
% load alignment surfaces
load('fitfn_file');

load('intensity_ratio.mat')

% load cih data (camera information header) - creates class: imagedata
imagedata = readcih(data_folder);

% FROM asset.m NOT findspotfit.m:
% Specify the precision for reading images
imagedata.Precision = strcat('*ubit',num2str(imagedata.ColorBit));
% Specify the bits per pixel according to image data
bits_per_pixel = imagedata.ColorBit;
% Always a constant
bits_per_byte = 8;
% Specify the total bytes packed in a frame
imagedata.NBytesPerFrame = imagedata.Pixels * bits_per_pixel / bits_per_byte;

%number of frames being processed in total:
end_frame = imagedata.TotalFrames - 1;

% number to images to be process per block (all at once crashes computer?)
blockSize = 1000;
% tnmwfi!!!!

% Calculate time vector for images up to speficied end frame
start_frame = 1;
vdt = 1/imagedata.FrameRate;
t_frame = (0:(end_frame-1)) * vdt; %video frame time
t_frame = round(t_frame,9); %round to ns to avoid rounding errors

% create meshgrid
[X,Y]=meshgrid(1:imagedata.Width,1:imagedata.Height);


%% load Laser data for Processing
% Read laser data file:
[t_daq,Diode,~,~,x,y,~] = importfile(strcat(data_folder,mat_folder_name));


% Create a kaiser filter to smooth DAQ data:
% Hard coded settings
F_PASS = 0; % Pass frequency
F_STOP = 2500; % Stop frequency
RIPPLE = 0.0001; % Max bandpass ripple
% Create filter
dt = t_daq(3) - t_daq(2); % DAQ sample interval
Fs = 1/dt; % Sampling freq
[n, w, beta, ftype] = kaiserord([F_PASS,F_STOP], [1,0], [RIPPLE,RIPPLE], Fs);
b = fir1(n,w,ftype,kaiser(n+1,beta),'noscale');

% interpolate
% Filter the data and then get xy mirror pos values at video frame times
y_pos = interp1(t_daq,filtfilt(b,1,y(1:end)),t_frame,'linear','extrap');
x_pos = interp1(t_daq,filtfilt(b,1,x(1:end)),t_frame,'linear','extrap');


%% load Image data for Processing

% frame number
frame_number=1;
%stored images 3D array
processed_images = zeros(128,128,imagedata.TotalFrames);
%1d vectors of spot stats in each from
midspottemp = zeros(1,end_frame);
meanmidspottemp = zeros(1,end_frame);
peakspottemp = zeros(1,end_frame);
spotradius = zeros(1,end_frame);

% load blocks loop
for j=1:blockSize:end_frame
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
    % process frames loop
    for i=1:blockSize
        frame_number = frame_number+1;
        
        %CAM1 find offset
        %find spot offsets (in pixels?) for current xy scan pos
        Xe=feval(IC1_x_fit,x_pos(frame_number),y_pos(frame_number)); % Evaluate obtained polynomial function to find offset
        Ye=feval(IC1_y_fit,x_pos(frame_number),y_pos(frame_number)); % Evaluate obtained polynomial function to find offset
        %apply offset, X1 and y1 contain grids of pixel numbers with 0,0 at the
        %laser spot centre
        X1=X-Xe+imagedata.Width/2;
        Y1=Y-Ye+imagedata.Height/2;
        
        %CAM2 find offset
        % Find offset in cam 2 xy positions
        Xe=feval(IC2_x_fit,x_pos(frame_number),y_pos(frame_number)); % Evaluate obtained polynomial function to find offset
        Ye=feval(IC2_y_fit,x_pos(frame_number),y_pos(frame_number)); % Evaluate obtained polynomial function to find offset
        X2=X-Xe+imagedata.Width/2;
        Y2=Y-Ye+imagedata.Height/2;
        
        % Find resized images without shifts
        IC1f=interp2(X1,Y1,double(IC1(:,:,i)),X,Y);
        IC2f=interp2(X2,Y2,double(IC2(:,:,i)),X,Y);
        
%%%%% RATIO TO TEMP NEEDS TO BE RE_MADE %%%%%%%%%        %calculate temp image from two images
        % Find ratio between camera images intensities
        imratio = IC1f./IC2f;
        imratio(isnan(imratio)) = 0;
        imratio(isinf(imratio)) = 0;
        
        image_temp = interp1(intensity_ratio,REF_TEMPERATURE,imratio);
        image_temp(isnan(image_temp)) = 293;
        
        % Multiply two camera image intensities together and compare against a
        % threshold value. If lower, it means the pixel is noisy and will be
        % registered as 293K
        immultiply = IC1f .* IC2f;
        image_temp(immultiply < DEFAULT_FILTER_THRESHOLD) = NaN;

        % save image to file (append?? - cannot be stored in memory)
        % for now simply append - to 3D array (small enough for sample data
        processed_images(:,:,frame_number) = image_temp;
        
        midspottemp(frame_number) = findMidSpotTemp(image_temp);
        meanmidspottemp(frame_number) = findMeanMidSpotTemp(image_temp);
        peakspottemp(frame_number) = findPeakTemp(image_temp);
        spotradius(frame_number) = findSpotRadius(image_temp);
        
    end
end
end


function midtemp = findMidSpotTemp(imgtemp)
    %midtemp = imgtemp(size(imgtemp,1)/2,size(imgtemp,2)/2);
    midtemp = imgtemp(38,72);
end

% Find spot temperature at middle position
function meanmidtemp = findMeanMidSpotTemp(imgtemp)
    %midtemp = imgtemp(size(imgtemp,1)/2,size(imgtemp,2)/2);
    meanmidtemp = mean(mean(imgtemp(37:39,71:73)));
end

% Find spot size
function spotsize = findSpotSize(immultiply)
%need to check what SIZE_SCALE is
DEFAULT_FILTER_THRESHOLD = 100;
% from IMAGE DATA class
% Temperature image relative to camera image scale
RESIZED_SCALE = 4;
% In micrometer
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
RESIZED_PIXEL_SIZE = ORIGINAL_PIXEL_SIZE / RESIZED_SCALE;
%from SPOT class
SIZE_SCALE = (RESIZED_PIXEL_SIZE * 1e3).^2;
spotsize = sum(sum(immultiply > DEFAULT_FILTER_THRESHOLD)) * SIZE_SCALE;
end

% Find peak temp in image
function peaktemp = findPeakTemp(imgtemp)
peaktemp = max(max(imgtemp));
end


% Find spot size
function spotradius = findSpotRadius(tempimg)
%need to check what SIZE_SCALE is
MELT_THRESHOLD = 1370;
% from IMAGE DATA class
% Temperature image relative to camera image scale
RESIZED_SCALE = 4;
% In micrometer
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
RESIZED_PIXEL_SIZE = ORIGINAL_PIXEL_SIZE / RESIZED_SCALE;
%from SPOT class
SIZE_SCALE = (RESIZED_PIXEL_SIZE * 1e3).^2;

spotarea = sum(sum(tempimg > MELT_THRESHOLD)) * SIZE_SCALE;

spotradius = sqrt(spotarea/pi);
end