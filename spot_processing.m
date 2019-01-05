%% spot_processing
% this function takes calibration data and file locations and processes raw
% output into thermal images, as well as doing some minor analysis on each
% frame (eg mean temp, spot radius etc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB: this is for spots only and does not properly calibrate X and Y
% alignment

% This should only be run for SMALL file sizes - this outputs a lot of
% extra data and is NOT OPTIMISED FOR LARGE FILES!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% intensity_ratio - list of wavelength ratios (9415 elements long) 
% corresponding to the REF_TEMPERATURE values, normally 
% load('intensity_ratio.mat') is enough, otherwise use Filter.m to generate
% the list.

% Outputs:

% a single structure:
% results
% which will include pyro images: .pyro
% image frame time: .t_frame
% meta data: .meta : .frame_rate .length and .date

function [results]=spot_processing(folder_path,folder_name,intensity_ratio)
    %% Hard Coded stuff

    DEFAULT_FILTER_THRESHOLD = 50;
    REF_TEMPERATURE = (293:0.5:5000)';

    MELT_THRESHOLD = 1370;

    % from IMAGE DATA class
    % Temperature image relative to camera image scale
    ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
    SIZE_SCALE = (ORIGINAL_PIXEL_SIZE)^2;

    %% Load calibration for MRAWs

    % denote folder with data (should contain .mat *1, .cih *2, .mraw *2)
    
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
    
    %number of frames being processed in total:
    end_frame = imagedata.TotalFrames - 1;

    % number to images to be process per block (all at once crashes computer?)
    blockSize = 50;

%% Frame timings
    
    % Calculate time vector for images up to specified end frame
    results.meta.framerate = 1/imagedata.FrameRate;
    t_frame = (0:(end_frame-1)) * results.meta.framerate; %video frame time
    results.t_frame = round(t_frame,9); %round to ns to avoid rounding errors
    %% load Image data for Processing
    % load images in blocks of n, average offset to align first n then apply
    % offset to rest of images assuming it is the same.

    

    results.pyro = zeros(128,128,imagedata.TotalFrames);

    results.IC1s = results.pyro;
    results.IC2s = results.pyro;

    %midspottemp = zeros(1,end_frame);
    % = zeros(1,end_frame);
    %peakspottemp = zeros(1,end_frame);
    %spotradius = zeros(1,end_frame);
    %boundarytempgradient = zeros(1,end_frame);

    frame_number=1;
    filt=(fspecial('gaussian', 7,1));
    odds = [1:2:1000];
    
    % create meshgrid
    [X,Y]=meshgrid(1:imagedata.Width,1:imagedata.Height);

    % load blocks loop
    for j=1:blockSize:end_frame
        if frame_number+1+blockSize > end_frame
            frange = [frame_number+1,end_frame];
            blockSize=range(frange);
        else
            frange = [frame_number+1,frame_number+1+blockSize];
        end

        %read in frame from mraw files for each camera
        IC1.raw=readmraw(imagedata.folderName,imagedata,frange,1);
        IC1.raw(isnan(IC1.raw)) = 0;
        IC2.raw=readmraw(imagedata.folderName,imagedata,frange,2);
        IC2.raw(isnan(IC1.raw)) = 0;
        if j==1
            % read store average first block
            
            IC1.x_store = [];
            IC1.y_store = [];
            IC2.x_store = [];
            IC2.y_store = [];
            
            for i=1:blockSize
                

                %process camera 1

                intensities = [];
                p1=FastPeakFind(IC1.raw(:,:,frame_number),DEFAULT_FILTER_THRESHOLD,filt,3,2);
                for q = 1:2:size(p1)
                    intensities = [intensities,IC1.raw(round(p1(q+1)),round(p1(q)),frame_number)];
                end
                true_peak = find(intensities==max(intensities));
                p1=p1(odds(true_peak):odds(true_peak)+1);
                IC1.x_store = [IC1.x_store,p1(1)];
                IC1.y_store = [IC1.y_store,p1(2)];
                X1=X-p1(1)+imagedata.Width/2;
                Y1=Y-p1(2)+imagedata.Height/2;
                IC1_new=interp2(X1,Y1,double(IC1.raw(:,:,i)),X,Y);
                results.IC1s(:,:,frame_number) = IC1_new;
                
                %process camera 2
                
                intensities = [];                
                p2=FastPeakFind(IC2.raw(:,:,frame_number),DEFAULT_FILTER_THRESHOLD,filt,3,2);
                for q = 1:2:size(p2)
                    intensities = [intensities,IC2.raw(round(p2(q+1)),round(p2(q)),frame_number)];
                end
                true_peak = find(intensities==max(intensities));
                p2=p2(odds(true_peak):odds(true_peak)+1);
                IC2.x_store = [IC2.x_store,p2(1)];
                IC2.y_store = [IC2.y_store,p2(2)];
                X2=X-p2(1)+imagedata.Width/2;
                Y2=Y-p2(2)+imagedata.Height/2;
                IC2_new=interp2(X2,Y2,double(IC2.raw(:,:,i)),X,Y);
                results.IC2s(:,:,frame_number) = IC2_new;
                
                
                % temperature comparison
                imratio = IC1_new./IC2_new;
                imratio(isnan(imratio)) = 0;
                imratio(isinf(imratio)) = 0;
                immultiply = IC1_new.* IC2_new;
                imratio(immultiply < DEFAULT_FILTER_THRESHOLD) = 0;
                
                image_temp = interp1(intensity_ratio,REF_TEMPERATURE,imratio);
                image_temp(isnan(image_temp)) = 0;
                
                results.pyro(:,:,frame_number) = image_temp;
                
                frame_number = frame_number+1;               
            end
            
            IC1.x_mean = mean(IC1.x_store);
            X1=X-IC1.x_mean+imagedata.Width/2;
            IC1.y_mean = mean(IC1.y_store);
            Y1=Y-IC1.y_mean+imagedata.Width/2;
            IC2.x_mean = mean(IC2.x_store);
            X2=X-IC2.x_mean+imagedata.Width/2;
            IC2.y_mean = mean(IC2.y_store);
            Y2=Y-IC2.y_mean+imagedata.Width/2;
            
        else
            % process frames loop
            for i=1:blockSize
                
                
                % camera 1
                IC1_new=interp2(X1,Y1,double(IC1.raw(:,:,i)),X,Y);
                results.IC1s(:,:,frame_number) = IC1_new;

                % camera 2
                IC2_new=interp2(X2,Y2,double(IC2.raw(:,:,i)),X,Y);
                results.IC2s(:,:,frame_number) = IC2_new;
                
                
                % temperature comparison
                imratio = IC1_new./IC2_new;
                imratio(isnan(imratio)) = 0;
                imratio(isinf(imratio)) = 0;
                immultiply = IC1_new.* IC2_new;
                imratio(immultiply < DEFAULT_FILTER_THRESHOLD) = 0;
                
                image_temp = interp1(intensity_ratio,REF_TEMPERATURE,imratio);
                image_temp(isnan(image_temp)) = 0;
                
                results.pyro(:,:,frame_number) = image_temp;
                
                frame_number = frame_number+1;   
                
            end
        end
            

%         %%%%% RATIO TO TEMP NEEDS TO BE RE_MADE %%%%%%%%%        %calculate temp image from two images
%             % Find ratio between camera images intensities
%         imratio = IC1f./IC2f;
%         imratio(isnan(imratio)) = 0;
%         imratio(isinf(imratio)) = 0;
% 
%         image_temp = interp1(intensity_ratio,REF_TEMPERATURE,imratio);
%         image_temp(isnan(image_temp)) = 293;
% 
%         % Multiply two camera image intensities together and compare against a
%         % threshold value. If lower, it means the pixel is noisy and will be
%         % registered as 293K
%         immultiply = IC1f .* IC2f;
%         %image_temp(immultiply < DEFAULT_FILTER_THRESHOLD) = NaN;
% 
%         % save image to file (append?? - cannot be stored in memory)
%         % for now simply append - to 3D array (small enough for sample data
%         results.processed_images(:,:,frame_number) = image_temp;
% 
%         %midspottemp(frame_number) = findMidSpotTemp(image_temp);
%         %meanmidspottemp(frame_number) = findMeanMidSpotTemp(image_temp);
%         %peakspottemp(frame_number) = findPeakTemp(image_temp);
%         %spotradius(frame_number) = findSpotRadius(image_temp, MELT_THRESHOLD);
%         %temp_diff = 150;    
%         %dr = findSpotRadius(image_temp, MELT_THRESHOLD-150)-findSpotRadius(image_temp, MELT_THRESHOLD+150);
%         %boundarytempgradient(frame_number) = (2*300)/dr;
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
% In meters
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
SIZE_SCALE = (ORIGINAL_PIXEL_SIZE)^2;
spotsize = sum(sum(immultiply > DEFAULT_FILTER_THRESHOLD)) * SIZE_SCALE;
end

% Find peak temp in image
function peaktemp = findPeakTemp(imgtemp)
peaktemp = max(max(imgtemp));
end


% Find spot size
function spotradius = findSpotRadius(tempimg, tempthreshold)
%need to check what SIZE_SCALE is
% from IMAGE DATA class
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
SIZE_SCALE = (ORIGINAL_PIXEL_SIZE)^2;

spotarea = sum(sum(tempimg > tempthreshold)) * SIZE_SCALE;

spotradius = sqrt(spotarea/pi);
end