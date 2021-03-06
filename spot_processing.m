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
    blockSize = 20;

%% Frame timings
    
    % Calculate time vector for images up to specified end frame
    results.meta.framerate = 1/imagedata.FrameRate;
    t_frame = (0:(end_frame-1)) * results.meta.framerate; %video frame time
    results.t_frame = round(t_frame,9); %round to ns to avoid rounding errors
    %% load Image data for Processing
    % load images in blocks of n, average offset to align first n then apply
    % offset to rest of images assuming it is the same.

    
    resamplelevel = 3;
    
    framesize = imagedata.Width*resamplelevel; %assume square frame
    
    results.pyro = zeros(framesize,framesize,imagedata.TotalFrames);

    results.IC1s = results.pyro;
    results.IC2s = results.pyro;
    results.Otsu_image = results.pyro;
    results.spotradius = zeros(1,imagedata.TotalFrames);
    results.G = zeros(1,imagedata.TotalFrames);
    results.ave = zeros(1,imagedata.TotalFrames);
    results.max = zeros(1,imagedata.TotalFrames);

    frame_number=1;
    
    % create meshgrid
    [X,Y]=meshgrid(1:(framesize),1:(framesize));

    
    % Flag:
    results.too_short=0;
    
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
                %boost resolution
                IC1.resample(:,:,frame_number)=ResBoost(IC1.raw(:,:,i),resamplelevel);
                % find peak, translate image and store offset
                [x1_pos,y1_pos,IC1.x_store,IC1.y_store] = PeakFinder(IC1.resample(:,:,frame_number),IC1.x_store,IC1.y_store);
                
                if isempty(x1_pos) | isempty(y1_pos)
                    results.too_short=1;
                    break
                end
                
                X1=X-x1_pos+framesize/2;
                Y1=Y-y1_pos+framesize/2;
                results.IC1s(:,:,frame_number)=interp2(X1,Y1,double(IC1.resample(:,:,i)),X,Y);

                
                %process camera 2
                %boost resolution
                IC2.resample(:,:,frame_number)=ResBoost(IC2.raw(:,:,i),resamplelevel);
                % find peak, translate image and store offset
                [x2_pos,y2_pos,IC2.x_store,IC2.y_store] = PeakFinder(IC2.resample(:,:,frame_number),IC2.x_store,IC2.y_store);
                
                if isempty(x2_pos) | isempty(y2_pos)
                    results.too_short=1;
                    break
                end
                
                X2=X-x2_pos+framesize/2;
                Y2=Y-y2_pos+framesize/2;
                results.IC2s(:,:,frame_number)=interp2(X2,Y2,double(IC2.resample(:,:,i)),X,Y);
                
                % temperature comparison
                results.pyro(:,:,frame_number) = Pyro(results.IC1s(:,:,frame_number),results.IC2s(:,:,frame_number),intensity_ratio);
                
                %analysis:                
                [results.Otsu_image(:,:,frame_number),results.spotradius(frame_number),results.G(frame_number),results.roundness(frame_number),results.ave(frame_number),results.max(frame_number)]=frame_analysis(results.pyro(:,:,frame_number));
                                
                frame_number = frame_number+1;               
            end
            
            if results.too_short
                break
            end
            
            IC1.x_mean = mean(IC1.x_store);
            X1=X-IC1.x_mean+framesize/2;
            IC1.y_mean = mean(IC1.y_store);
            Y1=Y-IC1.y_mean+framesize/2;
            IC2.x_mean = mean(IC2.x_store);
            X2=X-IC2.x_mean+framesize/2;
            IC2.y_mean = mean(IC2.y_store);
            Y2=Y-IC2.y_mean+framesize/2;
            
        else
            % process frames loop
            for i=1:blockSize
                
                % camera 1
                
                %boost resolution
                IC1.resample(:,:,frame_number)=ResBoost(IC1.raw(:,:,i),resamplelevel);
                results.IC1s(:,:,frame_number)=interp2(X1,Y1,double(IC1.resample(:,:,frame_number)),X,Y);

                % camera 2
                
                %boost resolution
                IC2.resample(:,:,frame_number)=ResBoost(IC2.raw(:,:,i),resamplelevel);
                results.IC2s(:,:,frame_number)=interp2(X2,Y2,double(IC2.resample(:,:,frame_number)),X,Y);
                
                % temperature comparison
                results.pyro(:,:,frame_number) = Pyro(results.IC1s(:,:,frame_number),results.IC2s(:,:,frame_number),intensity_ratio);
                
                %analysis:                
                [results.Otsu_image(:,:,frame_number),results.spotradius(frame_number),results.G(frame_number),results.roundness(frame_number),results.ave(frame_number),results.max(frame_number)]=frame_analysis(results.pyro(:,:,frame_number));

                frame_number = frame_number+1;

            end
        end
    end
    
    % extract R
    results.R = BoundarySpeed(results.spotradius,results.t_frame);
    
end

function [Otsu_image,spotradius,G,roundness,Ave,Max]=frame_analysis(thermal_image)

    MELT_THRESHOLD = 1370+273;
    bw = imbinarize(thermal_image);
    bw = bwareaopen(bw,5);
    cc = bwconncomp(bw,4);
    graindata = regionprops(cc,'Area','Perimeter','BoundingBox');
    [img_size,~] = size(thermal_image);
    idx = 0;
    for i = 1:length(graindata)
        boundingbox = graindata(i).BoundingBox;
        centre = img_size/2;
        if boundingbox(1)<centre & (boundingbox(1)+boundingbox(3))>centre & boundingbox(2)< centre & (boundingbox(2)+boundingbox(4))>centre
            idx = i;
            break
        end
    end
    if idx==0
        Otsu_image = thermal_image;
        spotradius = 0;
        G=0;
        roundness=0;
        Ave = 0;
        Max = 0;
    else
        grain = false(size(bw));
        grain(cc.PixelIdxList{idx}) = true;
        Otsu_image = thermal_image;
        Otsu_image(grain==0)=0;
        % further analysis:
        % measure roundness
        roundness = 4*pi*graindata(idx).Area/graindata(idx).Perimeter^2;
        % spot radius (at melting temp)
        spotradius = findSpotRadius(Otsu_image, MELT_THRESHOLD);
        % temperature gradient (using Phoooper histogram method)
        G = BoundaryGrad(Otsu_image,MELT_THRESHOLD,300);
        % find max of spot
        Max = max(max(Otsu_image));
        Otsu_image(Otsu_image == 0) = nan;
        Ave=mean(mean(Otsu_image,'omitnan'),'omitnan');
    end
end

function [x_pos,y_pos,x_store,y_store] = PeakFinder(img,x_store,y_store)
    % find peaks in image using Fast Peak find and returns the x and y
    % position of the highest one. Append to x and y lists is they are
    % entered.
    ave = mean(mean(img(img~=0)));
    FILTER_THRESHOLD = 25;
    thresh = ave*FILTER_THRESHOLD;
    filt=(fspecial('gaussian', 7,1));
    odds = 1:2:1000;
    intensities = [];

    peaks=FastPeakFind(img,thresh,filt,3,2);
    %assert(length(peaks)>=2,'No peaks found')
    for q = 1:2:size(peaks)
        intensities = [intensities,img(round(peaks(q+1)),round(peaks(q)))];
    end
    true_peak = find(intensities==max(intensities));
    x_pos = peaks(odds(true_peak));
    y_pos = peaks(odds(true_peak)+1);
    
    if(nargin>1)
        x_store = [x_store,x_pos];
        y_store = [y_store,y_pos];
    end
end

function thermal_image = Pyro(img1,img2,intensity_ratio)
    % temperature from ratio and look up table(intensity_ratio and REF
    % TEMPERATURE)
    
    FILTER_THRESHOLD = 300;
    REF_TEMPERATURE = (293:0.5:5000)';
    
    imratio = img1./img2;
    imratio(isnan(imratio)) = 0;
    imratio(isinf(imratio)) = 0;
    immultiply = img1.* img2;
    ave = mean(mean(immultiply));
    imratio(immultiply < FILTER_THRESHOLD) = 0;

    thermal_image = interp1(intensity_ratio,REF_TEMPERATURE,imratio);
    thermal_image(isnan(thermal_image)) = 0;
end

function resample = ResBoost(raw_image,resamplelevel)
    assert(isa(raw_image,'uint16'),'raw image is not type: uint16')
    F = griddedInterpolant(double(raw_image));
    [sx,sy,~] = size(raw_image);
    xq = (0:1/resamplelevel:sx)';
    yq = (0:1/resamplelevel:sy)';
    F.Method = 'linear';
    resample = uint16(F({xq,yq}));
    height = sx*resamplelevel;
    resample = resample(1:height,1:height);
end

% Find spot size
function spotradius = findSpotRadius(tempimg, tempthreshold)
%need to check what SIZE_SCALE is
% from IMAGE DATA class
ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
[img_size,~] = size(tempimg);
RESAMPLED_PIXEL_SIZE = 128/img_size * ORIGINAL_PIXEL_SIZE;
SIZE_SCALE = (RESAMPLED_PIXEL_SIZE)^2;

spotarea = sum(sum(tempimg > tempthreshold)) * SIZE_SCALE;

spotradius = sqrt(spotarea/pi);
end

function G = BoundaryGrad(img,thresh,div)
    [img_size,~] =  size(img);
    ORIGINAL_PIXEL_SIZE = 26.3e-6; % this is actually the source pixel size (not the sensor pixel size)
    RESAMPLED_PIXEL_SIZE = 128/img_size * ORIGINAL_PIXEL_SIZE;
    SIZE_SCALE = (RESAMPLED_PIXEL_SIZE)^2;
    
    img_as_line = reshape(img,[img_size^2,1]);
    temp = sort(img_as_line(~isnan(img_as_line)),'descend');

    x = sqrt([1:length(temp)]*SIZE_SCALE/pi);
    p_index = find(temp>=(thresh-div) & temp<=(thresh+div/2));
    if length(p_index)<=10
        G=0;
    else
        p_boundary = sort(temp(p_index),'descend');
        f=fit(p_boundary,x(p_index).','poly2');
        fx = differentiate(f, 1600);
        G=1/fx;
    end
end

function R = BoundarySpeed(spotradius,t_frame)

    Fs = 1000;
    fc = 10;
    Wn = (2/Fs)*fc;
    order = 20;
    b = fir1(order,Wn,'low',kaiser(order+1,2));
    
    y = filtfilt(b,1,spotradius);

    R = abs(diff(y(2:end))./diff(t_frame));
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

% Find peak temp in image
function peaktemp = findPeakTemp(imgtemp)
peaktemp = max(max(imgtemp));
end
