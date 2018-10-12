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

%% Inputs:

%-----------------------
%load image data - needed to read mraw (NBytesPerFrame, Pixels, Precision,
% Hight, Width)

% data folder (should contain .mat *1, .cih *2, .mraw *2)
data_folder = "example_data";


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

%% --------------------------- image processing

% ----------------------- load mraw

%number of frames being processed in total:
end_frame = imagedata.TotalFrames - 1;

% number to images to be process per block (all at once crashes computer?)
blockSize = end_frame;
% tnmwfi!!!!

%first and last frame to be loaded
frange = [0,blockSize];

%read in frame from mraw files for each camera
IC1=readmraw(imagedata.folderName,imagedata,frange,1);
IC2=readmraw(imagedata.folderName,imagedata,frange,2);

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

%----------------------- finding peaks

% frame number
frame_number=1;

for i=1:blockSize
    
    frame_number = frame_number+1;
    
    % Calculate threshold params for peak finder function
    thresh=(max([min(max(IC1(:,:,i),[],1))  min(max(IC1(:,:,i),[],2))]));
    %find some peaks for cam 1
    p = FastPeakFind(IC1(:,:,i),thresh,filt,3,2);
    % Check to see if peaks were found
    % Obtain peak in terms of pixel value (not real size)
    % Save only coords of first found peak
    
%     % *
%     % if zero, save as NaN - allows frames which do not have laser in them
%     % to be detected - surely this is already done by the Diode data???
    %if zero, save as nan
    if (size(p) > 0)
        IC1_p(:,frame_number) = p(1:2);
    else
        IC1_p(:,frame_number) = [NaN NaN];
    end

    % Repeat above but for other camera
    thresh = (max([min(max(IC2(:,:,i),[],1))  min(max(IC2(:,:,i),[],2))]));
    p = FastPeakFind(IC2(:,:,i),thresh,filt,3,2);
%     % *
    if (size(p) > 0)
        IC2_p(:,frame_number) = p(1:2);
    else
        IC2_p(:,frame_number) = [NaN NaN];
    end
end

%% ---------------------------------- diode processing

% Read laser data file
[t_daq,Diode,~,~,x,y,~] = importfile(laserdatafn,laserdata.startrow,laserdata.endrow);








