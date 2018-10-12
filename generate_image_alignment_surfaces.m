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

%% Input??

%load image data - needed to read mraw (NBytesPerFrame, Pixels, Precision,
% Hight, Width)

% load cih data (camera information header) - creates class: imagedata
imagedata = readcih("example_data");


% FROM asset.m NOT findspotfit.m:
% Specify the precision for reading images
imagedata.Precision = strcat('*ubit',num2str(imagedata.ColorBit));
% Specify the bits per pixel according to image data
bits_per_pixel = imagedata.ColorBit;
% Always a constant
bits_per_byte = 8;
% Specify the total bytes packed in a frame
imagedata.NBytesPerFrame = imagedata.Pixels * bits_per_pixel / bits_per_byte;




%first and last frame to be loaded
frange = [0,1000];

%read in frame from mraw files for each camera
IC1=readmraw(imagedata.folderName,imagedata,frange,1);
IC2=readmraw(imagedata.folderName,imagedata,frange,2);

