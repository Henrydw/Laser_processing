function ImageMatrix = readmraw(folderName,ImageData,numimgs,cameraNum)
% readmraw.m
% READMRAW Read monchrome Photron images with any bit depth into a m x n
% array according to the pixel width and pixel height
% -------------------------------------------------------------
% Remarks
% Please ensure that readcih(fileID) is being executed first to obtained
% the appropriate image data to be used in readmraw function
% .mraw files cannot be read using imread function. This function serves to
% read .mraw files using low level function fread
% This function has been modified from the original readmraw function to 
% serve the purpose of FYP
% -------------------------------------------------------------
% 
% -------------------------------------------------------------
% Example
% Load frame 7 from camera 2
% I=readmraw(fileID_mraw,ImageData_from_cih,7,2);
% Load frame 1200 from camera 1
% I=readmraw(fileID_mraw,ImageData_from_cih,1200,1);
% -------------------------------------------------------------
%
% Author: Adrian Kai Chwen Eu
% Project: Thermal Imaging of metal 3D printing process
% Date created: 01/11/16
% Latest update: 12/11/16

fileName=[folderName filesep 'C00' int2str(cameraNum) 'H001S0001.mraw'];
fileID = fopen(fileName);

% If no file opened, set image array value to 0
if fileID < 1
    ImageMatrix = 0;
    error('.mraw file not found - %s',fileName);
elseif fileID >= 1
    % Load only one frame
    if (length(numimgs) == 1)
        first_frame = numimgs;
        frames = 1;
    % Load a specified range of frames
    else
        first_frame = numimgs(1,1);
        last_frame = numimgs(1,2);
        frames = last_frame - first_frame + 1;
    end
    % Reset to 1st byte position in file
    %frewind(fileID);
    % Offset to specified frame from begingin of file
    fseek(fileID,(first_frame - 1) * ImageData.NBytesPerFrame, 'bof');
    % Initialise array for image values
    %I = zeros(ImageData.Pixels,frames,ImageData.ArrayInitPrecision);
    % Load image frame by frame into array
    % Bits ordering needs to be big-endians to allow bytes to be read in correct order   

    I = fread(fileID,ImageData.Pixels*frames,ImageData.Precision,'b');
  
    %I=gpuArray(I);
    
    % Algorithm to reshape the 1 x k array to m x n array
    ImageMatrix = permute(reshape(I,[ImageData.Width ImageData.Height frames]),[2 1 3]);
    
    % If it is first camera, just rotate then vert flip to align image XY with scanhead XY
    if cameraNum == 1
           ImageMatrix(:,:,1:1:frames)=flipud(imrotate(ImageMatrix(:,:,1:1:frames),-22,'bilinear','crop'));
    end
    
    % If it is second camera, flip the image horizontally
    % (extra mirror) then rotate and vert flip
    if cameraNum == 2
        ImageMatrix(:,:,1:1:frames)=flipud(imrotate(fliplr(ImageMatrix(:,:,1:1:frames)),-22,'bilinear','crop'));
    end
end

%ImageMatrix=gather(ImageMatrix);

fclose(fileID);

end



