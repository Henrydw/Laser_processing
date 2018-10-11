function imagedata = readcih(folderName)
% readcih.m
% READCIH Read Photron image data into MATLAB from .cih file
% -------------------------------------------------------------
% Remarks
% This function is to be executed before readmraw function.
% This function has been modified from the original readmraw function to 
% serve the purpose of FYP
% -------------------------------------------------------------
% 
% -------------------------------------------------------------
% Examples
% Load image data
% imageData = readcih(fileID_cih);
% -------------------------------------------------------------
% 
% Author: Adrian Kai Chwen Eu
% Project: Thermal Imaging of metal 3D printing process
% Date: 12/11/16

fileName=strcat(folderName,filesep,'C001H001S0001.cih');
fileID = fopen(fileName);

if fileID < 1
    error('.cih file not found - %s', fileName);
end

% Read necessary information from .cih file
    % Just need to read in information from 1 .cih file since both cih
    % files are the same
    % Algorithm to find the header title and obtain the corresponding value
    % Critical information including
    % a) Total Frame
    % b) Image Width
    % c) Image Height
    % d) FPS
    % e) Color Bit
    headerCount = 0;
    Header=textscan(fileID,'%s','delimiter',':');
    Length = length(Header{1});
    for i = 1:Length
        if ischar(Header{1}{i})
            if strcmpi('Total Frame ',Header{1}{i})
                Total_Frames = str2double(cell2mat(Header{1}(i+1)));
                headerCount = headerCount + 1;
            end
            if strcmpi('Image Width ',Header{1}{i})
                Width = str2double(cell2mat(Header{1}(i+1)));
                headerCount = headerCount + 1;
            end
            if strcmpi('Image Height ',Header{1}{i})
                Height = str2double(cell2mat(Header{1}(i+1)));
                headerCount = headerCount + 1;
            end
            if strcmpi('Record rate(fps) ',Header{1}{i})
                FrameRate = str2double(cell2mat(Header{1}(i+1)));
                headerCount = headerCount + 1;
            end
            if strcmpi('Color Bit ',Header{1}{i})
                ColorBit = str2double(cell2mat(Header{1}(i+1)));
                headerCount = headerCount + 1;
            end
        end
        % If header count reaches 5, break the loop
        if headerCount == 5
            break;
        end
    end


fclose(fileID);


% Assign to corresponding image data
imagedata.Header = Header;
imagedata.Width = Width;
imagedata.Height = Height;
imagedata.Pixels = Width * Height;
imagedata.ColorBit = ColorBit;
imagedata.TotalFrames = Total_Frames;
imagedata.FrameRate = FrameRate;

end