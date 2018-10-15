function [Time,Diode,Trig,Sync,x_pos,y_pos,Len] = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [TIME,DIODE,TRIG,SYNC,X_POS,Y_POS] = IMPORTFILE(FILENAME) Reads data
%   from text file FILENAME for the default selection.
%
%   [TIME,DIODE,TRIG,SYNC,X_POS,Y_POS] = IMPORTFILE(FILENAME, STARTROW,
%   ENDROW) Reads data from rows STARTROW through ENDROW of text file
%   FILENAME.
%
% .txt Example:
%   [Time,Diode,Trig,Sync,x_pos,y_pos] = importfile('correction_scan_200W_50us.txt',2, 2000001);

%
% .m Example:
%   [Time,Diode,Trig,Sync,x_pos,y_pos] = importfile('correction_scan_200W_50us.m');
%

%  NB: start and end rows are only needed for .txt files!

%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/08/24 15:57:38

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%check if file is csv (old style NI) or mat (new style picoscope data)
if strfind(filename,'.txt')>0
    
    %% Format string for each line of text:
    %   column1: double (%f)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    %	column6: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%f%f%f%f%f%f%[^\n\r]';
    
    %% Open the text file.
    fileID = fopen(filename,'r');
    
    %% Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1,...
        'Delimiter', delimiter, 'EmptyValue' ,NaN,...
        'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec,...
            endRow(block)-startRow(block)+1,...
            'Delimiter', delimiter, 'EmptyValue',NaN,...
            'HeaderLines', startRow(block)-1,...
            'ReturnOnError', false);
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    %% Close the text file.
    fclose(fileID);
    
%if .mat (picoscope) file
elseif strfind(filename,'.mat')>0
    %make importmatfile funtion??
    ldata=load(filename);
    if endRow==Inf
        endRow=ldata.Length+1; %if endrow not an argument, get from file
    end     
    %calculate time vector
    T=[ldata.Tstart:ldata.Tinterval:(ldata.Tstart+ldata.Tinterval*(ldata.Length-1))]';
    %this dac is too precise, round to 0.1 us
    T=round(T,9);
    
    % take 1 off start/end row because this data doesn't have a header row
    % like the csv data, convert to int type to avoid warnings about floats
    % for indexes
    startRow=int32(startRow-1);
    endRow=int32(endRow-1);
    
    %if sample rate is higher that 1MHz, downsample data to speed up
    %program
    % this is a bit of a bodge, will this cause issues???
    % yup - causing loads of issues.....
% %     if ldata.Tinterval > 1e-6
% %         ds_ratio= floor(1e-6/ldata.Tinterval);
% %         dataArray{:,1}=downsample(T(startRow:endRow*ds_ratio),ds_ratio);
% %         dataArray{:,2}=downsample(ldata.A(startRow:endRow*ds_ratio),ds_ratio);
% %         dataArray{:,6}=downsample(ldata.B(startRow:endRow*ds_ratio),ds_ratio);
% %         dataArray{:,5}=downsample(ldata.C(startRow:endRow*ds_ratio),ds_ratio);
% %     else
% %         dataArray{:,1}=T(startRow:endRow);
% %         %other data channels
% %         dataArray{:,2}=ldata.A(startRow:endRow);
% %         dataArray{:,6}=ldata.B(startRow:endRow);
% %         dataArray{:,5}=ldata.C(startRow:endRow);    
% %     end


    dataArray{:,1}=T(startRow:endRow);
    %other data channels
    dataArray{:,2}=ldata.A(startRow:endRow);
    dataArray{:,6}=ldata.B(startRow:endRow);
    dataArray{:,5}=ldata.C(startRow:endRow);

else
    disp('File format not valid, use .csv of .mat');
end

%% Allocate imported array to column variable names
Time = dataArray{:, 1};
Diode = dataArray{:, 2};
Trig = dataArray{:, 3};
Sync = dataArray{:, 4};
x_pos = dataArray{:, 5};
y_pos = -dataArray{:, 6}; %flip y axis
Len = length(dataArray{:,1});