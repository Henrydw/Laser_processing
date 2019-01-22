% Boundry analysis
% a script to iterate though the single exposure data and caluclate:
% G - temp gradient across liquid-solid interface
% R - interface velocity

 run_list = ["100W_400us_1000fps","100W_400us_10000fps","100W_400us_50000fps","100W_400us_100000fps","100W_400us_262500fps",...
    "100W_6400us_1000fps","100W_6400us_10000fps","100W_6400us_50000fps","100W_6400us_100000fps","100W_6400us_262500fps",...
    "200W_100us_1000fps","200W_100us_10000fps","200W_100us_50000fps","200W_100us_100000fps","200W_100us_262500fps",...
    "200W_1600us_1000fps","200W_1600us_10000fps","200W_1600us_50000fps","200W_1600us_100000fps","200W_1600us_262500fps"];

for j=1:length(run_list)
    %% Load info
    folder_name = char(run_list(j));

    folder_path = 'A:\Imperial College London\Hooper, Paul A - spots_v3';

    save_path = 'A:\thermal_highspeed\Data';

    load('intensity_ratio.mat');

    %% process
    
    results=spot_processing(folder_path,folder_name,intensity_ratio);
    
    %% present
    
    %plotting R
    R = figure;
    yyaxis left

    plot(results.t_frame,results.spotradius(2:end), "Color", [0, 0.4470, 0.7410])
    ylabel('Interface radius, r (mm)')
    xlabel('time (ms)')
    title(strcat(folder_name,': Interface radius and speed'))
    %xlim([0.0064,0.0164])

    yyaxis right

    plot(results.t_frame(2:end),results.R)
    ylabel('Interface speed, dr/dt (mm/s)')
    
    %plotting G
    G = figure;
    yyaxis left

    plot(results.t_frame,results.spotradius(2:end), "Color", [0, 0.4470, 0.7410])
    ylabel('Spot radius (mm)')
    xlabel('time (ms)')
    title(strcat(folder_name,': Interface radius and temperature gradient'))
    %xlim([0.0064,0.0164])

    yyaxis right

    plot(results.t_frame,results.G(2:end))
    ylabel('Temperature gradient (at interface), dT/dr (K/mm)')
    
        %simple script to scroll through layers
%     and produce a video of the result

%     vidObj = VideoWriter(strcat('A:\thermal_highspeed\Data\200W_1600us_100000fps','\','tempogram','.avi'));
%     open(vidObj);
% 
%     figure;
%     axis equal;
%     
% for i=1:700
%         cla;
%         hold on;
%         %temp = pixels_by_temperature(:,i);
%         %subplot(2,1,1);
%         %plot(temp,[1:length(temp)]);
%         imagesc(results.Otsu_image(:,:,i));
%         truesize([600 600]);
%         %useful_pixels = length(temp(temp~=0 & temp~=293));
%         %if useful_pixels==0
%         %    continue
%         %end
%         %ylim([0,length(temp(temp~=0 & temp~=293))]);
%         %xlim([293,5000]);
%         %subplot(2,2,1);
%         %imagesc(processed_images(:,:,i));
%         pause(0.01);
%         currFrame = getframe;
%         writeVideo(vidObj,currFrame);
%         hold off;
%     end

    %% save
    
    save_folder = strcat(save_path,'\',folder_name);
    
    if ~exist(save_folder)
        mkdir(save_folder);
    end
    
    final.t_frame = results.t_frame;
    final.R = results.R;
    final.spotradius = results.spotradius;
    final.max = results.max;
    final.ave = results.ave;
    
    save(strcat(save_folder,'\','final.mat'),'final');
    savefig(R,strcat(save_folder,'\','R.fig'));
    savefig(G,strcat(save_folder,'\','G.fig'));
    
    clear results

end




















