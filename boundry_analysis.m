% Boundry analysis
% a script to iterate though the single exposure data and caluclate:
% G - temp gradient across liquid-solid interface
% R - interface velocity

run_list = ["200W_1600us_50000fps"]%["100W_6400us_50000fps","200W_1600us_50000fps",;

for j=1:length(run_list)

    folder_name = char(run_list(j));

    folder_path = 'A:\Imperial College London\Hooper, Paul A - spots_v3';

    save_path = 'A:\thermal_highspeed\Data';

    load('intensity_ratio.mat');

    fits=generate_image_alignment_surfaces(folder_path,folder_name);

    [processed_images, t_frame, spotradius, meanmidspottemp, boundarytempgradient,IC1_full,IC2_full]=process_laser_data(folder_path,folder_name,intensity_ratio,fits);

    graph = plot(t_frame,spotradius);
    
    save(strcat(save_path,'\',folder_name,'\','processed_images.mat'),'processed_images');
    

    %%
    l=length(spotradius);

    images_as_lines = reshape(processed_images,[128*128,l+1]);
    %pixels_by_temperature = zeros([128*128,l+1]);
    G = zeros(1,l+1);

    for i = 1:(l+1)
        line = images_as_lines(:,i);
        if length(line)<=10
            continue
        end
        temp = sort(line(~isnan(line)),'descend');
        %pixels_by_temperature(1:length(temp),i) = temp;
        x = sqrt([1:length(temp)]*SIZE_SCALE/pi);
        p_index = find(temp>=1400 & temp<=temp);
        if length(p_index)<=10
            continue
        end
        p_boundary = sort(temp(p_index),'descend');
        f=fit(p_boundary,x(p_index).','poly2');
        fx = differentiate(f, 1600);
        G(i)=1/fx;
    end
    
    save(strcat(save_path,'\',folder_name,'\','G.mat'),'G');
    
% %     %% Finding G
% %     x = sqrt([1:length(pixels_by_temperature(:,75))]*SIZE_SCALE/pi);
% %     plot(x,pixels_by_temperature(:,75));
% %     
% %     p_index = find(pixels_by_temperature(:,75)>=1400 & pixels_by_temperature(:,75)<=1800);
% %     p_boundary = sort(pixels_by_temperature(p_index,75),'descend');
% %     f=fit(x(p_index).',p_boundary,'poly2');
% %     plot(f,x(p_index),p_boundary)
% %     title('200W 1600us 5kfps frame 75')
% %     legend('raw data','2nd order polynomial')
% %     ylabel('Temperature (K)')
% %     xlabel('Radius m')
% %     
% %     f=fit(p_boundary,x(p_index).','poly2');
% %     fx = differentiate(f, 1600);
% %     fx=1/fx;
% %     
% %     
% % % %     hold on
% % % %     f=fit(x,pixels_by_temperature(:,75),'poly2');
% % % %     plot(f);
% %     
% %     %%
% % 
% %     %simple script to scroll through layers
% %     % and produce a video of the result
% %     %vidObj = VideoWriter(strcat(save_path,'\',folder_name,'\','tempogram','.avi'));
% %     %open(vidObj);
% % 
% %     figure;
% %     axis equal;
% % 
% %     for i=2:length(pixels_by_temperature)
% %         cla;
% %         %hold on;
% %         temp = pixels_by_temperature(:,i);
% %         %subplot(2,1,1);
% %         plot(temp,[1:length(temp)]);
% %         
% %         useful_pixels = length(temp(temp~=0 & temp~=293));
% %         %if useful_pixels==0
% %         %    continue
% %         %end
% %         ylim([0,length(temp(temp~=0 & temp~=293))]);
% %         xlim([293,5000]);
% %         %subplot(2,2,1);
% %         %imagesc(processed_images(:,:,i));
% %         pause(0.5);
% %         %currFrame = getframe;
% %         %writeVideo(vidObj,currFrame);
% %         %hold off;
% %     end
% % 
% %     close(vidObj);
end




















