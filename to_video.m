%a function to write a 3D matrix to a video
%file name should include vidoe format 
% eg 'trial.avi'

%3dmatrix data to be made into a video
% first to axes are frame dimensions, 3rd axis is frame number.
% (pixel x, pixel y, frame number)
function to_video(file_name,matrix)
v = VideoWriter(file_name);
open(v);
%generate data
imagesc(matrix(:,:,1)); 
axis tight manual 
set(gca,'nextplot','replacechildren');
colorbar
%create frames and write to file
for k = 1:length(matrix) 
   imagesc(matrix(:,:,k));
   frame = getframe;
   writeVideo(v,frame);
end
close(v);
%