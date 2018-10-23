%a function to write a 3D matrix to a video
%file name should include vidoe format 
% eg 'trial.avi'

%3dmatrix data to be made into a video
% first to axes are frame dimensions, 3rd axis is frame number.
% (pixel x, pixel y, frame number)
function to_video()
v = VideoWriter('peaks.avi');
open(v);
%generate data
Z = peaks;
surf(Z); 
axis tight manual 
set(gca,'nextplot','replacechildren');
%create frames and write to file
for k = 1:20 
   surf(sin(2*pi*k/20)*Z,Z)
   frame = getframe;
   writeVideo(v,frame);
end
close(v);
%