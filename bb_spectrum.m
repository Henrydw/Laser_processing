function I = bb_spectrum( wl, T)
% bb_spectrum.m
% BB_SPECTRUM calculate black body spectrum intensity vaule at wavelength wl and temperature T 
% function will work with scalar input values, vector values of one input, or 2
% matricies of the same size i.e. output from meshgrid
% -------------------------------------------------------------
% 
% Author: Adrian Kai Chwen Eu
% Project: Thermal Imaging of metal 3D printing process
 
    % Physical constants
    h = 6.626e-34;
    c = 299792458;
    k = 1.381e-23;

    % Planck's law
    % Intensity SI unit: W.sr-1.m-3 where sr is the SI unit and called
    % steradian
    % A steradian can be defined as the solid angle subtended at the center
    % of a unit sphere by a unit area on its surface
    I=2*h*c^2./((wl*10^-9).^5)./(exp((c*h/k)./(wl*10^-9.*T))-1);

end

