function [DCM] = quatToDCM(quat)
%QUATTODCM Summary of this function goes here
%   Detailed explanation goes here
w = quat(1);
x = quat(2);
y = quat(3);
z = quat(4);

DCM(1,1) = w^2 + x^2 - y^2 - z^2;
DCM(1,2) = 2 * (x*y - w*z);
DCM(1,3) = 2 * (x*z + y*w);
DCM(2,1) = 2 * (x*y + z*w);
DCM(2,2) = w^2 - x^2 + y^2 - z^2;
DCM(2,3) = 2 * (y*z - x*w);
DCM(3,1) = 2 * (x*z - y*w);
DCM(3,2) = 2 * (y*z + x*w);
DCM(3,3) = w^2 - x^2 - y^2 + z^2;

end

