function [Omega] = skew_quat(gyro)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Form the right quaternion product matrix
Omega = [0, -gyro(1), -gyro(2), -gyro(3); ...
    gyro(1), 0, gyro(3), -gyro(2); ...
    gyro(2), -gyro(3), 0, gyro(1); ...
    gyro(3), gyro(2), -gyro(1), 0];
end

