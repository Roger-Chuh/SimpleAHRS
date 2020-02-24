%{
Error state Kalman Filter for attitude estimation
%}

clc; clear all; close all;

addpath('data');
data = importAPDM('data/', 'howard_arm_cal_646.h5', 'SI-000646');
freq = 128; 
dt = freq/128;
acc = data(:,1:3);
gyr = data(:,4:6);
truth_q = data(:,10:13); % Output of IMU built in EKF


% Initialize
X = truth_q(1,:)';
P = 1e-4 * eye(4);

est_q(1,:) = X';

noise_gyro = 0.02;
noise_accel = 0.15;

Q = noise_gyro .^2 * eye(3);
R = noise_accel .^2 * eye(3);

gravity = 9.81;

% Run!
N = length(acc);

for i = 2:N
    %% Time Update
    w = gyr(i,:);
    
    
    %% Measurement Update
    a = acc(i,:);

    
end

%% Plotting

eul_est = quat2eul(est_q);
eul_ekf = quatToEuler(truth_q);

figure(1)
plot([eul_est(:,2), eul_ekf(:,2)]); 
legend('kf', 'me')
title("Pitch");

figure(2)
plot([eul_est(:,3), eul_ekf(:,3)]); 
legend('kf', 'me')
title("Roll");