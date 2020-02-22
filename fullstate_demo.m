clc; clear all; close all;

addpath('FullState');
addpath('data');

%% Process Data
data = importAPDM('data\', 'howard_arm_cal_646.h5', 'SI-000646');
freq = 128; 
dt = freq/128;
acc = data(:,1:3);
gyr = data(:,4:6);
mag = data(:,7:9); 
quat_est = data(:,10:13); % Output of IMU built in EKF

%% Tuning Parameters

stationary_range = 1250:2250;
sigma_a = mean(std(acc(stationary_range, :)));
sigma_w = mean(std(gyr(stationary_range, :)));

sigma_p = .5;

Q = sigma_w^2 .* eye(3);
R = sigma_a^2 .* eye(3);
P_init = sigma_p^2 * eye(4);

%% Set up 

q_i_b = [1; 0; 0; 0];
state(1,:) = quatConj(q_i_b)';

EKF = FullStateEKF(Q, R, P_init, q_i_b, dt);

%% Run

for i = 1:length(acc)
   gyro = gyr(i,:);
   %disp("Time update");
   EKF.TimeUpdate(gyro);
   measurement = acc(i,:);
   %disp("Measurement update");
   EKF.MeasurementUpdate(measurement);
   state(i,:) = quatConj(EKF.GetState())';
end

eul_est = quat2eul(quat_est);
eul_acc = quatToEuler(incAccel(acc));
eul_ekf = quatToEuler(state);

%% Plot

figure(1)
plot([eul_acc(:,2),eul_est(:,2), eul_ekf(:,2)]); 
legend('acc','kf', 'me')
title("Pitch");

figure(2)
plot([eul_acc(:,3),eul_est(:,3), eul_ekf(:,3)]); 
legend('acc','kf', 'me')
title("Roll");
