%{
My OO version is killing me so here we are.
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
    
    % State transition matrix
    Omega = skew_quat(w);
    F = eye(4) + dt * Omega / 2;
    
    % Process noise matrix
    W = [-X(2), -X(3), -X(4); ...
          X(1), -X(4),  X(3); ...
          X(4),  X(1), -X(2); ...
         -X(3),  X(2),  X(1)] * (-dt/2);
     
    X = F * X;
    X = X ./ norm(X);
    
    P = F*P*F' + W*Q*W';
    
    %% Measurement Update
    a = acc(i,:);
    % Normalize to match expected value?
    z = (a ./ (norm(a) / gravity))';
    
    H = [-X(3),  X(4), -X(1),  X(2); ...
          X(2),  X(1),  X(4),  X(3); ...
          X(1), -X(2), -X(3),  X(4)] * 2 * gravity;
    
    V = eye(3);
      
    z_pred = H*X;
    
    residual = z - z_pred;
    
    S = H*P*H' + V*R*V';
    
    K = P*H' / S;
    
    X = X + K*residual;
    
    P = (eye(4) - K*H) * P;
    
    X = X / norm(X);
    
    P = (P + P')/2;
    
    est_q(i,:) = X';
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