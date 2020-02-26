%{
Error state Kalman Filter for attitude estimation
%}

clc; clear all; close all;

addpath('data');
data = importAPDM('data/', 'howard_arm_cal_646.h5', 'SI-000646');
freq = 128; 
dt = 1 / freq;
acc = data(:,1:3);
gyr = data(:,4:6);
truth_q = data(:,10:13); % Output of IMU built in EKF


% Initialize
X = truth_q(1,:)';
delta_theta = [0;0;0];
P = 1e-4 * eye(3); % state covariance corresponding to error quat

est_q(1,:) = X';

noise_gyro = 0.02;
noise_accel = 0.2;

Q = noise_gyro .^2 * eye(3);
R = noise_accel .^2 * eye(3);

gravity = 9.81;

% Run!
N = length(acc);

for i = 2:N
    %% Time Update
    w = gyr(i,:);
    
    % Using q \times q{w} = [w]_R * q
    Omega = skew_quat(w);
    F = eye(4) + (dt/2)*Omega;
     
    X = F*X;
    % Normalize to get nominal state
    X = X ./ norm(X);
    
    %% Error State Update
    q_w = gyroToQuat(w, dt);
    R_w = quatToDCM(q_w);
    
    F_delta = R_w';
    W_delta = eye(3) * dt;
    
    delta_theta = F_delta * delta_theta;
    P = F_delta * P * F_delta' + W_delta * Q * W_delta';
    
    %% Measurement Update
    a = acc(i,:);
    z = a';
    
    H_x = [-X(3),  X(4), -X(1),  X(2); ...
          X(2),  X(1),  X(4),  X(3); ...
          X(1), -X(2), -X(3),  X(4)] * 2 * gravity;
    X_theta = [-X(2), -X(3), -X(4); ...
        X(1), -X(4), X(3); ...
        X(4), X(1), -X(2); ...
        -X(3), X(2), X(1)] * (1/2);
    
    H = H_x * X_theta;
    
    V = eye(3);
    
    S = (H*P*H' + V*R*V');
    K = P*H' / S;
    delta_theta = K*(z - H_x * X);
    P = (eye(3) - K*H)*P*(eye(3) - K*H)' + K * (V*R*V') * K';
    
    delta_q = [1; delta_theta ./ 2];
    
    X = quatMult(X, delta_q);
    
    X = X ./ norm(X);
    
    % EsKF reset
    G = eye(3) - skew(0.5 .* delta_theta);
    P = G*P*G';

    delta_theta = [0;0;0];
    
    est_q(i,:) = X';
end

%% Plotting

eul_est = quat2eul(truth_q);
eul_eskf = quatToEuler(est_q);

figure(1)
plot([eul_est(:,2), eul_eskf(:,2)]); 
legend('kf', 'EsKF')
title("Pitch");

figure(2)
plot([eul_est(:,3), eul_eskf(:,3)]); 
legend('kf', 'EsKF')
title("Roll");