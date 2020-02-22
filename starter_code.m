[file,path] = uigetfile('*.h5','hdf5');
data = importAPDM(path,file,'SI-000646'); 

addpath('data');
%% Process Data
freq = 128; 
acc = data(:,1:3);
gyr = data(:,4:6);
mag = data(:,7:9); 
quat = data(:,10:13);

% eul_kf is the output of the IMU's internal kalman filter
eul_kf = quatToEuler(quat);
eul_acc = quatToEuler(incAccel(acc)); 

% Find sigmas for Q and R
stationary_range = 1250:2250;
sigma_a = std(acc(stationary_range, :));
sigma_g = std(gyr(stationary_range, :));


%% Plots
%pitch
figure
plot([eul_acc(:,2),eul_kf(:,2)]); 
legend('acc','kf')
title("Pitch");

%roll
figure
plot([eul_acc(:,3),eul_kf(:,3)]); 
legend('acc','kf')
title("Roll");