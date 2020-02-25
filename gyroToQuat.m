function q = gyroToQuat(gyro, dt)
%GYROTOQUAT Summary of this function goes here
%   Detailed explanation goes here

if size(gyro,1) == 1
    gyro = gyro';
end

magnitude = norm(gyro);

normalized_gyro = (gyro ./ magnitude);

q = zeros(4,1);

q(1) = cos(magnitude * dt / 2);
q(2:4) = normalized_gyro .* sin(magnitude * dt / 2);

end

