function [conj] = quatConj(quat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
w = quat(1);
x = quat(2);
y = quat(3);
z = quat(4);

conj = [w; -x; -y; -z];
end

