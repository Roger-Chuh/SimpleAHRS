function [q_new] = quatMult(q1,q2)
%QUATMULT Summary of this function goes here
%   Detailed explanation goes here

q_new = zeros(4,1);

q_new(1) = q1(1) * q2(1) - dot(q1(2:4), q2(2:4));
q_new(2:4) = q1(1)*q2(2:4) + q2(1)*q1(2:4) + cross(q1(2:4), q2(2:4));

end

