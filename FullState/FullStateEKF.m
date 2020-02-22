classdef FullStateEKF < handle
    %FULLSTATEEKF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        q_i_b % Rotation from body to inertial
        Q % Process noise covariance
        R % Measurement noise covariance
        P % Estimate error covariance
        dt % Time step
        
    end
    
    methods
        function obj = FullStateEKF(Q, R, P, q_i_b, dt)
            %FULLSTATEEKF Construct an instance of this class
            %   Detailed explanation goes here
            obj.Q = Q;
            obj.R = R;
            obj.P = P;
            obj.q_i_b = q_i_b;
            obj.dt = dt;
        end
        
        function state = GetState(obj)
            % We return state in euler angles
            state = obj.q_i_b;
        end
        
        function [x_hat_prior, P_prior] = TimeUpdate(obj, gyr)
            %disp(["State before update: "]);
            %obj.q_b_i
            
            %disp(["New angular velocity: "]);
            %gyr
            
            F = obj.GetStateTransitionJacobian(gyr);
            %disp(["State transition jacobian: "]);
            %F
            
            W = obj.GetProcessNoiseJacobian(gyr);
            %disp(["Process noise jacobian: "]);
            %W
            
            x_hat_prior = F * obj.q_i_b;
            
            norm_quat = norm(x_hat_prior);
            
            x_hat_prior = x_hat_prior ./ norm_quat;
            %disp("State after time update: ");
            %x_hat_prior
            
            P_prior = F*obj.P*F' + W*obj.Q*W';
            
            obj.q_i_b = x_hat_prior;
            obj.P = P_prior;
        end
        
        function [x_hat_post, P_post] = MeasurementUpdate(obj, accel)
            H = obj.GetMeasurementJacobian(accel);
            V = obj.GetMeasurementNoiseJacobian(accel);
            
            P_prior = obj.P;
            R = obj.R;
            Q = obj.Q;
            x_hat_prior = obj.q_i_b;
            
            q_b_i_prior = quatConj(x_hat_prior); %rotate inertial to body
            
            %disp("Measured acceleration")
            %accel
            
            %disp("Predicted acceleration")
            R_b_i = quatToDCM(q_b_i_prior);
            g_i = [0; 0; 9.81];
            z_pred = R_b_i * g_i;
            %z_pred'
            
            residual = accel' - z_pred
            
            K = (P_prior*H') / (H*P_prior*H' + V*R*V');
            delta_x = K*(residual);
            x_hat_post = x_hat_prior + delta_x;
            P_post = (eye(4) - K*H)*P_prior;
            
            norm_quat = norm(x_hat_post);
            
            x_hat_post = x_hat_post ./ norm_quat;
            
            obj.q_i_b = x_hat_post;
            obj.P = P_post;
            
        end
        
        function [F] = GetStateTransitionJacobian(obj, gyro)
            omega = skew_quat(gyro);
            F = eye(4) + (obj.dt / 2)*omega;
        end
        
        function [H] = GetMeasurementJacobian(obj, accel)
            w = obj.q_i_b(1);
            x = obj.q_i_b(2);
            y = obj.q_i_b(3);
            z = obj.q_i_b(4);
            
            g = 9.81;
            
            H = zeros(3,4);
            H(1,1) = -2*y;
            H(1,2) = 2*z;
            H(1,3) = -2*w;
            H(1,4) = 2*x;
            
            H(2,1) = 2*x;
            H(2,2) = 2*w;
            H(2,3) = 2*z;
            H(2,4) = 2*y;
            
            H(3,1) = 2*w;
            H(3,2) = -2*x;
            H(3,3) = -2*y;
            H(3,4) = 2*z;
            
            H = g.*H;
        end
        
        function [W] = GetProcessNoiseJacobian(obj, gyro)
            q0 = obj.q_i_b(1);
            qv = obj.q_i_b(2:4);
            W = (-obj.dt/2)*[-qv'; (q0 * eye(3) + skew(qv))];
        end
        
        function [V] = GetMeasurementNoiseJacobian(obj, accel)
            V = eye(3);
        end
    end
end

