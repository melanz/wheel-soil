classdef WheelSoil < handle
    % Encapsulates Bekker-Wong Theory for Vehicle-Wheel Interaction
    %   Input: Soil parameters, Output: Wheel performance
    
    properties
        % properties with temporary values
        phi_degree = 19.735;        % angle of internal shearing resistance
        phi = 0;
        c = 298.123;                % cohesion, Pa
        K = 0.00068509;             % shear deformation modulus, m
        k_1 = 0.141e6;              % pressure sinkage modulus 1
        k_2 = 0;                    % pressure sinkage modulus 2
        n = 0.89;                   % exponent of sinkage to width ratio
        gamma = 0.0575;             % density, N/m^3
        
        % Wheel properties
        r = .13;                    % radius, m
        b = .16;                    % width, m
        W = 135;                    % vertical axle load, N
        
        % Coefficients for the relative position of max. radial stress
        c_1 = 0.5;                  % 0.43;
        c_2 = 0.5;                  % 0.32;
        
        slip = 0.7;                 % slip
        
        % Important angles along the wheel
        theta_1 = 0.1;
        theta_2 = 0;                % assume no t
        theta_m = 0;
        
        z_0 = 0;                    % maximum sinkage, m
    end
    
    methods
        % Constructor for the WheelSoil class (contains the variables)
        function TS = WheelSoil(i,keq,n,c,phi,K,c1,c2,Fz,r,b,theta2)
            TS.phi_degree = phi;    % angle of internal shearing resistance
            TS.phi = 0;
            TS.c = c;               % cohesion, Pa
            TS.K = K;               % shear deformation modulus, m
            
            TS.k_1 = keq;           % pressure sinkage modulus 1
            TS.k_2 = 0;             % pressure sinkage modulus 2
            TS.n = n;               % exponent of sinkage to width ratio
            
            TS.gamma = 0.048;       % density, N/m^3 (doesn't matter)
            
            % Wheel properties
            TS.r = r;               % radius, m
            TS.b = b;               % width, m
            TS.W = Fz;              % vertical axle load, N
            
            % Coefficients for the relative position of max. radial stress
            TS.c_1 = c1;
            TS.c_2 = c2;
            
            TS.slip = i;            % slip
            
            TS.theta_1 = 0.1;
            TS.theta_2 = theta2*pi/180;    % assume no t
            TS.theta_m = 0;
            
            TS.z_0 = 0;
            TS.slip = i;
            TS.phi = TS.phi_degree*pi/180;
        end
        
        % Get the normal stress along the front of the wheel for a given
        % angle
        function sigma_1 = getSigma1(TS,theta)
            z = (cos(theta)-cos(TS.theta_1))*TS.r;
            sigma_1 = TS.k_1.*(z./TS.b).^TS.n; % k_1 is k_eq here
        end
        
        % Get the normal stress along the back of the wheel for a given
        % angle
        function sigma_2 = getSigma2(TS,theta)
            z = (cos(TS.theta_1-((theta-TS.theta_2)/...
                (TS.theta_m-TS.theta_2))*(TS.theta_1-TS.theta_m))...
                -cos(TS.theta_1))*TS.r;
            sigma_2 = TS.k_1.*(z./TS.b).^TS.n; % k_1 is k_eq here
        end
        
        % Get the normal stress at any angle on the wheel
        function sigma = getSigma(TS,theta)
            sigma = zeros(length(theta),1);
            for i=1:length(theta) 
                %sigma(i) = getSigma1(TS,theta(i));
                if(theta(i)>=TS.theta_m && theta(i)<=TS.theta_1)
                    sigma(i) = getSigma1(TS,theta(i));
                elseif(theta(i)>=TS.theta_2 && theta(i)<TS.theta_m)
                    sigma(i) = getSigma2(TS,theta(i));
                else
                    sigma(i) = 0;
                end
            end
            sigma = sigma';
        end
        
        % Get the shear stress at any angle on the wheel
        function tau = getTau(TS,theta)
            tau = zeros(length(theta),1);
            for i=1:length(theta) 
                if(theta(i)>=TS.theta_2 && theta(i)<=TS.theta_1)
                    sigma = TS.getSigma(theta(i));
                    j = ((TS.theta_1-theta(i))-(1-TS.slip)*(sin(TS.theta_1)...
                        -sin(theta(i))))*TS.r;
                    tau(i) = (TS.c+sigma*tan(TS.phi))*(1-exp(-j/TS.K));
                else
                    tau(i) = 0;
                end
            end
            tau = tau';
        end
        
        % Error function, including derivative (used for optimization)
        function [error,derrordthetac] = errorFunction(TS,theta_c)
            error = TS.getError(theta_c);
            if (nargout > 1)
                % use complex differentiation
                derrordthetac = imag(TS.getError(theta_c+sqrt(-1)*1e-10))/1e-10;
            end
        end
        
        % Calculate the normalized error in weight (used for optimization)
        function error = getError(TS,theta_c)
            W_prime = TS.getWeight(theta_c);
            error =  (W_prime-TS.W)/TS.W;
        end
        
        % Get the weight for a given entrance angle
        function Wcheck = getWeight(TS,theta_c)
            radius = TS.r;
            width = TS.b;
            TS.theta_1 = theta_c;
            TS.theta_m = (TS.c_1+TS.c_2*TS.slip)*theta_c;
            
            % TEST
            %Wcheck = 2*integral(@(theta)(TS.getSigma(theta).*cos(theta)+TS.getTau(theta).*sin(theta)),0,TS.theta_1).*radius.*width;
            % END TEST

            Wcheck = integral(@(theta)(TS.getSigma(theta).*cos(theta)),...
                TS.theta_2,TS.theta_1);
            Wcheck = Wcheck + integral(@(theta)(TS.getTau(theta).*...
                sin(theta)),TS.theta_2,TS.theta_1);
            Wcheck = Wcheck.*radius.*width;
        end
        
        % Calculate the wheel performance for a given set of soil
        % parameters
        function [H,R,D,T,z] = calculateWheelPerformance(TS,slip,keq,n,c,...
                phi,K,c1,c2,Fz,r,b,theta2,guess)
            TS.slip = slip;
            TS.k_1 = keq;
            TS.n = n;
            TS.c = c;
            TS.phi = phi;
            TS.K = K;
            TS.c_1 = c1;
            TS.c_2 = c2;
            TS.W = Fz;
            TS.r = r;
            TS.b = b;
            TS.theta_2 = theta2;

            radius = r;
            width = b;
            
            % Calculate the entrance angle (assuming an exit angle)
            options = optimset('Jacobian','off','Display','on');
            f = @(theta_c) TS.errorFunction(theta_c);
            theta_c = fsolve(f,guess,options);
            TS.theta_1 = theta_c;
            TS.theta_m = (TS.c_1+TS.c_2*TS.slip)*TS.theta_1;
            TS.z_0 = (1-cos(TS.theta_1))*TS.r;

            % Calculate the total thrust
            H = integral(@(theta)(TS.getTau(theta).*cos(theta)),...
                TS.theta_2,TS.theta_1);
            H = H.*radius.*width;

            % Calculate the total motion resistance
            R = integral(@(theta)(TS.getSigma(theta).*sin(theta)),...
                TS.theta_2,TS.theta_1);
            R = R.*radius.*width;
            
            % Calculate the total drawbar pull
            D = H-R;
            
            % Calculate the total torque
            T = integral(@(theta)(TS.getTau(theta)),TS.theta_2,TS.theta_1);
            T = T.*(radius^2*width);
            
            %Calculate the maximum sinkage
            z = (1-cos(TS.theta_1)).*radius;
        end
    end
end

