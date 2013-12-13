function [H, R, D,T,z,D_unc,T_unc,z_unc] = WheelSoilExampleDriver(slip,weight)
    % Example code for using the WheelSoil class
    %   Input:  wheel slip
    %           weight of the wheel
    %   Output: drawbar pull, D
    %           torque, T
    %           maximum sinkage, z
    %           drawbar pull uncertainty, D_unc
    %           torque uncertainty, T_unc
    %           sinkage uncertainty, z_unc

%     % Wheel/Soil Parameters
%     phi = 33.083;       % angle of friction, degrees       
%     c = 714.971;        % cohesion         
%     K = 5.571e-4;     % shear modulus
%     k_eq = 0.141e6;     % pressure-sinkage modulus
%     n = 0.89;           % exponent of sinkage to width
%     theta_2 = -5;       % exit angle (assumed)
%     r = .13;            % wheel radius
%     b = .16;            % wheel width
%     W = weight;         % wheel weight
    
    phi = 34.74;       % angle of friction, degrees       
    c = 139.28;        % cohesion         
    K = 5.15e-4;     % shear modulus
    k_eq = 2.54e5;     % pressure-sinkage modulus
    n = 1.387;           % exponent of sinkage to width
    theta_2 = 0;       % exit angle (assumed)
    r = .13;            % wheel radius
    b = .16;            % wheel width
    W = weight;         % wheel weight
    
    % Coefficients for determining the relative position of max. radial stress
    c_1 = 0.43;
    c_2 = 0.32;
    
%         % Coefficients for determining the relative position of max. radial stress
%     c_1 = 0.5;
%     c_2 = 0.5;

    %slip = 0.7;
    guess = 0.4058;

    % soil parameter variability
    phi_unc = 2.875*pi/180;
    c_unc = 415.491;
    K_unc = 7.525e-5;
    k_eq_unc = 13000;
    n_unc = 0.14;
    theta_2_unc = 0.03;
    W_unc = weight*.05;
    c_1_unc = 0.2;
    c_2_unc = 0.2;
    
    D_unc = 0;
    T_unc = 0;
    z_unc = 0;

    % Create the WheelSoil object
    test = WheelSoil(slip,k_eq,n,c,phi,K,c_1,c_2,W,r,b,theta_2);
    
%     % Calculate the uncertainty!
%     %   Need the partial derivatives w.r.t. soil parameters
% 
%     % Calculate the partial derivatives w.r.t. pressure sinkage modulus
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq+sqrt(-1)*1e-10,...
%         n,c,phi*pi/180,K,c_1,c_2,W,r,b,theta_2*pi/180,guess);
%     dDdk_eq = imag(D)/1e-10;
%     dTdk_eq = imag(T)/1e-10;
%     dzdk_eq = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. exponent of sinkage to width
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,...
%         n+sqrt(-1)*1e-10,c,phi*pi/180,K,c_1,c_2,W,r,b,theta_2*pi/180,guess);
%     dDdn = imag(D)/1e-10;
%     dTdn = imag(T)/1e-10;
%     dzdn = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. cohesion
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,...
%         c+sqrt(-1)*1e-10,phi*pi/180,K,c_1,c_2,W,r,b,theta_2*pi/180,guess);
%     dDdc = imag(D)/1e-10;
%     dTdc = imag(T)/1e-10;
%     dzdc = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. angle of friction
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
%         phi*pi/180+sqrt(-1)*1e-10,K,c_1,c_2,W,r,b,theta_2*pi/180,guess);
%     dDdphi = imag(D)/1e-10;
%     dTdphi = imag(T)/1e-10;
%     dzdphi = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. shear modulus
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
%         phi*pi/180,K+sqrt(-1)*1e-10,c_1,c_2,W,r,b,theta_2*pi/180,guess);
%     dDdK = imag(D)/1e-10;
%     dTdK = imag(T)/1e-10;
%     dzdK = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. constant 1
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
%         phi*pi/180,K,c_1+sqrt(-1)*1e-10,c_2,W,r,b,theta_2*pi/180,guess);
%     dDdc_1 = imag(D)/1e-10;
%     dTdc_1 = imag(T)/1e-10;
%     dzdc_1 = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. constant 2
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
%         phi*pi/180,K,c_1,c_2+sqrt(-1)*1e-10,W,r,b,theta_2*pi/180,guess);
%     dDdc_2 = imag(D)/1e-10;
%     dTdc_2 = imag(T)/1e-10;
%     dzdc_2 = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. wheel weight
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
%         phi*pi/180,K,c_1,c_2,W+sqrt(-1)*1e-10,r,b,theta_2*pi/180,guess);
%     dDdW = imag(D)/1e-10;
%     dTdW = imag(T)/1e-10;
%     dzdW = imag(z)/1e-10;
% 
%     % Calculate the partial derivatives w.r.t. exit angle
%     [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
%         phi*pi/180,K,c_1,c_2,W,r,b,theta_2*pi/180+sqrt(-1)*1e-10,guess);
%     dDdtheta_2 = imag(D)/1e-10;
%     dTdtheta_2 = imag(T)/1e-10;
%     dzdtheta_2 = imag(z)/1e-10;
% 
%     % Calculate the total uncertainty in wheel performance
%     D_unc = sqrt(dDdk_eq^2*k_eq_unc^2+dDdn^2*n_unc^2+dDdc^2*c_unc^2+...
%         dDdphi^2*phi_unc^2+dDdK^2*K_unc^2+dDdc_1^2*c_1_unc^2+...
%         dDdc_2^2*c_2_unc^2+dDdW^2*W_unc^2+dDdtheta_2^2*theta_2_unc^2);
%     T_unc = sqrt(dTdk_eq^2*k_eq_unc^2+dTdn^2*n_unc^2+dTdc^2*c_unc^2+...
%         dTdphi^2*phi_unc^2+dTdK^2*K_unc^2+dTdc_1^2*c_1_unc^2+...
%         dTdc_2^2*c_2_unc^2+dTdW^2*W_unc^2+dTdtheta_2^2*theta_2_unc^2);
%     z_unc = sqrt(dzdk_eq^2*k_eq_unc^2+dzdn^2*n_unc^2+dzdc^2*c_unc^2+...
%         dzdphi^2*phi_unc^2+dzdK^2*K_unc^2+dzdc_1^2*c_1_unc^2+...
%         dzdc_2^2*c_2_unc^2+dzdW^2*W_unc^2+dzdtheta_2^2*theta_2_unc^2);

    % Calculate the wheel performance
    [H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,...
        phi*pi/180,K,c_1,c_2,W,r,b,theta_2*pi/180,guess);
    z = z*1000;
end