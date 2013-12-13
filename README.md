wheel-soil
==========
The rigid wheel free-body diagram based on the work of Wong and Reece [1], shown in the figure below, is used to model the interaction between the wheel and the soil. Using this model, the drawbar pull D, torque T, and sinkage z, can be estimated for a wheel of weight W, radius r, and wheel width b, travelling at a linear velocity v.

![wheelsoil](https://f.cloud.github.com/assets/5438923/1746579/5dd61358-6444-11e3-9ba1-7be07cad128d.png)

[1] J. Wong and A. Reece, “Prediction of rigid wheel performance based on the analysis of soil-wheel stresses part I. Performance of driven rigid wheels,” J. Terramechanics, vol. 4, no. 1, 1967.

## Example

```javascript
phi = 34.74;       % angle of friction, degrees       
c = 139.28;        % cohesion         
K = 5.15e-4;     % shear modulus
k_eq = 2.54e5;     % pressure-sinkage modulus
n = 1.387;           % exponent of sinkage to width
theta_2 = 0;       % exit angle (assumed)
r = .13;            % wheel radius
b = .16;            % wheel width
W = 130;         % wheel weight

% Coefficients for determining the relative position of max. radial stress
c_1 = 0.43;
c_2 = 0.32;

% Create the WheelSoil object
test = WheelSoil(slip,k_eq,n,c,phi,K,c_1,c_2,W,r,b,theta_2);

guess = 0.4058;
[H,R,D,T,z] = test.calculateWheelPerformance(slip,k_eq,n,c,phi*pi/180,K,c_1,c_2,W,r,b,theta_2*pi/180,guess);
```
