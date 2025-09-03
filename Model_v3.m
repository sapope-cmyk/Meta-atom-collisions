%%% Sets up the parameters for the spheres and resonators

%% Constants
g = 9.81; % Acceleration due to gravity

%% Newton's cradle model parameters - Grade 316 steel spheres
rho = 7850; % Sphere density
L = 0.1; % String length
R = 0.01; % Sphere radius
E = 193e9; % Sphere material Youngs modulus
v = 0.27; % Sphere material Poisson's ratio
Rh = 0.1*R; % Shell wall thickness
m = rho*(4/3)*pi*R^3-rho*(4/3)*pi*(R-Rh)^3; % Sphere mass for a shell

kg = m*g/L; % Gravitational spring constant - system modelled for small angle displacements to create a 1D system.
k = E*Rh^2/(2*R); % Sphere spring constant for shell
c = 0; % Damping

mr = m;