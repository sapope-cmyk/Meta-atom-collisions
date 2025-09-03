%% Sets the resonator initial conditions for the simulation of the initial collision and subsequent motion for different resonator parameters

%%% Run "Model_v3.m" first

%% Resonantor frequency for the different simualtions to create the time series plots
% Set 1
Omega = 0.01;

% % Set 2
% Omega = 0.1;

% % Set 3 a and b
% Omega = 10;

% % Set 4
% Omega = 100;

kr = Omega*mr*(2*k/m);
cr = 0;

%% Initial conditions for the outer sphere 1
x1_initial = 0; % Initial displacement for mass 1
dx1_initial = 0.05*sqrt(kg/m); % Initial velocity for mass 1

%% Resonator initial conditions for the different simulations to create the time series plots
% Set 1
AE = 1;
phir = pi;
x1r_initial = sqrt(m*AE/kr)*dx1_initial*sin(phir); % Initial displacement for resonantor mass 1
dx1r_initial = sqrt(m*AE/kr)*dx1_initial*sqrt(kr/mr)*cos(phir); % Initial velocity for resonantor mass 1          

% % Set 2
% AE = 0.75;
% phir = 1.5*pi;
% x1r_initial = sqrt(m*AE/kr)*dx1_initial*sin(phir); % Initial displacement for resonantor mass 1
% dx1r_initial = sqrt(m*AE/kr)*dx1_initial*sqrt(kr/mr)*cos(phir); % Initial velocity for resonantor mass 1 

% % Set 3_a
% AE = 1.5;
% phir = 0*pi;
% x1r_initial = sqrt(m*AE/kr)*dx1_initial*sin(phir); % Initial displacement for resonantor mass 1
% dx1r_initial = sqrt(m*AE/kr)*dx1_initial*sqrt(kr/mr)*cos(phir); % Initial velocity for resonantor mass 1 

% % Set 3_b
% AE = 1.5;
% phir = 0.5*pi;
% x1r_initial = sqrt(m*AE/kr)*dx1_initial*sin(phir); % Initial displacement for resonantor mass 1
% dx1r_initial = sqrt(m*AE/kr)*dx1_initial*sqrt(kr/mr)*cos(phir); % Initial velocity for resonantor mass 1 

% % Set 4
% AE = 2;
% phir = pi;
% x1r_initial = sqrt(m*AE/kr)*dx1_initial*sin(phir); % Initial displacement for resonantor mass 1
% dx1r_initial = sqrt(m*AE/kr)*dx1_initial*sqrt(kr/mr)*cos(phir); % Initial velocity for resonantor mass 1 