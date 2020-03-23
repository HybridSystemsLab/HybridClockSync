%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: run_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all
close all
global flag K v N G L mu sigmaS mag sigmaM H


flag = 0;                                   % Flag to include noise
                                            % Set 'flag = 1;' for
                                            % measurement noise, set 'flag
                                            % = 2;' for noise on the clock
                                            % rate reference, set 'flag =
                                            % 0;' for no noise. 

% System Parameters
K = 0.125;                                  % Controller Jump Gain
H = -1.3;                                   % Controller Flow Gain
v = [0.01,0.1];                             % Communication Interval
G = [0 1 1 0 1; ....                        % Adjacency matrix
     1 0 1 0 0; .... 
     1 0 0 1 0; ....
     0 0 1 0 1; .... 
     1 0 1 1 0];
 
L = diag(sum(G')) - G;                      % Laplacian
N = size(G,2);                              % Number of Agents
mu = 3;                                     % Estimator gain
mag = 1;
sigmaS =  [1;    1;    1;    1;    1;   ];
sigmaM = [0.85 1.15];

run('conditions_script.m')

t0 =      [1;   -1;    2;   -2;    0;   ];
e0 =      [1;   -1;    2;   -2;    0;   ];
u0 =      [0;    0;    0;    0;    0;   ];
eta0 =    [0;   -3;    1;   -4;   -1;   ];
t_star0 = [0;    0;    0;    0;    0;   ];
t_hat0 =  [0;    0;    0;    0;    0;   ];
a0 =      [0.98; 0.99; 1.01; 1.07; 0.87;];
a_hat0 =  [1;    1;    1;    1;    1;   ];
sigma0 = sigmaS;
tau0 = v(1) + (v(2)-v(1))*rand(1);

% state initial condition
x0 = [e0; u0; eta0; t_star0; t_hat0; a0; a_hat0; sigma0; tau0; t0];

% simulation horizon
TSPAN=[0 15];
JSPAN = [0 800];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

% simulate
[t,j,x] = HyEQsolver( @f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

run('post_run_script.m')