function xdot = f(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: f_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system
% Description: Flow map
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

    global mu sigmaS sigmaM H flag

    %e =      [x(1);  x(2);  x(3);  x(4);  x(5)];
    u =      [x(6);  x(7);  x(8);  x(9);  x(10)];
    eta =    [x(11); x(12); x(13); x(14); x(15)];
    t_star = [x(16); x(17); x(18); x(19); x(20)];
    t_hat =  [x(21); x(22); x(23); x(24); x(25)];
    a =      [x(26); x(27); x(28); x(29); x(30)];
    a_hat =  [x(31); x(32); x(33); x(34); x(35)];
    sigma =  [x(36); x(37); x(38); x(39); x(40)];
    %tau = x(41);

    sigma0 =  [sigmaM(1) + (sigmaM(2)-sigmaM(1))*rand(1);...
               sigmaM(1) + (sigmaM(2)-sigmaM(1))*rand(1);...
               sigmaM(1) + (sigmaM(2)-sigmaM(1))*rand(1);...
               sigmaM(1) + (sigmaM(2)-sigmaM(1))*rand(1);...
               sigmaM(1) + (sigmaM(2)-sigmaM(1))*rand(1);   ];
    
    % differential equations
    
    
    if (flag == 0) || (flag == 1)
        xdot = [a + u - sigma + (sigma - sigmaS); zeros(size(u)); H*eta; a; a_hat - (t_hat - t_star); zeros(size(a)); -mu*(t_hat - t_star); zeros(size(sigma)); -1; a + u];
    elseif flag == 2
        %xdot = [a + u - sigma + (sigma - (sigmaS + sigma0)); zeros(size(u)); zeros(size(eta)); a; a_hat - (t_hat - t_star); zeros(size(a)); -mu*(t_hat - t_star); zeros(size(sigma)); -1; a + u];
        xdot = [a + u - sigma + (sigma - sigma0); zeros(size(u)); zeros(size(eta)); a; a_hat - (t_hat - t_star); zeros(size(a)); -mu*(t_hat - t_star); zeros(size(sigma)); -1; a + u];
    elseif flag == 3
        xdot = [a + u - sigma + (sigma - sigmaS); zeros(size(u)); zeros(size(eta)); a; a_hat - (t_hat - (t_star + 0.1*rand(5,1))); zeros(size(a)); -mu*(t_hat - (t_star + 0.1*rand(5,1))); zeros(size(sigma)); -1; a + u];
    end
            
end