function xplus = g(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: g_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
% Description: Jump map
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

    % state
    %x0 = [tauP0; tauO0; tauT0; tauM0; tauS0; q0; M_m0; M_s0; p0];
    global K v N G L flag mag

    e =      [x(1);  x(2);  x(3);  x(4);  x(5)];
    u =      [x(6);  x(7);  x(8);  x(9);  x(10)];
    eta =    [x(11); x(12); x(13); x(14); x(15)];
    t_star = [x(16); x(17); x(18); x(19); x(20)];
    t_hat =  [x(21); x(22); x(23); x(24); x(25)];
    a =      [x(26); x(27); x(28); x(29); x(30)];
    a_hat =  [x(31); x(32); x(33); x(34); x(35)];
    sigma =  [x(36); x(37); x(38); x(39); x(40)];
    tau = x(41);
    t =  [x(42); x(43); x(44); x(45); x(46)];
    
    
    tau_plus = v(1) + (v(2)-v(1))*rand(1);
    
    if flag == 1
        m = mag*rand(N,1);
        u_plus = -K*L*e - a_hat + sigma - K*L*m;
        eta_plus =-K*L*e - K*L*m;
    else
        u_plus = -K*L*e - a_hat + sigma;
        %u_plus = -K*L*e - a_hat;
        eta_plus = -K*L*e;
    end

    xplus = [e; u_plus; eta_plus; t_star; t_hat; a; a_hat; sigma; tau_plus; t];
    
end