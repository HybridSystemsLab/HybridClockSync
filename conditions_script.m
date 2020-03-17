[T,Q] = eig(L,'vector');
indx = find(abs(Q) < .0001);
T(:,[1,indx]) = T(:,[indx,1]); % If indx neq 0, then it will switch columns

F = T^-1*L*T;

barL = F(2:end,2:end);

Z = zeros(length(barL));
I = eye(length(barL));

Af = [Z I; Z H*I];
Ag = [I Z; -K*barL Z];


A1low = expm(Af*v(1))*Ag;
A1high = expm(Af*v(2))*Ag;
eigA1high = max(abs(eig(A1high)))
eigA1low = max(abs(eig(A1low)))

if eigA1high > 1
    error('eigA1high > 1')
end

Af3 = [0 mu; -1 -1;];
Af4 = [zeros(N-1) mu*eye(N-1); -eye(N-1) -eye(N-1);];

cvx_begin sdp
    variable P2(size(Af3)) symmetric
    Af3'*P2 + P2*Af3 <= -zeros(size(Af3));
    P2 >= eye(size(Af3));
cvx_end

cvx_begin sdp
    variable P3(size(Af4)) symmetric
    variable Q2(size(Af4)) symmetric
    Af4'*P3 + P3*Af4 <= -1*Q2;
    P3 >= eye(size(Af4));
    Q2 >= eye(size(Af4));
cvx_end

cvx_begin sdp 
    variable P1(size(Ag)) symmetric 
    variable Q1(size(Ag)) symmetric 
    minimize (0)
    subject to
    A1high'*P1*A1high - P1 <= -0*Q1;
    A1low'*P1*A1low - P1 <= -0*Q1;
    P1 >= eye(size(Ag));
    Q1 >= eye(size(Ag));
cvx_end

beta1 = norm(P2*Af3 + Af3'*P2);
beta2 = norm(P3*Af4 + Af4'*P3);

%%
N = 5;
T1 = v(1);
T2 = v(2);

[T,Q] = eig(L,'vector');
indx = find(abs(Q) < 0.0001);
T(:,[1,indx]) = T(:,[indx,1]);
F = (T^-1)*L*T;
barL = F(2:end,2:end);

Af1 = [0 1; 0 H];
Af2 = [zeros(N-1) eye(N-1); zeros(N-1) H*eye(N-1)];
Af3 = [0 mu; -1 -1;];
Af4 = [zeros(N-1) mu*eye(N-1); -eye(N-1) -eye(N-1);];
Ag1 = [1 0; 0 0];
Ag2 = [eye(N-1) zeros(N-1); -K*barL zeros(N-1)];
Bf1 = [1 0];
Bf2 = [eye(N-1) zeros(N-1)];
  
v_n=0:0.01:T2;
lambda_min = zeros(size(v_n));
lambda_max = zeros(size(v_n));
e2hv = zeros(size(v_n));
a2 = zeros(size(v_n));
for i=1:max(size(v_n))
tt=v_n(i);
Matrix=expm(Af2'*tt)*P1*expm(Af2*tt);
a_2(i)=max(sqrt(eig(Matrix'*Matrix)));
lambda_min(i)=min(sqrt(eig(Matrix'*Matrix)));
lambda_max(i)=max(sqrt(eig(Matrix'*Matrix)));
e2hv(i) = exp(2*H*tt);
end

a2 = max([max(e2hv) max(lambda_max) max(eig(P2)) max(eig(P3))]);
a2_alt = max([max(e2hv) max(lambda_max) max(sqrt(eig(P2'*P2))) max(sqrt(eig(P3'*P3)))]);
a1 = min([min(e2hv) min(lambda_min) min(eig(P2)) min(eig(P3))]);

v_n=T1:0.01:T2;
lambda_min = zeros(size(v_n));
lambda_max = zeros(size(v_n));
for i=1:max(size(v_n))
tt=v_n(i);
Matrix=Ag2'*expm(Af2'*tt)*P1*expm(Af2*tt)*Ag2 - P1;
lambda_min(i)=min(eig(Matrix));
lambda_max(i)=max(eig(Matrix));
end

k1 = 2*max(a_2)*norm(Bf2);

k2_rge = [0,-min(lambda_min)]

epsilon1 = min(roots([k1 -2*beta2  -k1]));
epsilon2 = max(roots([k1 -2*beta2  -k1]));

% if (epsilon1 < 0) && (epsilon2 < 0)    
%     epsilon1 = eps;
%     epsilon2 = eps;
% end
% 
% 
% if eps > 0
%     epsilon1 = eps;
%     epsilon2 = eps;
% end

a2bar = max([1 max(eig(P1))]);

kbar1 = max([k1/(2*epsilon1) ((k1*epsilon1)/2 - beta2)]);
kbar2 = max([k1/(2*epsilon2) ((k1*epsilon2)/2 - beta2)]);
kbar = max([kbar1 kbar2]);

a1_til = a2bar/2;
a2_til = a2bar/2;

k2 = max(k2_rge);

display('Condition:')
exp(T2*kbar/a2)*(1 - k2/a2)

display('T2 <')
(log(1/(1 - k2/a2)))*a2/kbar

display('Condition alt:')
exp(T2*kbar/a2)*(1 - k2/a2_alt)

if exp(T2*kbar/a2)*(1 - k2/a2) > 1    
    error('exp(T2*kbar)*(1 - k2/a2) > 1')
end
