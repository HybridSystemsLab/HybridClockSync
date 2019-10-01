clear N T1 T2 H mu gamma G L T Q indx F barL sum

N = 5;
T1 = 0.5;
T2 = 1.5;
H = 0;
mu = 400;
gamma = 0.3;
G = [0 1 1 0 1; 1 0 1 0 0; 1 0 0 1 0; ....
    0 0 1 0 1; 1 0 1 1 0];                         % Adjacency matrix
L = diag(sum(G')) - G;
[T,Q] = eig(L,'vector');
indx = find(abs(Q) < 0.0001);
T(:,[1,indx]) = T(:,[indx,1]);
F = (T^-1)*L*T;
barL = F(2:end,2:end);

Af1 = [0 1; 0 H];
Af2 = [zeros(N-1) eye(N-1); zeros(N-1) H*eye(N-1)];
Af3 = [0 -mu; 1 -1;];
Af4 = [zeros(N-1) -mu*eye(N-1); eye(N-1) -eye(N-1);];
Ag1 = [1 0; 0 0];
Ag2 = [eye(N-1) zeros(N-1); -gamma*barL zeros(N-1)];
Bf1 = [1 0];
Bf2 = [eye(N-1) zeros(N-1)];

P1 = [33.6103    0.0000   -0.0000    0.0000    4.2030   -0.0000    0.0000   -0.0000;...
       0.0000   28.6082    0.0000    0.0000   -0.0000    5.7271   -0.0000    0.0000;...
      -0.0000    0.0000   25.3450   -0.0000    0.0000    0.0000    4.7508   -0.0000;...
       0.0000    0.0000   -0.0000   28.6082    0.0000   -0.0000    0.0000    5.7271;...
       4.2030   -0.0000    0.0000    0.0000    7.0216   -0.0000    0.0000    0.0000;...
      -0.0000    5.7271    0.0000   -0.0000   -0.0000   11.1316    0.0000   -0.0000;...
       0.0000   -0.0000    4.7508    0.0000    0.0000    0.0000   14.9554   -0.0000;...
      -0.0000    0.0000   -0.0000    5.7271    0.0000   -0.0000   -0.0000   11.1316];

P2 = [5.2642   -2.2417;... 
     -2.2417    7.5398];

P3 = [5.5041   -0.0000   -0.0000   -0.0000   -2.3546    0.0000    0.0000    0.0000;...
     -0.0000    5.5041   -0.0000   -0.0000    0.0000   -2.3546    0.0000    0.0000;...
     -0.0000   -0.0000    5.5041   -0.0000    0.0000    0.0000   -2.3546    0.0000;...
     -0.0000   -0.0000   -0.0000    5.5041    0.0000    0.0000    0.0000   -2.3546;...
     -2.3546    0.0000    0.0000    0.0000    7.9037   -0.0000   -0.0000   -0.0000;...
      0.0000   -2.3546    0.0000    0.0000   -0.0000    7.9037   -0.0000   -0.0000;...
      0.0000    0.0000   -2.3546    0.0000   -0.0000   -0.0000    7.9037   -0.0000;...
      0.0000    0.0000    0.0000   -2.3546   -0.0000   -0.0000   -0.0000    7.9037];
  
v=0:0.01:T2;
lambda_min = zeros(size(v));
lambda_max = zeros(size(v));
e2hv = zeros(size(v));
a2 = zeros(size(v));
for i=1:max(size(v))
tt=v(i);
Matrix=expm(Af2'*tt)*P1*expm(Af2*tt);
a_2(i)=max(eig(Matrix));
Matrix=Ag2'*expm(Af2'*tt)*P1*expm(Af2*tt)*Ag2 - P1;
lambda_min(i)=min(eig(Matrix));
lambda_max(i)=max(eig(Matrix));
e2hv(i) = exp(2*H*tt);
end

a2 = max([max(e2hv) max(lambda_max) max(eig(P1)) max(eig(P2))]);

k1 = 2*max(a_2)*norm(Bf2);

k2_rge = [1,max(lambda_max)];

beta1 = norm(P2*Af3 + Af3'*P2);
beta2 = norm(P3*Af4 + Af4'*P3);

epsilon1 = min(roots([1 (beta2 - k1) k1]));
epsilon2 = max(roots([1 (beta2 - k1) k1]));

a2bar = max([1 max(eig(P1))]);
kbar1 = max([k1/(2*epsilon1) ((k1*epsilon1)/2 - beta2)]);
kbar2 = max([k1/(2*epsilon2) ((k1*epsilon2)/2 - beta2)]);
kbar = max([kbar1 kbar2]);

a1_til = a2bar/2;
a2_til = a2bar/2;

e = [x(:,1)  x(:,2)  x(:,3)  x(:,4)  x(:,5)]';
e_bar = zeros(size(e));
e_a = [x(:,26) x(:,27) x(:,28) x(:,29) x(:,30)]' - [x(:,31) x(:,32) x(:,33) x(:,34) x(:,35)]';
e_t = [x(:,21) x(:,22) x(:,23) x(:,24) x(:,25)]' - [x(:,16) x(:,17) x(:,18) x(:,19) x(:,20)]';
ea_bar = zeros(size(e_a));
et_bar = zeros(size(e_a));
eta = [x(:,11) x(:,12) x(:,13) x(:,14) x(:,15)]';
eta_bar = zeros(size(eta));

k2 = max(k2_rge);
%z_bar = zeros(size(z));

% while norm(exp((T2*kbar))*((1 - k2/a2))) > 1
%     k2 = k2_rge(1) + (k2_rge(2)-k2_rge(1))*rand(1);
% end

%\eps^2 + (\beta2-\kappa1) \eps + \kappa1 = 0

for i = 1:length(e_a)
    ea_bar(:,i) = (T^-1)*e_a(:,i);
    et_bar(:,i) = (T^-1)*e_t(:,i);
    eta_bar(:,i) = (T^-1)*eta(:,i);
    e_bar(:,i) = (T^-1)*e(:,i);
end

z_bar_1 = [e_bar(1,:); eta_bar(1,:);];
w_bar_1 = [ea_bar(1,:); et_bar(1,:);];
z_bar_2 = [e_bar([2,3,4,5],:); eta_bar([2,3,4,5],:);];
w_bar_2 = [ea_bar([2,3,4,5],:); et_bar([2,3,4,5],:);];

z_bar = [z_bar_1; z_bar_2];
w_bar = [w_bar_1; w_bar_2];

V_bound = zeros(length(e_a),1);
V_0 = exp(2*H*x(1,41))*norm(eta_bar(1,1))^2 + z_bar_2(:,1)'*expm(Af2'*x(1,41))*P1*expm(Af2*x(1,41))*z_bar_2(:,1) + w_bar_1(:,1)'*P2*w_bar_1(:,1) + w_bar_2(:,1)'*P3*w_bar_2(:,1);

%%

for i = 1:length(j)
    if i > 1
        if j(i) - j(i-1) > 0
            w_k(:,j(i)) = w_bar(:,i);
        end
    end
end


%%
epsilon = 1000;
V_bound(1) = V_0;
V(1) = V_0;
sum = 0;
for i = 1:(length(j)-1)
    if i > 1
        if j(i) - j(i-1) > 0
            jsum = j(i);
            sum = 0;
            for k=1:1:jsum
                sum = sum + exp((T2*kbar))*(norm(w_bar(:,i))^2);
            end
        end
    end
    %disp('exp')
    %exp((T2*(j(i)+1)*kappa)/2*epsilon)
    %disp('1-sigma')
    %(1 - sigma/a1_til)^(j(i))
    %disp('sigma*sum')
    %sigma*sum
    V_bound(i+1) = exp((T2*(j(i)+1)*kbar))*((1 - k2/a2)^(j(i)))*V_0 + k2*sum;
    V(i+1) = exp(2*H*x(i,41))*norm(eta_bar(1,i))^2 + z_bar_2(:,i)'*expm(Af2'*x(i,41))*P1*expm(Af2*x(i,41))*z_bar_2(:,i) + w_bar_1(:,i)'*P2*w_bar_1(:,i) + w_bar_2(:,i)'*P3*w_bar_2(:,i);
end

%%

figure(1) % position
clf
plot(t,x(:,1)-x(:,2));
hold on
plot(t,x(:,1)-x(:,3));
hold on
plot(t,x(:,1)-x(:,4));
hold on
plot(t,x(:,1)-x(:,5));
hold on
plot(t,x(:,2)-x(:,3));
hold on
plot(t,x(:,2)-x(:,4));
hold on
plot(t,x(:,2)-x(:,5));
hold on
plot(t,x(:,3)-x(:,4));
hold on
plot(t,x(:,3)-x(:,5));
hold on
plot(t,x(:,4)-x(:,5));
grid on
ylabel('$e_i - e_k$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)

%%

modificatorF{1} = 'b';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;
modificatorJ{1} = '--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '.';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'b';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'b';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 5;

modificatorV{1} = 'r';
modificatorV{2} = 'LineWidth';
modificatorV{3} = 1;
modificatorM{1} = '--';
modificatorM{2} = 'LineWidth';
modificatorM{3} = 1;
modificatorM{4} = 'Marker';
modificatorM{5} = '.';
modificatorM{6} = 'MarkerEdgeColor';
modificatorM{7} = 'r';
modificatorM{8} = 'MarkerFaceColor';
modificatorM{9} = 'r';
modificatorM{10} = 'MarkerSize';
modificatorM{11} = 5;

e_diff = [x(:,1)-x(:,2) x(:,1)-x(:,3) x(:,1)-x(:,4) x(:,1)-x(:,5) x(:,2)-x(:,3) x(:,2)-x(:,4) x(:,2)-x(:,5) x(:,3)-x(:,4) x(:,3)-x(:,5) x(:,4)-x(:,5)];

ymin1 = min(min(e_diff)) - round(max(mean(e_diff)),1);
ymax1 = max(max(e_diff)) + round(max(mean(e_diff)),1);

ymin2 = min(min([x(:,26)-x(:,31) x(:,27)-x(:,32) x(:,28)-x(:,33) x(:,29)-x(:,34) x(:,30)-x(:,35)]))...
        - round(max(mean([x(:,26)-x(:,31) x(:,27)-x(:,32) x(:,28)-x(:,33) x(:,29)-x(:,34) x(:,30)-x(:,35)])),2);
ymax2 = max(max([x(:,26)-x(:,31) x(:,27)-x(:,32) x(:,28)-x(:,33) x(:,29)-x(:,34) x(:,30)-x(:,35)]))...
        + round(max(mean([x(:,26)-x(:,31) x(:,27)-x(:,32) x(:,28)-x(:,33) x(:,29)-x(:,34) x(:,30)-x(:,35)])),2);

figure(4)
clf
subplot(4,1,1), plot(t,x(:,1)-x(:,2));
hold on
subplot(4,1,1), plot(t,x(:,1)-x(:,3));
hold on
subplot(4,1,1), plot(t,x(:,1)-x(:,4));
hold on
subplot(4,1,1), plot(t,x(:,1)-x(:,5));
hold on
subplot(4,1,1), plot(t,x(:,2)-x(:,3));
hold on
subplot(4,1,1), plot(t,x(:,2)-x(:,4));
hold on
subplot(4,1,1), plot(t,x(:,2)-x(:,5));
hold on
subplot(4,1,1), plot(t,x(:,3)-x(:,4));
hold on
subplot(4,1,1), plot(t,x(:,3)-x(:,5));
hold on
subplot(4,1,1), plot(t,x(:,4)-x(:,5));
grid on
ylabel('$e_i - e_k$','Interpreter','latex','FontSize',30)
grid on
set(gca,'FontSize',20)
axis([0 t(end) ymin1 ymax1])
h = findobj(gca,'Type','line');

subplot(4,1,2), plot(t,x(:,26)-x(:,31));
hold on
subplot(4,1,2), plot(t,x(:,27)-x(:,32));
hold on
subplot(4,1,2), plot(t,x(:,28)-x(:,33));
hold on
subplot(4,1,2), plot(t,x(:,29)-x(:,34));
hold on
subplot(4,1,2), plot(t,x(:,30)-x(:,35));
grid on
set(gca,'FontSize',20)
axis([0 t(end) ymin2 ymax2])
h = findobj(gca,'Type','line');
i = legend([h(1) h(2) h(3) h(4) h(5)],'$$\varepsilon_{a_1}$$','$$\varepsilon_{a_2}$$','$$\varepsilon_{a_3}$$','$$\varepsilon_{a_4}$$','$$\varepsilon_{a_5}$$');
set(i,'Interpreter','latex','FontSize',15)
ylabel('$\varepsilon_{a_i}$','Interpreter','latex','FontSize',30)
y2h = get(gca,'ylabel');
gy2 = get(y2h);                                                         % Object Information
y2p = get(y2h, 'Position');
set(y2h, 'Rotation',0, 'Position',y2p, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
subplot(4,1,3), plotHarc(t,j,x(:,41));
grid on
axis([0 t(end) 0 T2])
set(gca,'FontSize',20)
ylabel('$\tau$','Interpreter','latex','FontSize',40)
y3h = get(gca,'ylabel');
gy3 = get(y3h);                                                         % Object Information
y3p = get(y3h, 'Position');
set(y3h, 'Rotation',0, 'Position',y3p, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(y3h, 'Rotation',0, 'Position',y3p-[0 0 0], 'VerticalAlignment','middle', 'HorizontalAlignment','right')
subplot(4,1,4), plotHarc(t,j,V_bound(:),[],modificatorV,modificatorM);
hold on
subplot(4,1,4), plotHarc(t,j,V(:),[],modificatorF,modificatorJ);
grid on
h = findobj(gca,'Type','line');
i = legend([h(10) h(21)],'$$V$$','$$V_b$$');
set(i,'Interpreter','latex','FontSize',35)
xlabel('$t,j$','Interpreter','latex','FontSize',40)
set(gca,'FontSize',20)

figure(6) % position
clf
plot(t,x(:,42));
hold on
plot(t,x(:,43));
hold on
plot(t,x(:,44));
hold on
plot(t,x(:,45));
hold on
plot(t,x(:,46));
grid on
ylabel('$\tilde{\tau}_i$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)