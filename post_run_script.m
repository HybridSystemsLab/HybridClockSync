%%

e = [x(:,1)  x(:,2)  x(:,3)  x(:,4)  x(:,5)]';
e_bar = zeros(size(e));
e_a = [x(:,26) x(:,27) x(:,28) x(:,29) x(:,30)]' - [x(:,31) x(:,32) x(:,33) x(:,34) x(:,35)]';
e_t = [x(:,21) x(:,22) x(:,23) x(:,24) x(:,25)]' - [x(:,16) x(:,17) x(:,18) x(:,19) x(:,20)]';
ea_bar = zeros(size(e_a));
et_bar = zeros(size(e_a));
eta = [x(:,11) x(:,12) x(:,13) x(:,14) x(:,15)]';
eta_bar = zeros(size(eta));
sum_e = zeros(size(x(1,:)));

for i = 1:length(e_a)
    ea_bar(:,i) = (T^-1)*e_a(:,i);
    et_bar(:,i) = (T^-1)*e_t(:,i);
    eta_bar(:,i) = (T^-1)*eta(:,i);
    e_bar(:,i) = (T^-1)*e(:,i);
    sum_e(i) = ones(1,N)*L*e(:,i);
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

a1_w = min([min(eig(P2)) min(eig(P3))]);
a2_w = max([max(eig(P2)) max(eig(P3))]);
betaT = max([beta1 beta2]);
V_bound(1) = V_0;
V(1) = V_0;
sum_1 = 0;
for i = 1:(length(j)-1)
    if i > 1
        if j(i) - j(i-1) > 0
            jsum = j(i);
            sum_1 = 0;
            for k=1:1:jsum
                sum_1 = sum_1 + exp((T2*kbar))*(norm(w_bar(:,i))^2);
            end
        end
    end
    V_bound(i+1) = exp(T2*kbar/a2)*((exp(T2*kbar/a2)*(1 - k2/a2))^(j(i)))*V_0 + k2*exp((T2*kbar/a2))*norm(sqrt((a2_w/a1_w)*exp(-(betaT/2*a2_w)*(t(i))))*norm(w_bar(:,1)))^2;
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

figure(2)
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

figure(3) % position
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

%%

figure(4)
clf
subplot(3,1,1), plot(t,x(:,1)-x(:,2),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,1)-x(:,3),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,1)-x(:,4),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,1)-x(:,5),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,2)-x(:,3),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,2)-x(:,4),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,2)-x(:,5),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,3)-x(:,4),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,3)-x(:,5),'LineWidth',1);
hold on
subplot(3,1,1), plot(t,x(:,4)-x(:,5),'LineWidth',1);
grid on
ylabel('$e_i - e_k$','Interpreter','latex','FontSize',40)
grid on
set(gca,'FontSize',20)
axis([0 t(end) ymin1 ymax1])
h = findobj(gca,'Type','line');
set(gca,'linewidth',1)
subplot(3,1,2), plot(t,x(:,26)-x(:,31),'LineWidth',1);
hold on
subplot(3,1,2), plot(t,x(:,27)-x(:,32),'LineWidth',1);
hold on
subplot(3,1,2), plot(t,x(:,28)-x(:,33),'LineWidth',1);
hold on
subplot(3,1,2), plot(t,x(:,29)-x(:,34),'LineWidth',1);
hold on
subplot(3,1,2), plot(t,x(:,30)-x(:,35),'LineWidth',1);
grid on
set(gca,'linewidth',1)
set(gca,'FontSize',20)
axis([0 t(end) ymin2 ymax2])
h = findobj(gca,'Type','line');
i = legend([h(1) h(2) h(3) h(4) h(5)],'$$\varepsilon_{a_1}$$','$$\varepsilon_{a_2}$$','$$\varepsilon_{a_3}$$','$$\varepsilon_{a_4}$$','$$\varepsilon_{a_5}$$');
set(i,'Interpreter','latex','FontSize',15)
ylabel('$\varepsilon_{a_i}$','Interpreter','latex','FontSize',40)
y2h = get(gca,'ylabel');
gy2 = get(y2h);                                                         % Object Information
y2p = get(y2h, 'Position');
set(y2h, 'Rotation',0, 'Position',y2p, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
subplot(3,1,3), plotHarc(t,j,x(:,41));
grid on
axis([0 t(end) 0 T2])
set(gca,'FontSize',20)
ylabel('$\tau$','Interpreter','latex','FontSize',40)
xlabel('$t$','Interpreter','latex','FontSize',40)
set(gca,'FontSize',20)
set(gca,'linewidth',1)

%%

figure(5)
clf
plotHarc(t,j,V_bound(:),[],modificatorV,modificatorM);
hold on
plotHarc(t,j,V(:),[],modificatorF,modificatorJ);
grid on
h = findobj(gca,'Type','line');
i = legend([h(10) h(21)],'$$V$$','$$V_b$$');
set(i,'Interpreter','latex','FontSize',35)
xlabel('$t,j$','Interpreter','latex','FontSize',40)
set(gca,'FontSize',20)