clc;
clear;
close all
% Xm=3.0;
index1=2;
index2=3;
a1 = 0.003  :0.002:0.097;
% a_1="C15";
dif =[0.0000017;0.0000017];
ms = 200;
xmin = [0;0];
xmax = [0.27;1];
num = length(a1);
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);
file_path = 'pp%d%d_w=%0.5f.txt';
for j = 1:num
%    sample=sprintf(file_path,index1,index2,a_1,a1(j),Xm);   
     sample=sprintf(file_path,index1,index2,a1(j)); 
%method one
px = load(sample);
p = reshape(px(:,3),ms,ms);
FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));
Pi = p/z;
PP = eq(Pi,0)+Pi;   
P_eps=min(min(PP));  
P = P_eps*eq(Pi,0)+Pi;
U = -log(P);
dx = x(2)-x(1);
dy = y(2)-y(1);
[GPx,GPy] = gradient(P,dx,dy);
Jx = FPx.*P - dif(1)*GPx ;
Jy = FPy.*P - dif(2)*GPy ;
JJP=(Jx.^2)./(dif(1)*P)+(Jy.^2)./(dif(2)*P);
EPR(j)=trapz(y,trapz(x,JJP));
FJ=(FPx.*Jx+FPy.*Jy)/dif(1);
HDR(j)=trapz(y,trapz(x,FJ));
JJ=sqrt(Jx.^2+Jy.^2);
Flux(j) = trapz(y,trapz(x,JJ)); 

%method two
% P_eps =1.0e-4;
% px = load(sample);
% p= reshape(px(:,3),ms,ms);
% FPx = reshape(px(:,4),ms,ms);
% FPy= reshape(px(:,5),ms,ms);
% z = trapz(y,trapz(x,p));
% P = p/z;
% U = -log(P+P_eps);
% dx = x(2)-x(1);
% dy = y(2)-y(1);
% [GPx,GPy] = gradient(P,dx,dy);
% Jx = FPx.*P - dif(1)*GPx ;
% Jy = FPy.*P - dif(2)*GPy ;
% JJP=(Jx.^2)./(dif(1)*(P+P_eps))+(Jy.^2)./(dif(2)*(P+P_eps));
% EPR(j)=trapz(y,trapz(x,JJP));
% FJ=(FPx.*Jx+FPy.*Jy)/dif(1);
% HDR(j)=trapz(y,trapz(x,FJ));
% JJ=sqrt(Jx.^2+Jy.^2);
% Flux(j) = trapz(y,trapz(x,JJ));
end

h=figure(1)
plot(a1,EPR,'k-o','LineWidth',1,'markersize',10);
hold on
plot(a1,EPR,'k.','LineWidth',1,'markersize',10);
hold on

xlabel('\fontsize{27} w');
ylabel('\fontsize{27} EPR');


set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% axis([0 0.00395,-0.0035 0.04])


% axis([0 0.00395,0 0.00025])
set(gca,'xtick',0:0.02:0.1)
ylim([-0.007179399671724,0.065420600328276])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
ax = gca();
ax.YRuler.Exponent = -2;
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度




% plot([0.7 0.7],[4.5*10^(-4) 6.7*10^(-4)],'r--','LineWidth',1)
% hold on
% plot([2.2 2.2],[4.5*10^(-4) 6.7*10^(-4)],'b--','LineWidth',1)

% 
% xlabel('\fontsize{27} C_{15}');
% ylabel('\fontsize{27} EPR');
% 
% set(gca,'xtick',0:0.0005:0.0040)
% set(gca,'LineWidth',1.2,'Fontsize',27)
% set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])

% print(h, '-r600', '-depsc', ['EPR_C15', num2str(1),'.eps']);
% print(h, '-r600', '-dpdf', ['EPR_C15', num2str(1),'.pdf']);

h=figure(2)
plot(a1,Flux,'ko-','LineWidth',1,'markersize',10);
hold on
plot(a1,Flux,'k.','LineWidth',1,'markersize',10);
hold on

xlabel('\fontsize{27} w');
ylabel('\fontsize{27} Flux');

set(gca,'xtick',0:0.02:0.1)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ylim([-0.00003181546064835954,0.0002866998011184508])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
ax = gca();
ax.YRuler.Exponent = -4;
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度



figure(3)

plot(a1,HDR,'ko-','LineWidth',1,'markersize',10);
hold on
plot(a1,HDR,'k.','LineWidth',1,'markersize',10);
hold on

xlabel('\fontsize{27} w');
ylabel('\fontsize{27} HDR');

set(gca,'xtick',0:0.02:0.1)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ylim([-0.008560944456821,0.064039055543179])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
ax = gca();
ax.YRuler.Exponent = -2;
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度


