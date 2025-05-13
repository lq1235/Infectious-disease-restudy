clc;
clear;
seed=1;rng(seed); 
dim=7;
a1=0.5:0.5:12;
threshold = 0.4;

h=0.1;
dc=0.0000017;
steps=10000;               
iters=10000; 
n_traj=1000;
n_varpara=length(a1);

% Sets the number of times the data is simulated under the parameters stored
k1=11;  i1=10;
% Sets the sequence to be stored and calculates the two sequences used for cross-correlation
series1=2; series2=3;


for k=1:n_varpara             
    k
    i = 1; 
func=@(x)force_flower(x,a1(k)); 
x0=load('fixpoint1.txt');    
x0=x0(k,2:end);          
x0=x0';



%设置方差和的初始值为0
variance1=zeros(7,1);

%设置振幅初始值为0
amplitude1=zeros(7,1);






while i<=n_traj   
[t,x]=sode_Heue(func,h,x0,iters,dc);

if any(x(1,:)<threshold)                                        
   continue;                                                  
end



        % 计算矩阵每一行的方差
variance = var(x,1, 2);

% 计算数据的方差
% variance = var(x(1,:), 1); % 使用 "n" 的估计值
variance1=variance+variance1;

amplitude=max(x,[],2)-min(x,[],2);

amplitude1=amplitude1+amplitude;









if k == k1 && i == i1      % Store the number of sequences of track % store track % Here can add an if statement to judge, and then specify the calculated parameters of the track save, otherwise the memory consumption is too large
x1=x';
dlmwrite('series.txt', [x1(:,[series1,series2]),], 'delimiter', '\t');% stores the trajectories of x1 and x2
end

% Calculated autocorrelation
traj_1=x(series1,:);
[acf, lags]=autocorr(traj_1,'NumLags',steps-1);
t0=lags(1:1:end)*h;
index(k,i)=min(find(acf<0));
tau_ik(k,i)=trapz(t(1:index(k,i)-1),acf(1:index(k,i)-1)); 
tau=mean(tau_ik,2); 

traj_2=x(series2,:);
[acf1, lags]=autocorr(traj_2,'NumLags',steps-1);
t0=lags(1:1:end)*h;
index1(k,i)=min(find(acf1<0));
tau_ik1(k,i)=trapz(t(1:index1(k,i)-1),acf1(1:index1(k,i)-1)); 
tau1=mean(tau_ik1,2); 

% Autocorrelation storage
if k == k1 && i == i1 
dlmwrite('autocorrelation.txt', [t0',  acf'], 'delimiter', '\t');
dlmwrite('autocorrelation1.txt', [t0',  acf1'], 'delimiter', '\t');
end


% calculates cross correlatio
traj1=x([series1,series2],1:end);  
[xcf,lags1]=crosscorr(traj1(1,:),traj1(2,:),'NumLags',steps-1); 
CXY=xcf(steps:1:2*steps-1);
CYX=xcf(steps:-1:1);
deltaC=CXY-CYX;
t1=lags1(steps:2*steps-1)*h;

n=5000;% averages t approaching 0
DeltaC(k,i)=sqrt(trapz( t1(1:n) ,deltaC(1:n).^2 )/(h*n));
DeltaC1=mean(DeltaC,2);
i=i+1;


if k == k1 && i == i1 
dlmwrite('crosscorrelation.txt', [t1',  CXY', CYX'], 'delimiter', '\t');
dlmwrite('crosscorrelation-difference.txt', [t1',  deltaC'], 'delimiter', '\t');
end

end



    data(k,:)=variance1/n_traj;
% amplitude1/n_traj

    data1(k,:)=amplitude1/n_traj;





end 

dlmwrite('tau-tau1-DeltaC.txt', [a1', tau, tau1, DeltaC1], 'delimiter', '\t');

filename='CSD_DCC1.mat';
save(filename);

% critical slowing down
h=figure(1)
plot(a1,tau,'k-o','LineWidth',1,'markersize',10)
hold on
plot(a1,tau,'k.','LineWidth',1,'markersize',10)
hold on
title('series1')
xlabel('\delta','FontSize',27);
ylabel('\tau','FontSize',27);



h1=figure(2)
plot(a1,tau1,'k-o','LineWidth',1,'markersize',10)
hold on
plot(a1,tau1,'k.','LineWidth',1,'markersize',10)
hold on
title('series2')
xlabel('\delta','FontSize',27);
ylabel('\tau','FontSize',27);


%time irreversibility
h2=figure(3)
plot(a1,DeltaC1,'k-o','LineWidth',1,'markersize',10)
hold on
plot(a1,DeltaC1,'k.','LineWidth',1,'markersize',10)
hold on
xlabel('\delta','FontSize',27);
ylabel('\DeltaCC','FontSize',27);


figure(4)
plot(traj1(1,:),traj1(2,:))




h5=figure(5)
plot(a1,data(:,3),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a1,data(:,3),'k.','LineWidth',1,'markersize',10)
hold on

xlabel('\delta','FontSize',27);
ylabel('variance','FontSize',27);

set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ax = gca();
ax.YRuler.Exponent = -4;
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
ax = gca();
ax.XRuler.Exponent = 0;
xlim([0 15])
% set(gca,'xtick',5:5:20)
% set(gca,'ytick',0.3:0.05:5.5)
ylim([-0.000124897959176,0.001855102040824])
print(h5, '-r600', '-dpdf', ['variance', num2str(1),'.pdf']);


h6=figure(6)
plot(a1,data1(:,3),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a1,data1(:,3),'k.','LineWidth',1,'markersize',10)
hold on

xlabel('\delta','FontSize',27);
ylabel('amplitude','FontSize',27);
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ax = gca();
ax.YRuler.Exponent = -2;
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
ax = gca();
ax.XRuler.Exponent = 0;
xlim([0 15])
% set(gca,'xtick',5:5:20)
% set(gca,'ytick',0.3:0.05:5.5)
ylim([0.018726989195729,0.194726989195729])
print(h6, '-r600', '-dpdf', ['amplitude', num2str(1),'.pdf']);


   
   











function f=force_flower(x,a1)
f=zeros(length(x),1);   
%N = 10^4;
%K = 10^5;
w = 0.04;
r = 0.002;
q = 0.0064;%q = 0.0064;
p=0.0018;

f(1) = q*x(2) - p*(a1)*x(4);
f(2) = r*(1-x(1)-x(2)) - q*x(2);
f(3) = q*x(6) + w*x(1)*x(4)/(x(1)+x(2)) - 2*p*(a1)*(x(3)*x(4))/x(1);
f(4) = 2*p*(a1)*(x(3)*x(4))/x(1) + q*x(7) - r*x(4) - w*x(4) - p*(x(4)+(a1)*x(4)*x(4)/x(1));
f(5) = p*(x(4)+(a1)*x(4)*x(4)/x(1)) - 2*r*x(5);
f(6) = r*x(4) + w*x(2)*x(4)/(x(1)+x(2)) + 2*q*(1-x(3)-x(4)-x(5)-x(6)-x(7)) - q*x(6) - p*(a1)*x(4)*x(6)/x(1) + w*x(1)*x(7)/(x(1)+x(2));
f(7) = 2*r*x(5) + p*(a1)*x(4)*x(6)/x(1) - q*x(7) - r*x(7) - w*x(7);
% N = 10^4;
% K = 10^5;
% %w = 0.04;
% r = 0.002;
% q = 0.0064;%q = 0.0064;
% p=0.0018;
% 
% 
% f(1) = q*x(2) - p*(K/N)*x(4);
% f(2) = r*(1-x(1)-x(2)) - q*x(2);
% f(3) = q*x(6) + a1*x(1)*x(4)/(x(1)+x(2)) - 2*p*(K/N)*(x(3)*x(4))/x(1);
% f(4) = 2*p*(K/N)*(x(3)*x(4))/x(1) + q*x(7) - r*x(4) - a1*x(4) - p*(x(4)+(K/N)*x(4)*x(4)/x(1));
% f(5) = p*(x(4)+(K/N)*x(4)*x(4)/x(1)) - 2*r*x(5);
% f(6) = r*x(4) + a1*x(2)*x(4)/(x(1)+x(2)) + 2*q*(1-x(3)-x(4)-x(5)-x(6)-x(7)) - q*x(6) - p*(K/N)*x(4)*x(6)/x(1) + a1*x(1)*x(7)/(x(1)+x(2));
% f(7) = 2*r*x(5) + p*(K/N)*x(4)*x(6)/x(1) - q*x(7) - r*x(7) - a1*x(7);
end

function [t,x] = sode_Heue(func,h,x0,steps,dc)
dim = length(x0);
t = 0:h:(steps-1)*h;
h_sqrt=sqrt(h);
gx = sqrt(2*dc);
  x(:,1) = x0;
  Xm=1.0;% Boundary
%  Xm1=ones(dim,1)*Xm;
for j = 1:steps-1 
    noise = h_sqrt*randn(dim,1);
    gxn = gx*noise;
    aux = h*feval(func,x(:,j))+gxn;
    fxh =feval(func,x(:,j)+aux);
    xn = x(:,j)+ 0.5*(aux+h*fxh+gxn);

    if any(xn<0)
        x(:,j+1) = abs(xn);  % minimum boundary takes absolute value, which can only correspond to C++ when the minimum boundary is 0
% maximum boundary processing               
    elseif any(xn > Xm)
            indices = find(xn > Xm);
            xn(indices) = 2 * Xm - xn(indices);
            x(:,j+1) = xn;
    else
       x(:,j+1) = xn;
    end
end

end
