clc
clear
a=load('data.txt');

% a(:,2)-a(:,3)

a(:,2)-a(:,4)
%endemic
figure(1)
plot(a(:,1),a(:,2)-a(:,3),'s-r','LineWidth',1,'Markersize',10)
hold on
plot(a(:,1),a(:,2)-a(:,3),'r.','LineWidth',1,'Markersize',10)
hold on
% axis([-0.2 3.2,-0.1 15])

xlabel('\fontsize{27} \delta')
ylabel('\fontsize{27}\DeltaU')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% set(gca,'xtick',0:0.5:10)
set(gca,'xtick',0:0.5:13)
xlim([7.5 10])
ylim([-0.437497300766564,8.374007125813787])

%disease-free
figure(2)
plot(a(:,1),a(:,2)-a(:,4),'s-k','LineWidth',1,'Markersize',10)
hold on
plot(a(:,1),a(:,2)-a(:,4),'k.','LineWidth',1,'Markersize',10)
hold on
% axis([-0.2 3.2,-0.1 15])

xlabel('\fontsize{27} \delta')
ylabel('\fontsize{27}\DeltaU')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])

set(gca,'xtick',0:0.5:13)
xlim([7.5 10])
ylim([-1.010458569618501,12.299541430381515])

