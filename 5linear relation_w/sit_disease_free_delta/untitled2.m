

clc
clear
a=xlsread('data.xlsx');

plot(a(:,2),log(a(:,3)),'k-o','LineWidth',1.5,'Markersize',10);

% axis([0 15,3 16])
xlabel('\fontsize{27}\DeltaU')
ylabel('\fontsize{27}logT_{MFPT}')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlim([-0.040762803870582,13.136137196129425])
ylim([8.25078120955605,15.27846120955606])
set(gca,'xtick',0:2:13)

set(gca,'XTickLabelRotation',0);%46是字体的旋转