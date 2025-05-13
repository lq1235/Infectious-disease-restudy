

clc
clear
a=xlsread("data.xlsx")

plot(a(:,2),log(a(:,3)),'r-o','LineWidth',1.5,'Markersize',10);

% axis([0 15,3 16])
xlabel('\fontsize{27}\DeltaU')
ylabel('\fontsize{27}logT_{MFPT}')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% xlim([0.683725974535662,6])
% ylim([7.525260623202302,13.211388621425415])
set(gca,'xtick',0:1:13)

