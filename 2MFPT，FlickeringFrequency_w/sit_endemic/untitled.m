clc
clear
a=load('MFPT_h=0.002_r=1.00.txt');
%MFPT
h=figure(1)
plot(a(:,1),a(:,2),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
hold on

xlabel("\fontsize{27} w");
ylabel("\fontsize{27} MFPT");

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0.001 0.0025,-0.00005 0.001000001])
% xlim([0.0012 0.0018])
% ylim([-543080 1343300])
% set(gca,'xtick',0.0012:0.0002:0.0018)


h=figure(2)
plot(a(:,1),1./a(:,2),'sk-','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),1./a(:,2),'k.','LineWidth',1,'markersize',10)
hold on

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0 0.00395,-0.00005 0.001000001])
% set(gca,'xtick',0:0.0005:0.0040)
xlim([0.03 0.063])
ylim([-0.00048257167225,0.00484142832775])
xlabel("\fontsize{27} w");
ylabel("\fontsize{27}\nu");
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'ytick',0:0.001:0.01)


h=figure(3)
plot(a(:,1),log(a(:,2)),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),log(a(:,2)),'k.','LineWidth',1,'markersize',10)
hold on

h=figure(4)
plot(a(:,1),log(1./a(:,2)),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),log(1./a(:,2)),'k.','LineWidth',1,'markersize',10)
hold on
