clc
clear
% a=load('MFPT.txt');
a=xlsread('MFPT.xlsx');
%MFPT
h=figure(1)
plot(a(:,1),a(:,2),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
hold on

xlabel("\fontsize{27} \delta");
ylabel("\fontsize{27} MFPT");

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0.001 0.0025,-0.00005 0.001000001])
xlim([7 11])
ylim([-57309.7662515773,1042690.233748423])
set(gca,'ytick',0:200000:1000000)
ax = gca();
ax.YRuler.Exponent = 5;

h=figure(2)
plot(a(:,1),1./a(:,2),'sk-','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),1./a(:,2),'k.','LineWidth',1,'markersize',10)
hold on

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
axis([7 11,-0.00007030817779763656,0.0005457416503489927])
set(gca,'ytick',0:0.0001:0.0040)
% xlim([0.001 0.00175])
% ylim([-0.00021343 0.00054994])
xlabel("\fontsize{27} \delta");
ylabel("\fontsize{27}\nu");
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
% set(gca,'xtick',0:0.00025:0.003)


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
