clc
clear
a=load('data.txt');



a(:,2)-a(:,3)


% a(:,2)-a(:,4)

%endemic
figure(1)
plot(a(:,1),a(:,2)-a(:,3),'s-r','LineWidth',1,'Markersize',10)
hold on
plot(a(:,1),a(:,2)-a(:,3),'r.','LineWidth',1,'Markersize',10)
hold on
% axis([-0.2 3.2,-0.1 15])

xlabel('\fontsize{27} w')
ylabel('\fontsize{27}\DeltaU')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
xlim([0.039553497576352,0.054553497576352])
ylim([-1.210443681629535,7.024040708370472])
% set(gca,'xtick',0:0.5:10)
set(gca,'xtick',0.021:0.004:0.056)
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
hold on




%disease-free
figure(2)
plot(a(:,1),a(:,2)-a(:,4),'s-k','LineWidth',1,'Markersize',10)
hold on
plot(a(:,1),a(:,2)-a(:,4),'k.','LineWidth',1,'Markersize',10)
hold on
% axis([-0.2 3.2,-0.1 15])

xlabel('\fontsize{27} w')
ylabel('\fontsize{27}\DeltaU')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])

set(gca,'xtick',0.021:0.004:0.056)
xlim([0.039553497576352,0.054553497576352])
ylim([-1.879492155247867,16.49640274475214])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
% ax = gca();
% ax.YRuler.Exponent = -2;




