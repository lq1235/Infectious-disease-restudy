clc;
clear;
%Xm=1.0;
% index1=3;
% index2=4;
% dif =[0.05;0.05];
ms = 200;
xmin = [0;0];
xmax = [0.27;1];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

a1 =0.001:0.002:0.101;
num = length(a1);
file_path = 'pp23_w=%0.5f.txt';
for j = 1:num

h=figure(j);

sample=sprintf(file_path,a1(j));
% sample = 'pp23_w=0.03500.txt';
px = load(sample);
p = reshape(px(:,3),ms,ms);
sum(sum(p))

FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));

Pi = p/z;      

%第一种处理方法
PP = eq(Pi,0)+Pi;     
P_eps=min(min(PP));  
P = P_eps*eq(Pi,0)+Pi;
% %第二种方法
% eps=1.1e-0;
% P=Pi+eps;


U = -log(P);



% U_max =16.5;%U_max =8.6;
% U=U.*(U<U_max)+U_max.*(U>U_max);
% U(U==U_max)=NaN;
surf(X,Y,U)
% mesh(X,Y,U);
% pcolor(X,Y,U);
shading interp;
colormap([jet(256)]);

xlabel('\fontsize{25} P_{R}');
ylabel('\fontsize{25} P_{SS}');



% colorbar;
% caxis([-2,U_max]);%caxis([-2,8.6]);
% caxis([min(min(Ux)),U_max]);
view([-45.497542818583412,14.180563854604134]);
% al=min(min(U(88:155,1:52)))
% al1=min(min(U(1:52,88:155)))

zlabel('\fontsize{25} U')
axis([0 0.27,0 1])
% zlim([-5 20])
% set(gca,'xtick',0:1:4)
% set(gca,'ytick',0:1:4)
% set(gca,'ztick',0:5:15)
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])



filename = ['Landscape_',  'w=', num2str(a1(j), '%.5f'), '.pdf'];
print(h, '-r600', '-dpdf', filename);

% colorbar;
% caxis([-2,U_max]);%caxis([-2,8.6]);
% caxis([min(min(Ux)),U_max]);
% view([-44.54,13.855]);
% al=min(min(U(88:155,1:52)))
% al1=min(min(U(1:52,88:155)))




% zlabel('\fontsize{25} U')
axis([0 0.27 ,0 1])
% zlim([-5 20])
% set(gca,'xtick',0:1:4)
% set(gca,'ytick',0:1:4)
% % set(gca,'ztick',0:5:15)
% set(gca,'LineWidth',1.2,'Fontsize',25)
% set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])


% print(h, '-r600', '-dpdf', ['C15_', num2str(j),'.pdf']);
% print(h, '-r600', '-depsc', ['C15_', num2str(j),'.eps']);

end