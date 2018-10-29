%% Plot enclosed cylinder
clear; clc

% Sample values
len = 5;     % height
ra = 1;   % radius

% Create constant vectors
tht = linspace(0,2*pi,200); z = linspace(0,len,40);

% Create cylinder
xa = repmat(ra*cos(tht),40,1); ya = repmat(ra*sin(tht),40,1);

za = repmat(z',1,200);
% To close the ends

X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0];
Z = [za; flipud(za); za(1,:)];



% tht1=01:1/50:2*pi-1/50;
tht1=linspace(0,2*pi,30);
z1=linspace(0,len,10);
xa1=repmat(ra*cos(tht1),10,1);ya1=repmat(ra*sin(tht1),10,1);
za1=repmat(z1',1,30);
X1=[xa1*0; flipud(xa1); (xa1(1,:))*0];Y1 = [ya1*0; flipud(ya1); (ya1(1,:))*0];
Z1 = [za1; flipud(za1); za1(1,:)];

%% 无扭振
n=0:0.005:30-0.05;
N=length(n);
w=-0.5;
y=exp(j*2*pi*w*n);
n1=n(1:10:end);
y1=y(1:10:end);
%% figure1
figure
h=plot3(n,zeros(1,N),zeros(1,N),'k');
hold on;
h=plot3(ones(2,1)*n1,[zeros(1,length(n1));imag(y1)],[zeros(1,length(n1));real(y1)],'y');%'color',[0.1,0.1,0.1]
%   绘制螺旋面
    h=surf(ones(1,1)*n,[zeros(1,length(n));imag(y)],[zeros(1,length(n));real(y)],...
        'facecolor',[0.75 0.75 0.75],'edgecolor','none');%'color',[0.1,0.1,0.1]
    alpha(0.3)
    % hp- 绘制动点P的起点
    hp=plot3(0,imag(y(1)),real(y(1)),'y.','markersize',15);
    hp1=plot3(n(end),imag(y(end)),real(y(end)),'g.','markersize',15);
    hp2=quiver3(0,0,0,0,imag(y(1)),real(y(1)),'y','linewidth',1);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]
%     hp3=quiver3(n(end),0,0,0,imag(y(end)),real(y(end)),'g','linewidth',2);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]

%    h=plot3(n,imag(y),real(y),'r');
%   绘制螺旋线
    h=plot3(n,imag(y),real(y),'r');
%   绘制实部投影直线
%     h=plot3(ones(2,1)*n,1*ones(2,N),[zeros(1,N);real(y)],'color',[0.8,0.8,0.8]);
%   绘制实部投影曲线
    h=plot3(n,1*ones(1,N),1.1*real(y),'b','linewidth',0.5);
%     h=plot3(n,1.1*ones(1,N),real(y),'k');
%   绘制虚部投影直线
%     h=plot3(ones(2,1)*n,[zeros(1,N);imag(y)],-1*ones(2,N),'color',[0.8,0.8,0.8]);
%   绘制虚部投影曲线
    h=plot3(n,1.15*imag(y),-1*ones(size(n)),'k','linewidth',0.5);
%     h=plot3(n,imag(y),-1.1*ones(size(n)),'k');

% 
% 修饰figure
 ylabel('Imaginary Part');%,'Rotation',-6
 set(gca,'XTick',-5:6:30)
set(gca,'XTickLabel',{'Shaft ','0','0.2','0.4','0.6','0.8','1.0'})
zlabel('Real Part');
xlabel('time')
% view(27.5,22) %设置视角
mesh(Z-5,X*0.99,Y*0.99,'facecolor',[0.9 0.9 0.9],'edgecolor','none')
% hold on
surf(Z1-5,X1,Y1,'facecolor','none','edgecolor',[0.3 0.3 0.3]);
camlight;
lighting gouraud;
% light
% lighting phong
material metal
% axis equal;
% view(60,15) 

ylim([-1.5 1.5])
zlim([-1.5 1.5])
hold off
grid on

%% 扭振下的定点轨迹
% Plot enclosed cylinder
clear; clc

% Sample values
len = 5;     % height
ra = 1;   % radius
% Create constant vectors
tht = linspace(0,2*pi,200); z = linspace(0,len,40);
% Create cylinder
xa = repmat(ra*cos(tht),40,1); ya = repmat(ra*sin(tht),40,1);
za = repmat(z',1,200);
% To close the ends
X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0];
Z = [za; flipud(za); za(1,:)];

% tht1=01:1/50:2*pi-1/50;
tht1=linspace(0,2*pi,30);
z1=linspace(0,len,10);
xa1=repmat(ra*cos(tht1),10,1);ya1=repmat(ra*sin(tht1),10,1);
za1=repmat(z1',1,30);
X1=[xa1*0; flipud(xa1); (xa1(1,:))*0];Y1 = [ya1*0; flipud(ya1); (ya1(1,:))*0];
Z1 = [za1; flipud(za1); za1(1,:)];

%% 扭振
n=0:0.005:30-0.05;
N=length(n);
w=-0.5;
wt=1*sin(2*pi*0.1*n);

y=exp(j*2*pi*w*n+j*2.5*sin(2*pi*0.1*n));
n1=n(1:2:end);
y1=y(1:2:end);
%%
figure
h=plot3(n,zeros(1,N),zeros(1,N),'k');
hold on;
%   h=plot3(ones(2,1)*n1,[zeros(1,length(n1));imag(y1)],[zeros(1,length(n1));real(y1)],'r');%'color',[0.1,0.1,0.1]
%   绘制螺旋面
    h=surf(ones(1,1)*n,[zeros(1,length(n));imag(y)],[zeros(1,length(n));real(y)],...
        'facecolor',[0.75 0.75 0.75],'edgecolor','none');%'color',[0.1,0.1,0.1]
    alpha(0.3)
    % hp- 绘制动点P的起点
    hp=plot3(0,imag(y(1)),real(y(1)),'y.','markersize',15);
    hp1=plot3(n(end),imag(y(end)),real(y(end)),'g.','markersize',15);
    hp2=quiver3(0,0,0,0,imag(y(1)),real(y(1)),'y','linewidth',1);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]
%     hp3=quiver3(n(end),0,0,0,imag(y(end)),real(y(end)),'g','linewidth',2);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]

%    h=plot3(n,imag(y),real(y),'r');
%   绘制螺旋线
    h=plot3(n,imag(y),real(y),'r');
%   绘制实部投影直线
%     h=plot3(ones(2,1)*n,1*ones(2,N),[zeros(1,N);real(y)],'color',[0.8,0.8,0.8]);
%   绘制实部投影曲线
    h=plot3(n,1*ones(1,N),1.1*real(y),'b','linewidth',0.5);
%     h=plot3(n,1.1*ones(1,N),real(y),'k');
%   绘制虚部投影直线
%     h=plot3(ones(2,1)*n,[zeros(1,N);imag(y)],-1*ones(2,N),'color',[0.8,0.8,0.8]);
%   绘制虚部投影曲线
    h=plot3(n,1.15*imag(y),-1*ones(size(n)),'k','linewidth',0.5);
%     h=plot3(n,imag(y),-1.1*ones(size(n)),'k');

% 
% 修饰figure
 ylabel('Imaginary Part');%,'Rotation',-6
 set(gca,'XTick',-5:6:30)
set(gca,'XTickLabel',{'Shaft ','0','0.2','0.4','0.6','0.8','1.0'})
zlabel('Real Part');
xlabel('time')
% view(27.5,22) %设置视角
mesh(Z-5,X*0.99,Y*0.99,'facecolor',[0.9 0.9 0.9],'edgecolor','none')
% hold on
surf(Z1-5,X1,Y1,'facecolor','none','edgecolor',[0.3 0.3 0.3]);
camlight;
lighting gouraud;
% light
% lighting phong
material metal
% axis equal;
% view(60,15) 

ylim([-1.5 1.5])
zlim([-1.5 1.5])
hold off
grid on

%% 共振
% 轴段说明

% Plot enclosed cylinder
clear; clc

% Sample values
len = 5;     % height
ra = 1;   % radius
% Create constant vectors
tht = linspace(0,2*pi,200); z = linspace(0,len,40);
% Create cylinder
xa = repmat(ra*cos(tht),40,1); ya = repmat(ra*sin(tht),40,1);
za = repmat(z',1,200);
% To close the ends

X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0];
Z = [za; flipud(za); za(1,:)];
% 响应轨迹
% tht1=01:1/50:2*pi-1/50;
tht1=linspace(0,2*pi,30);
z1=linspace(0,len,10);
xa1=repmat(ra*cos(tht1),10,1);ya1=repmat(ra*sin(tht1),10,1);
za1=repmat(z1',1,30);
X1=[xa1*0; flipud(xa1); (xa1(1,:))*0];Y1 = [ya1*0; flipud(ya1); (ya1(1,:))*0];
Z1 = [za1; flipud(za1); za1(1,:)];

% 扭振
n=0:0.005:30-0.05;
N=length(n);
w=-0.5;
wt=1*sin(2*pi*0.1*n);

y=exp(j*2*pi*w*n+j*2.5*sin(2*pi*0.1*n));
n1=n(1:2:end);
y1=y(1:2:end);

figure
hold on;
%   h=plot3(ones(2,1)*n1,[zeros(1,length(n1));imag(y1)],[zeros(1,length(n1));real(y1)],'r');%'color',[0.1,0.1,0.1]
%   绘制螺旋面
%     h=surf(ones(1,1)*n,[zeros(1,length(n));imag(y)],[zeros(1,length(n));real(y)],...
%         'facecolor',[0.75 0.75 0.75],'edgecolor','none');%'color',[0.1,0.1,0.1]
   
    % hp- 绘制动点P的起点
    hp=plot3(0,0,0,'b.','markersize',15);
%     hp1=plot3(n(end),imag(y(end)),real(y(end)),'g.','markersize',15);
    hp2=quiver3(0,0,0,0,1,1,'b','linewidth',1);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]
%     hp3=quiver3(n(end),0,0,0,imag(y(end)),real(y(end)),'g','linewidth',2);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]

mesh(Z-5,X*0.99,Y*0.99,'facecolor',[0.9 0.9 0.9],'edgecolor','none')
 alpha(0.8)
% hold on
surf(Z1-5,X1,Y1,'facecolor','none','edgecolor',[0.3 0.3 0.3]);
camlight;
lighting gouraud;
% light
% lighting phong
material metal
% axis equal;
% view(60,15) 

ylim([-1.5 1.5])
zlim([-1.5 1.5])
hold off
grid on

%% 冲击
% 轴段说明

% Plot enclosed cylinder
clear; clc

% Sample values
len = 5;     % height
ra = 1;   % radius
% Create constant vectors
tht = linspace(0,2*pi,200); z = linspace(0,len,40);
% Create cylinder
xa = repmat(ra*cos(tht),40,1); ya = repmat(ra*sin(tht),40,1);
za = repmat(z',1,200);
% To close the ends

X = [xa*0; flipud(xa); (xa(1,:))*0]; Y = [ya*0; flipud(ya); (ya(1,:))*0];
Z = [za; flipud(za); za(1,:)];

% tht1=01:1/50:2*pi-1/50;
tht1=linspace(0,2*pi,30);
z1=linspace(0,len,10);
xa1=repmat(ra*cos(tht1),10,1);ya1=repmat(ra*sin(tht1),10,1);
za1=repmat(z1',1,30);
X1=[xa1*0; flipud(xa1); (xa1(1,:))*0];Y1 = [ya1*0; flipud(ya1); (ya1(1,:))*0];
Z1 = [za1; flipud(za1); za1(1,:)];

% 无扭振
n=0:0.005:30-0.05;
N=length(n);
w=-0.5;
wt=1*sin(2*pi*0.1*n);
% 测试用，不执行
y2x=exp(-0.05*n).*sin(2*pi*0.1*n);
y2y=exp(-0.05*n).*cos(2*pi*0.1*n);
y1=exp(j*2*pi*w*n+j*2.5*sin(2*pi*0.1*n));
%
% 绘图用
C=0.15;
wt=0.2;
A=2.5;
y2=A*exp(-C*n).*sin(2*pi*wt*n);
y=exp(j*2*pi*w*n+j*y2);

% 测试用
yc1=exp(j*2*pi*w*n);
% y=exp(j*2*pi*w*n+j*2.5*sin(2*pi*0.1*n));
y_1=exp(j*A*exp(-C*n).*sin(2*pi*wt*n));
y_r=y-yc1;
% 

figure
     h=plot3(n,imag(y),real(y),'r');
     hold on
     hc1=plot3(n,imag(yc1),real(yc1),'k')
      
 figure    
     h_r=plot3(n,imag(y_r),real(y_r))
     hold on
 figure
     h2=plot(n,y2)
 figure
     h_1=plot3(n,imag(y_1),real(y_1))
%

%y=exp(j*2*pi*w*n+j*2.5*exp(-3*n)*sin(2*pi*0.1*n));

n1=n(1:2:end);
y1=y(1:2:end);

figure
hold on;
%   h=plot3(ones(2,1)*n1,[zeros(1,length(n1));imag(y1)],[zeros(1,length(n1));real(y1)],'r');%'color',[0.1,0.1,0.1]
%   绘制螺旋面
     h=surf(ones(1,1)*n,[zeros(1,length(n));imag(y)],[zeros(1,length(n));real(y)],...
         'facecolor',[0.75 0.75 0.75],'edgecolor','none');%'color',[0.1,0.1,0.1]
alpha(0.3)   
    % hp- 绘制动点P的起点
    hp=plot3(0,0,0,'b.','markersize',15);
%     hp1=plot3(n(end),imag(y(end)),real(y(end)),'g.','markersize',15);
%    hp2=quiver3(0,0,0,0,1,1,'b','linewidth',1);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]
%     hp3=quiver3(n(end),0,0,0,imag(y(end)),real(y(end)),'g','linewidth',2);%[0,0],[0,0],[0 0],[0 1],[0 1],[0 1]
%    h=plot3(n,imag(y),real(y),'r');
%   绘制螺旋线
    h=plot3(n,imag(y),real(y),'r');
%   绘制螺旋面
%   绘制实部投影直线
%     h=plot3(ones(2,1)*n,1*ones(2,N),[zeros(1,N);real(y)],'color',[0.8,0.8,0.8]);
%   绘制实部投影曲线
    h=plot3(n,2*ones(1,N),1*real(y),'b','linewidth',0.5);
%     h=plot3(n,1.1*ones(1,N),real(y),'k');
%   绘制虚部投影直线
%     h=plot3(ones(2,1)*n,[zeros(1,N);imag(y)],-1*ones(2,N),'color',[0.8,0.8,0.8]);
%   绘制虚部投影曲线
    h=plot3(n,1*imag(y),-2*ones(size(n)),'k','linewidth',0.5);
%     h=plot3(n,imag(y),-1.1*ones(size(n)),'k');
mesh(Z-5,X*0.99,Y*0.99,'facecolor',[0.9 0.9 0.9],'edgecolor','none')

% hold on
surf(Z1-5,X1,Y1,'facecolor','none','edgecolor',[0.3 0.3 0.3]);
camlight;
lighting gouraud;
% light
% lighting phong
material metal
% axis equal;
% view(60,15) 
ylabel('Imaginary Part');%,'Rotation',-6
zlabel('Real Part');
xlabel('time')
xlim([-5 30])
ylim([-2 2])
zlim([-2 2])
hold off
grid on

%%
t=[0:1/100:1-1/100];
sin9=sin(2*pi*9*t);
sin10=sin(2*pi*9*t)+sin(2*pi*1*t);
sinx=sin10-sin9;
figure
subplot(311)
plot(sin10)
subplot(312)
plot(sin9)
subplot(313)
plot(sinx)


