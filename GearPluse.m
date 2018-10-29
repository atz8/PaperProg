%计算齿轮脉冲

%% 创建齿轮
%%
clear;clc;
z=10;  %齿数
r_h=100; %
r_l=80;  %
fs=360;
theta=fs/z;
theta_h=15;
% len_h=13;
% len_l=13;
theta_l=theta_h;
theta_m=theta-(theta_l+theta_h);
alpha=0:pi/180:2*pi-1/360;

x_h=r_h*cos(alpha);
y_h=r_h*sin(alpha);
rho=r_h*sin(alpha);
x_h08=0.8*x_h;
y_h08=0.8*y_h;

figure
hold on
plot(x_h(1:15),y_h(1:15),'.','markersize',3)
plot(x_h(16:30),y_h(16:30),'b.')

xlim([-110,110])
ylim([-110,110])

%
g_x=zeros(z,36);
g_y=zeros(z,36);

for c=0:z-1
    g_x(c+1,:)=x_h(c*36+1:(c+1)*36);
    g_y(c+1,:)=y_h(c*36+1:(c+1)*36);
end
% 加入旋转
M=[1,0,0;0,1,0;0,0,1];
Rot(1,:)=x_h;
Rot(2,:)=y_h;
Rot(3,:)=1;
Rot_30=M*Rot;
M1=[cos(pi/6)  sin(pi/6)  0;
        -sin(pi/6) cos(pi/6)  0;
             0            0           1];
Rot2_30=M1*Rot_30;
% 重绘齿轮
x_h1=Rot2_30(1,:);
y_h1=Rot2_30(2,:);
g_x1=zeros(z,36);
g_y1=zeros(z,36);

for c=0:z-1
    g_x1(c+1,:)=x_h1(c*36+1:(c+1)*36);
    g_y1(c+1,:)=y_h1(c*36+1:(c+1)*36);
end
%%
% 绘制齿轮
figure
for c1=1:z
    plot(g_x(c1,1:12),g_y(c1,1:12),'r','linewidth',1,'color','r')
    hold on
    plot(g_x(c1,1:2:12),g_y(c1,1:2:12),'.','markersize',10,'color','b')
    plot(g_x(c1,18:31)*0.8,g_y(c1,18:31)*0.8,'b','linewidth',1,'color','k')
    line([g_x(c1,12),g_x(c1,18)*0.8],[g_y(c1,12),g_y(c1,18)*0.8],'linewidth',1)
    
end
for c1=1:z-1
    line([g_x(c1,31)*0.8,g_x(c1+1,1)],[g_y(c1,31)*0.8,g_y(c1+1,1)],'linewidth',1)
end
line([g_x(10,30)*0.8,g_x(1,1)],[g_y(10,30)*0.8,g_y(1,1)],'linewidth',1)
xlim([-110,110])
ylim([-110,110])

% 齿轮旋转重绘
figure
for c1=1:z
    plot(g_x1(c1,1:12),g_y1(c1,1:12),'r','linewidth',2,'color','r')
    hold on
%     plot(g_x1(c1,1:2:12),g_y1(c1,1:2:12),'.','markersize',10,'color','b')
    plot(g_x1(c1,18:31)*0.8,g_y1(c1,18:31)*0.8,'b','linewidth',2,'color','k')
    line([g_x1(c1,12),g_x1(c1,18)*0.8],[g_y1(c1,12),g_y1(c1,18)*0.8],'linewidth',2)
    
end
for c1=1:z-1
    line([g_x1(c1,31)*0.8,g_x1(c1+1,1)],[g_y1(c1,31)*0.8,g_y1(c1+1,1)],'linewidth',2)
end
line([g_x1(10,30)*0.8,g_x1(1,1)],[g_y1(10,30)*0.8,g_y1(1,1)],'linewidth',2)
xlim([-110,110])
ylim([-110,110])

%%


