function [P1,Pn]=gp(Gear_shape,Amp_theta,N)
%% Generate pulse signal
%   本函数用于生成梯形脉冲波形
%   Gear_shape 参数为5个参数的向量，
%              用于表示一个周期的波形分段时间（用数据点表示）
%   Amp_theta  参数为5个参数的向量，
%              与Gear_shape对应的脉冲幅值
%   N          生成脉冲信号的周期
T_theta=cumsum(Gear_shape-1); %累加求和
% Amp_theta=[1,2,2,1,1];

lx=[];
ly=[];
P1=[];
figure
hold on

for c=1:length(Amp_theta)-1
    
    cx=[T_theta(c),T_theta(c+1)];
    cy=[Amp_theta(c),Amp_theta(c+1)];
    cl=polyfit(cx,cy,1);
    lx=[T_theta(c):T_theta(c+1)];
    ly=polyval(cl,lx);
    P1=[P1,ly];
    plot(lx,ly)
    
end
ylim([0.5,2.5])

Pn=[];
for c2=1:N
    Pn=[Pn,P1];
end

figure
plot(Pn)
ylim([0.5,2.5])


% l1=k(1)*t(1:T_theta(1));
% l2=k(2)*t(T_theta(1):T_theta(2))+l1(end);
% l3=(l2(end)-k(3)*T_theta(2)/fs)+k(3)*t(T_theta(2):T_theta(3));
% l4=k(4)*t(T_theta(3):T_theta(4));


% %%
% x1=[1,2];
% x2=[5,8];
% X=[x1(1),x2(1)];% 两点坐标的x值
% Y=[x1(2),x2(2)];% 两点坐标得y值
%  p=polyfit(X,Y,1);%多项式拟合，后面的1表示一阶多项式，即直线
% %函数返回的p是对应多项式按次数下降的系数，一阶显然有2个系数
%  x=1:10;
%  y=polyval(p,x)%通过p求对应x的y值
%  figure(2),plot(x,y)