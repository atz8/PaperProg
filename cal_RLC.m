%%
clear;format compact
%% 二阶系统参数与初值设置
L=0.5;C=0.05;
U_c0=1;i_L0=0;
%% 仿真参数设置
N=10;dt=0.01;T=1;
t=0:dt:T;
U_cn1=zeros(N,length(t));
i_Ln1=zeros(N,length(t)-1);
% mag=zeros(N,)
%% 仿真二阶系统零响应
for R=1:N
    alpha=R/2/L; %阻尼系数
    w_n=sqrt(1/(L*C));  %固有频率
    p1=-alpha-sqrt(alpha^2-w_n^2);  %特征根1
    p2=-alpha+sqrt(alpha^2-w_n^2);  %特征根2
    
    % 系统分子，分母
    num=[U_c0,R/L*U_c0+i_L0/C];
    den=[1,R/L,1/L/C];
    [r,p,k]=residue(num,den);
    H=tf(num,den)
%     [mag(R,:),phase(R,:),w]=bode(H);
    [mag,phase,w]=bode(H);
    % 
    U_cn=r(1)*exp(p(1)*t)+r(2)*exp(p(2)*t);
    i_Ln=C*diff(U_cn)/dt;

    U_cn1(R,:)=U_cn;
    i_Ln1(R,:)=i_Ln;
end

U_cn2=U_cn1;
U_cn3=zeros(size(U_cn2));
U_cn4=0.001*randn(size(U_cn2));

for c1=1:N
    for c2=1:length(t)
        if U_cn1(c1,c2)<=0;
            U_cn2(c1,c2:end)=0;
            U_cn3(c1,c2:end)=U_cn4(c1,c2:end);
        end
    
    end
end

U_cn2n=U_cn2+0.001*randn(size(U_cn2));
    figure(1)
    plot(t,U_cn1)

    figure(2)
    plot(t(2:end),i_Ln1)

    figure(3)
    plot(t,U_cn2)
    
    figure(4)
    plot(t,U_cn2n)

 U_cn3n=U_cn2+U_cn3;   
    figure(5)
    plot(t,U_cn3n)
    hold on
 
 U_cn4=U_cn3;  
 for c3=1:R
     for c4=1:length(t)
         if U_cn3(c3,c4)<0;
            U_cn4(c3,c4)=1;
         end 
     end
 end
 
 U_cn5=U_cn4+U_cn2;
%  x_axis=[1:10]*ones(length(t));
 figure(7)
 plot(t,U_cn5)
 figure(6)
 hold on
 for c5=1:N
    plot3(t,c5*ones(1,length(t)),U_cn5(c5,:))
 end
    ylim([-0.5,1.5]) 
 
%% 梯形脉冲，Uc响应，比较器输出
A_theta=[2,1,1,2,2];
G_shape=[1,15,35,15,35];
N=5;
[P1,P10]=gp(G_shape,A_theta,N);
% 增加正向充电和反向充电的点Vcc1 和 Vcc2
Vcc1_x=[4:100:404];
Vcc2_x=[54:100:454];
Vcc1_y=P10(Vcc1_x);
Vcc2_y=P10(Vcc2_x);
figure
plot(P10,'k')
hold on
% 增加正向充电和反向充电的点Vcc1 和 Vcc2
plot(Vcc1_x,Vcc1_y,'b.');
plot(Vcc2_x,Vcc2_y,'r.');
ylim([0.5,2.5])

%% cal_RLC
L=0.5;C=0.05;
U_c0=0.5;i_L0=0;
RLC=[L,C];
init_Value=[U_c0,i_L0];
[a,b,c]=rlc(RLC,init_Value,10,0.01,1);
figure
for c1=1:10
    subplot(5,2,c1)
    plot(c(c1,:));
    xlim([0,100])
    ylim([-0.5,1.5])
end


figure
for c1=1:10
    subplot(5,2,c1)
    plot(b(c1,:));
    xlim([0,100])
    ylim([-0.5,1.5])
end

fs_Gear=1000;
fs_Uc=2000;
temp_t=[0:1/fs_Uc:0.5-1/fs_Uc];
temp_t1=[0:1/fs_Gear:0.5-1/fs_Gear];

Vcc_x=[4:50:length(P10)];
Vcc_y=P10(Vcc_x);
dVcc=4/fs_Gear;


figure
%%% fig311
subplot(311)
plot(temp_t1,P10,'b','linewidth',2)
hold on
% 增加正向充电和反向充电的点Vcc1 和 Vcc2
plot(Vcc1_x/1005,Vcc1_y,'k.','markersize',14);
plot(Vcc2_x/1005,Vcc2_y,'r.','markersize',14);
ylim([0.9,2.1])
hold off

subplot(312)
%%% fig 312
% plot(temp_t,-c(10,:)+1.5)

plot(temp_t1,1.5*ones(1,length(temp_t1)),'-.')
hold on
for c2=1:2*N
%     plot([Vcc1_x(c2),Vcc1_x(c2)]/fs_Gear,[1,1.5])
    if mod(c2,2)~=0
        plot([Vcc(c2),Vcc(c2)]/fs_Gear,[1,c(10,end-1)+1.5],'k','linewidth',2)
        plot(temp_t(1+(c2-1)*100:c2*100)+dVcc,-c(10,1:end-1)+1.5,'k','linewidth',2)  %负向充电
    else
        plot([Vcc(c2),Vcc(c2)]/fs_Gear,[2,-c(10,end-1)+1.5],'r','linewidth',2)
        plot(temp_t(1+(c2-1)*100:c2*100)+dVcc,c(10,1:end-1)+1.5,'r','linewidth',2)     %正向充电    
    end
end
hold off
xlim([0,0.5])
ylim([0.9,2.1])

subplot(313)
%%% fig313
xlim([0,0.5])
ylim([-0.2,1.2])
hold on
for c2=1:2*N
%     plot([Vcc1_x(c2),Vcc1_x(c2)]/fs_Gear,[1,1.5])
    if mod(c2,2)~=0
        plot([Vcc(c2),Vcc(c2)]/fs_Gear,[0,1],'k','linewidth',2)
        plot(temp_t(1+(c2-1)*100:c2*100)+dVcc,zeros(length(c)-1),'k','linewidth',2)  %负向充电
    else
        plot([Vcc(c2),Vcc(c2)]/fs_Gear,[0,1],'r','linewidth',2)
        plot(temp_t(1+(c2-1)*100:c2*100)+dVcc,ones(length(c)-1),'r','linewidth',2)     %正向充电    
    end
end
hold off
box on

%%  响应特性













    
    