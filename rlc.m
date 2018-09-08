function [U_cn1,i_Ln1,U_cn5]=rlc(RLC,init_Value,N,dt,T)
format compact
% N=RLC(1);
L=RLC(1);
C=RLC(2);
U_c0=init_Value(1);
i_L0=init_Value(2);
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
%  x_axis=[1:10].*ones(length(t));
 figure(7)
 plot(t,U_cn5)
 figure(6)
 hold on
 for c5=1:N
    plot3(t,c5*ones(1,length(t)),U_cn5(c5,:))
 end
    ylim([-0.5,1.5]) 