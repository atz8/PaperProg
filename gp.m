function [P1,Pn]=gp(Gear_shape,Amp_theta,N)
%% Generate pulse signal
%   ���������������������岨��
%   Gear_shape ����Ϊ5��������������
%              ���ڱ�ʾһ�����ڵĲ��ηֶ�ʱ�䣨�����ݵ��ʾ��
%   Amp_theta  ����Ϊ5��������������
%              ��Gear_shape��Ӧ�������ֵ
%   N          ���������źŵ�����
T_theta=cumsum(Gear_shape-1); %�ۼ����
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
% X=[x1(1),x2(1)];% ���������xֵ
% Y=[x1(2),x2(2)];% ���������yֵ
%  p=polyfit(X,Y,1);%����ʽ��ϣ������1��ʾһ�׶���ʽ����ֱ��
% %�������ص�p�Ƕ�Ӧ����ʽ�������½���ϵ����һ����Ȼ��2��ϵ��
%  x=1:10;
%  y=polyval(p,x)%ͨ��p���Ӧx��yֵ
%  figure(2),plot(x,y)