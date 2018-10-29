% ≤‚ ‘shaft2D
% % Basic Parameters 
% ExRadius=ShaftData(:,1);    % Shaft InRadius 
% InRadius=ShaftData(:,2);  % Shaft ExRadius
% ShaftSectionLen=ShaftData(:,3);  %Shaft Section Length
% ShaftSectionRho=ShaftData(:,4);   %Shaft Section Density
% ShaftSectionE=ShaftData(:,5);     %Modulus of elasticity
% ShaftSectionMu=ShaftData(:,6);    %Poisson Ratio
ShaftData=[0.01 0 0.3 7900 2.1e11 0.3];
s=CrtShaftModel(ShaftData)
E=2.1e11;
mu=0.3;
T=0.01;
Ip=pi*ShaftData(1)^4/32;
G=E/(2*(1+mu));
G=8.0e10;
a=(180/3.14)*(T*ShaftData(3))/(G*Ip);
K=G*Ip/0.3*(pi/180)