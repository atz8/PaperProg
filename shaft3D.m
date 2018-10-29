function Tor_ang=shaft3D(ShaftData,T,G,I_p)
%% Script Information 
%--------------------------------------------------------------------------
% 20181011 Rotor Dynamic Torsional Vibration
%  writen by atz8
%  ref. Zhang Wen.Basic of Rotar Dynamic M.Science Press.1995
%--------------------------------------------------------------------------
InRadius=ShaftData(:,1);
ExRadius=ShaftData(:,2);
ShaftSectionLen=ShaftData(:,3);
ShaftNum=length(ShaftSectionLen);
for c=1:ShaftNum-1
    

end