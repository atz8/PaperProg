function Shaft3D_Data=PlotShaft3D(ShaftData)
%% Script Information 
%--------------------------------------------------------------------------
% 20181029 Rotor Dynamic Torsional Vibration
%  writen by atz8
%  Plot Shaft Model with 3D 
%  ref. [1].Zhang Wen.Basic of Rotar Dynamic M.Science Press.1995
%        [2].Tang Xikuan. Dynamic of Machinery M.Higher Eudcation Press.1984
%--------------------------------------------------------------------------
% Output Parameter Description
% ShaftModelData 
% Input Parameter Description
% ShaftData
%--------------------------------------------------------------------------
ExRadius=ShaftData(:,1);    % Shaft InRadius 
InRadius=ShaftData(:,2);  % Shaft ExRadius
ShaftSectionLen=ShaftData(:,3);  %Shaft Section Length






