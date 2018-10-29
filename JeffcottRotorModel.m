%% Script Information 
%--------------------------------------------------------------------------
% 20181007 Rotor Dynamic Jeffcott Model
%  writen by atz8
%  ref. Zhang Wen.Basic of Rotar Dynamic M.Science Press.1995
%--------------------------------------------------------------------------
%
% Create static coodinate system o-xyz
% Definition  variable :'delta' is deviation of centroid
% Definition variable :'theta' is the angular velocity of 'a' axis of dynamic system 
delta=[yt;zt];
P=complex(yt,zt);
% Create dynamic coodinate system o-abc
delta=bt*b+ct*c;
Rot_Matrix=[cos(theta) sin(theta);-sin(theta) cos(theta)];
[b;c]=Rot_Matrix*delta;

