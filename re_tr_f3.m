%% re_tr_f3.m
% Autor: Marcelo e Andre
% Last modified: 02/12/2020
% Funcao para o codigo retrsim.m

function dBdt = re_tr(t,Bi,flag,parm)
% param = [EEj PBj];
EE = parm (1);
PB = parm (2);
M0 = PB*(1 - EE)/Bi;      
dBdt = PB*Bi - M0*Bi;