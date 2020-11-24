%% re_tr_f.m
% Autor: Marcelo e Andre
% Last modified: 23/11/2020
% Funcao para o codigo retrsim.m

function dBdt = re_tr(t,Bi,flag,parm)
%[Ngrp, g[..], Q_ji [..], Q_ji [..], M0 [..]]

% por um filtro para evitar pop negativa
Ngrp = parm (1);
% for i = 1:Ngrp,
%   if ( Bi(i) <= 1.e-4 ), 
%       Bi(i) = 0.;
%   end; % if
% end;
g(1:Ngrp) = parm (2:Ngrp+1);
Q_ji(1:Ngrp) = parm (Ngrp+2:2*Ngrp+1);
Q_ij(1:Ngrp) = parm (2*Ngrp+2:3*Ngrp+1);
M0(1:Ngrp) = parm (3*Ngrp+2:4*Ngrp+1);
%
for i = 1:Ngrp,
    dBdt(i) = g(i)*Q_ji(i) - Q_ij(i) - M0(i)*Bi(i)
end;

dBdt = dBdt(:);