%% Solucao do Rede trofica
% Autor: Marcelo R Errera 23/11/2020 
%        Andre Roque
% Last revised: 02/12/2020
% Versao: 1.0 crescimento do fito-plancton
%
% Observacoes:
% 3.resolver pelo ODE45 MATLAB
% 4.Utiliza as funcoes: re_tr_f2.m (Apendice I)
% 5. vide artigo:  
%
% Alteracao 1.1 em 12.02.2020 
% adicao da mortalidade emperica pelo coeficiente de von Bertalanffi
% novos paramentros de k, Tc e Linf, adicao de DCji reparametrizado
%% -------------------------------------------------------------------
%% 
clc;
clear;

% Abertura
fprintf('\nSolucao do Sistema Rede Trofica Aquatica')
%
tmax = 15;
h = 0.1; % time step
Ntmax = tmax/h; 
%% Condicao Inicial, tempo de simulacao e passo de tempo:
% [B(t=0),tmax,h
Bj = zeros (Ntmax, 1); % densidade relativa [ton/km2]
%
Bj(1) = 31.66 % Fitoplancton / densidade relativa [ton/km2] 
CI = Bj(1);
%% Parametros:
% Taxa de crescimento por grupo [tonNB/tonB/ano]
%gj = 0.250066; % Fitoplancton
%
% Dieta / consumo fracao da dieta/consumo do predador (j) de presas (i)
%         quanto a presa (i) contribui para a dieta do predador (j)
%        de biomassa existente de (j) no instante t, [1/ano]
% QBj = 75.94;   % [1/ano] Fitoplancton
% Qji consumo do grupo (i) sobre o grupo (j) ton(i)/km2/ano
% Qij consumo de (j) por todos grupos (i)  ton(j)/km2/ano
% Qij = 

% % PBj  - produtividade do grupo (j) pela biomassa total tonP/tonB/ano
PBj = 183; % fitoplacnton
% PBj (1) = 1.5; % Hopmal J

% % EEj - eficiencia ecotrofica do grupo (j)  ( 0 <= EE <= 1)
EEj = 0.29; % Fitoplancton

%% Solucao do sistema (solver ODE45)
%[Ngrp, g[..], Q_ji [..], Q_ij, M0]
param = [EEj PBj];
[t, Bj] = ode45 ('re_tr_f2', [0:h:tmax], CI, [], param);

%% Resultados
%% figure 1
figure(1); clf;
semilogy(t,Bj,'b-.');
legend('fitoplancton');
title('Evolucao Temporal da Biomassa do Reservatorio','Fontsize',16)
xlabel('Tempo [anos]','fontsize',14)

%% figure 2

%output = [(0:tmax)',Iobs]; 
%save('Iobs.txt',output','-ascii')

