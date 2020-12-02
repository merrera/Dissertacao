%% Solucao do Rede trofica
% Autor: Marcelo R Errera 23/11/2020 
%        Andre Roque
% Last revised: 02/12/2020
% Versao: 1.1 
%
% Observacoes:
% 3.resolver pelo ODE45 MATLAB
% 4.Utiliza as funcoes: re_tr.m (Apendice I)
% 5. vide artigo:  
%
% Alteracao 1.1 em 25.11.2020 
% adicao da mortalidade emperica pelo coeficiente de von Bertalanffi
% novos paramentros de k, Tc e Linf, adicao de DCji reparametrizado
%% -------------------------------------------------------------------
%% 
clc;
clear;

% Abertura
fprintf('\nSolucao do Sistema Rede Trofica Aquatica')

%% Condicao Inicial, tempo de simulacao e passo de tempo:
% [B(t=0),tmax,h
%Nt= ; % tamanho da populacao
Ngrupos = 2;
Bj = zeros (1, Ngrupos); % densidade relativa [ton/km2]
gj = zeros (1, Ngrupos); % Taxa de crescimento por grupo [tonNB/tonB/ano]
DCji = zeros (Ngrupos, Ngrupos); % Dieta / consumo fracao da dieta/consumo do predador (j) de presas (i)
DCij = zeros (Ngrupos, Ngrupos); % Dieta / consumo fracao da dieta/consumo do predador (i) de presas (j)
Qji = zeros (1, Ngrupos);
Qij = zeros (1, Ngrupos);
QBj = zeros (1, Ngrupos);
M0j = zeros (1, Ngrupos);
PBj = zeros (1, Ngrupos);
%
Bj(1) = 1.66; % Hopmal J / densidade relativa [ton/km2]
Bj(2) = 3.46; % Hopmal A / densidade relativa [ton/km2]
CI = [Bj(1) Bj(2)]; % initial condition
%
tmax = 15;
h = 0.1; % time step
%% Parametros:
% Taxa de crescimento por grupo [tonNB/tonB/ano]
gj(1) = 0.119713; % Hopmal J
gj(2) = 0.157895; % Hopmal A
%
% Dieta / consumo fracao da dieta/consumo do predador (j) de presas (i)
%              quanto a presa (i) contribui para a dieta do predador (j)
DCji = [ 0.	0.;
         0. 0.167 ]; %  4 e 5
% % Adicao de DCji reparametrizado      
% DCji = [ 0.	0.;
%          0. 1. ]; %  4 e 5     
%
% DCji = [ 0	0	0	0	0	0	0	0	0;
%          0	0	0	0	0	0	0	0	0;
%          0	0	0	0	0	0	0	0	0;
%          0.11437	0	0	0	0	0	0	0	0;
%          0	0.130682	0	0.167	0	0	0	0	0;
%          0,539589	0,625	0,0128	0,667	0,25	0	0	0	0;
%          0,322581	0,125	0	0,083	0,0625	0	0	0	0;
%          0,02346	0,099432	0,6422	0,083	0,3125	0	0,5	0	0;
%          0	0,019886	0,345	0	0,375	1	0,5	1	0];
% (Q/B)j / taxa de consumo de biomassa da especie (j) pela quantidade total
%          de biomassa existente de (j) no instante
QBj(1) = 12.53;
QBj(2) =  5.70;
% Qji consumo do grupo (i) sobre o grupo (j) ton(i)/km2/ano
for j = 1:Ngrupos,   
    for i = 1:Ngrupos,
        Qji (j) = Qji(j)+ Bj(j)*QBj(j)*DCji(j,i);
    end; % i
end; % j
% Qij consumo de (j) por todos grupos (i)  ton(j)/km2/ano
for j = 1:Ngrupos,   
    for i = 1:Ngrupos,
        Qij (j) = Qij(j)+ Bj(i)*QBj(i)*DCji(i,j);
    end; % i
end; % j
% % PBj  - produtividade do grupo (j) pela biomassa total tonP/tonB/ano
% PBj (1) = 1.5; % Hopmal J
% PBj (2) = 0.9; % Hopmal A
% % EEj - eficiencia ecotrofica do grupo (j)  ( 0 <= EE <= 1)
% EEj (1) = 0.723; % Hopmal J
% EEj (2) = 0.138; % Hopma A
%
% M0j  - mortalidade espontenea de biomassa do grupo (j) tonBM/tonB/ano 
% for j = 1:Ngrupos,
%      M0j (j) = PBj(j)*(1 - EEj(j))/Bj(j);  
% end;
%
% coeficientes de crescimento de von Bertalanffi (VBGF)
k(1) = 0.15 % por ano para Hopmal J
k(2) = 0.15 % por ano para Hopmal A
%
Linf (1) = 68.25; % comprimento total em cm Hopmal J
Linf (2) = 101.85; % comprimento total em cm Hopmal A
%
Tc = 17; % temperatura media da agua em oC
%
% M0j  - mortalidade total do grupo (j) tonBM/tonB/ano 

for j = 1:Ngrupos,
    M0j (j) = (k(j)^0.65) * (Linf(j)^-0.279) * (Tc^0.463);
end;

%% Solucao do sistema (solver ODE45)
%[Ngrp, g[..], Q_ji [..], Q_ij, M0]
param = [Ngrupos, gj, Qji, Qij, M0j];
[t, Bj] = ode45 ('re_tr_f2', [0:h:tmax], CI, [], param);

%% Resultados
%% figure 1
figure(1); clf;
plot(t,Bj(:,1),'g-.', t,Bj(:,2),'r-');
legend('Hopmal J','Hopmal A');
title('Evolucao Temporal da Biomassa do Reservatorio','Fontsize',16)
xlabel('Tempo [anos]','fontsize',14)
ylabel('B [ton/km2]','fontsize',14)

%% figure 2

%% Espaco de Fase:
% figure(4);
% plot(N(:,1), Ninftotal); 
% xmax = max ( N(:,1));
% [Sc,i] = max( Ninftotal );
% axis ( [ 0 1.25*xmax 0 1.25*Sc ] );
% hold on;
% plot( N(1,1),Ninftotal(1),'k*');
% plot( N(end,1), Ninftotal(end),'ko');
% plot( N(i,1), Sc, 'o');
% title('SEIR - Espaco de Fase','Fontsize',16)
% xlabel('Susceptiveis (S)','fontsize',14)
% ylabel('Infectados Totais (E+I)','fontsize',14)
% hold off;

%output = [(0:tmax)',Iobs]; 
%save('Iobs.txt',output','-ascii')

