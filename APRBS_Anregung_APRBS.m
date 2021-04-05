clear all;
clc;

format long

figure;

%--------------------------------------------------------------------------
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultTextFontSize',20);
set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');
set(0,'DefaultAxesXLim',[-1 1]);
set(0,'DefaultAxesYLim',[-1 0.6]);
%--------------------------------------------------------------------------


%INPUT
              
T=0.01;

%Zeitbereich des APRBS-Signals in Sekunden
t_min = T;
t_max = 25*T; 

%Minimale und maximale Amplitude der Trainingsdaten

time_max = 10000;

tp = tp_input(time_max,[t_min,t_max],[0,0.4],T); %kg/s
tp2 = tp_input(time_max,[t_min,t_max],[0,1],T); 


sim('Thermal_Energy',[0.0 time_max]); %Hz (0->60Hz)

hold on;
plot3(data.signals.values(:,2),data.signals.values(:,3),data.signals.values(:,1),'+',...
                  'MarkerEdgeColor','k',...
                  'MarkerSize',8)



