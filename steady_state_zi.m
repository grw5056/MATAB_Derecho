%% Created by: Garett Warner
% Takes output of aver_time_series.F and 
clear
clc
close all
addpath('/glade/work/gwarner/Functions')
%% INPUT
ug=15; % geostrophic wind magnitude
tend=160000; %total timesteps
Lz=4000; % vertical domain
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/'; % data location
%% Time Domain
dt=0.1; % time step of simulation
base=1000; % intervels data are stored in aver.out files
time=(base*dt*(1:(tend/base)))/3600; % time domain (hours)
%% Reading in data
filename=["zi_timeseries.txt","w_star_timeseries.txt"]; %name of files
variables=["zi","w_star"]; % what to name variables
scale=[Lz,ug]; % dimensional numbers
figure_name='steady_state_zi'; %save file
for j=1:length(filename)
    location=strjoin([filepath,filename(j)],'');
    data=load(location)*scale(j);
    assignin('base',variables(j),data);
end
%% time tendency
dzidt=(zi(3:end)-zi(1:end-2))./(2*base*dt);
%% plotting
t=figure;
set(t,'PaperPositionMode','auto')
plot(time(2:end-1),dzidt./w_star(2:end-1),'k')
hold on
yline(0.1,'--r')
xlim([0 time(end-1)])
ylim([0 1])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$(\partial{z_{i}}/\partial{t})/w_{\ast}$','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
drawnow;
exportgraphics(t,[filepath,'figures/',figure_name,'.pdf'],'ContentType', 'vector', 'BackgroundColor', 'none');
