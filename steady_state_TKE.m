%% Created by: Garett Warner
% Takes output of aver_time_series.F and 
clear
clc
close all
addpath('/glade/work/gwarner/Functions')
%% INPUT
ug=15; % geopstrophic wind along x
vg=0; % geopstrophic wind aong y
ug_mag=sqrt((ug^2)+(vg^2)); % geostrophic wind magnitude
Lz=4000; % vertical domain
tend=160000; % total timestep
f=8.57733e-5; % coriolis term
level=40; % height index for analysis
order=1e-3; % order for y-axis
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/'; % data location 
%% Time Domain
dt=0.1; % time step of simulation
base=1000; % intervels data are stored in aver.out files
time=(base*dt*(1:(tend/base)))/3600; % time domain (hours)
%% Reading in data
filename=["TKE_profile_timeseries.txt","bouy_prod_TKE_timeseries.txt","ushear_prod_TKE_timeseries.txt","vshear_prod_TKE_timeseries.txt","dissip_prod_TKE_timeseries.txt"]; %name of files
variables=["TKE","B","ushear","vshear","dissip"];
scale=[ug_mag^2,(ug_mag^3)/Lz,(ug_mag^3)/Lz,(ug_mag^3)/Lz,(ug_mag^3)/Lz];
figure_name1=['TKE_timeseries_level_',num2str(level),'.pdf'];
figure_name2=['dTKEdt_timeseries_level_',num2str(level),'.pdf'];
figure_name3=['bouy_prod_TKE_timeseries_level_',num2str(level),'.pdf'];
figure_name4=['ushear_prod_TKE_timeseries_level_',num2str(level),'.pdf'];
figure_name5=['vshear_prod_TKE_timeseries_level_',num2str(level),'.pdf'];
figure_name6=['dissip_prod_timeseries_level_',num2str(level),'.pdf'];
%% Reading in data
for j=1:length(filename)
    location=strjoin([filepath,filename(j)],'');
    data=load(location)*scale(j);
    assignin('base',variables(j),data);
end
%% Calculations
dTKEdt=(TKE(:,3:end)-TKE(:,1:end-2))./(2*base*dt);
%% Plotting
t=figure;
plot(time,TKE(level,:))
xlim([0 time(end)])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$\frac{1}{2}\overline{u^{\prime}_{i}u^{\prime}_{i}}$ (m$^{2}$ s$^{-2}$)','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/',figure_name1],'ContentType', 'vector', 'BackgroundColor', 'none');
%%
t=figure;
plot(time(2:end-1),dTKEdt(level,:))
xlim([0 time(end-1)])
ylim([-order order])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$\frac{1}{2}\partial\overline{u^{\prime}_{i}u^{\prime}_{i}}/\partial{t}$ (m$^{2}$ s$^{-3}$)','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/',figure_name2],'ContentType', 'vector', 'BackgroundColor', 'none');
%%
t=figure;
plot(time(2:end-1),B(level,2:end-1))
xlim([0 time(end-1)])
ylim([-order order])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$\frac{g}{\overline{\theta}}\overline{w^{\prime}\theta^{\prime}}$ (m$^{2}$ s$^{-3}$)','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/',figure_name3],'ContentType', 'vector', 'BackgroundColor', 'none');
%%
t=figure;
plot(time(2:end-1),ushear(level,2:end-1))
xlim([0 time(end-1)])
ylim([-order order])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$\overline{u^{\prime}w^{\prime}}\frac{\partial{\overline{u}}}{\partial{z}}$ (m$^{2}$ s$^{-3}$)','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/',figure_name4],'ContentType', 'vector', 'BackgroundColor', 'none');
%%
t=figure;
plot(time(2:end-1),vshear(level,2:end-1))
xlim([0 time(end-1)])
ylim([-order order])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$\overline{v^{\prime}w^{\prime}}\frac{\partial{\overline{v}}}{\partial{z}}$ (m$^{2}$ s$^{-3}$)','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/',figure_name5],'ContentType', 'vector', 'BackgroundColor', 'none');
%%
t=figure;
plot(time(2:end-1),dissip(level,2:end-1))
xlim([0 time(end-1)])
ylim([-order order])
xlabel('$t$ (hrs)','fontsize',12,'Interpreter','latex')
ylabel('$\epsilon$ (m$^{2}$ s$^{-3}$)','fontsize',12,'Interpreter','latex')
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/',figure_name6],'ContentType', 'vector', 'BackgroundColor', 'none');