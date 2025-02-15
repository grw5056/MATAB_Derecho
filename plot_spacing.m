clc;
clear;
close all;
%%
tini=111000;
tend=160000;
time_100=(tini:100:tend-100)*0.1/60/60;
time_1000=(tini:1000:tend-1000)*0.1/60/60;
time_6000=(tini:6000:tend-6000)*0.1/60/60;
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/figures/r_snapshots/';
filenames=["correlation_spacing_z16_t1_k120.txt","correlation_spacing_z16_t1_k24.txt","correlation_spacing_z16_t1_k6.txt","correlation_spacing_z16_t1_k4.txt"];
name=["c_t1_k120","c_t1_k24","c_t1_k6","c_t1_k4"];
for n=1:length(filenames)
    data=load(strjoin([filepath,filenames(n)],''));
    assignin('base',name(n),data)
    variance=var(data);
    assignin('base',strjoin([name(n),'_var'],''),variance)
end
%
filenames=["correlation_spacing_z16_t10_k120.txt","correlation_spacing_z16_t10_k24.txt","correlation_spacing_z16_t10_k6.txt"];
name=["c_t10_k120","c_t10_k24","c_t10_k6"];
for n=1:length(filenames)
    data=load(strjoin([filepath,filenames(n)],''));
    assignin('base',name(n),data)
    variance=var(data);
    assignin('base',strjoin([name(n),'_var'],''),variance)
end
%
filenames=["correlation_spacing_z16_t100_k120.txt","correlation_spacing_z16_t100_k24.txt","correlation_spacing_z16_t100_k6.txt"];
name=["c_t100_k120","c_t100_k24","c_t100_k6"];
for n=1:length(filenames)
    data=load(strjoin([filepath,filenames(n)],''));
    assignin('base',name(n),data)
    variance=var(data);
    assignin('base',strjoin([name(n),'_var'],''),variance)
end
%
filenames=["correlation_spacing_z16_t1000_k120.txt","correlation_spacing_z16_t1000_k24.txt","correlation_spacing_z16_t1000_k6.txt"];
name=["c_t1000_k120","c_t1000_k24","c_t1000_k6"];
for n=1:length(filenames)
    data=load(strjoin([filepath,filenames(n)],''));
    assignin('base',name(n),data)
    variance=var(data);
    assignin('base',strjoin([name(n),'_var'],''),variance)
end
%
filenames=["correlation_spacing_z16_t6000_k120.txt","correlation_spacing_z16_t6000_k24.txt","correlation_spacing_z16_t6000_k6.txt"];
name=["c_t6000_k120","c_t6000_k24","c_t6000_k6"];
for n=1:length(filenames)
    data=load(strjoin([filepath,filenames(n)],''));
    assignin('base',name(n),data)
    variance=var(data);
    assignin('base',strjoin([name(n),'_var'],''),variance)
end
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/figures/r_snapshots/';
filenames=["correlation_spacing_rot_z16_t1_k120.txt","correlation_spacing_rot_z16_t1_k24.txt","correlation_spacing_rot_z16_t1_k6.txt","correlation_spacing_rot_z16_t1_k4.txt"];
name=["c_t1_k120_rot","c_t1_k24_rot","c_t1_k6_rot","c_t1_k4_rot"];
for n=1:length(filenames)
    data=load(strjoin([filepath,filenames(n)],''));
    assignin('base',name(n),data)
end

%%
t=figure;
plot(time_100,c_t1_k120,'r-')
hold on
plot(time_100,c_t1_k24,'g-')
hold on
plot(time_100,c_t1_k6,'b-')
hold on
plot(time_100,c_t1_k4,'k-')
xlim([time_100(1) time_100(end)])
ylim([2.8 3.1])
xlabel('Time (hrs)','Interpreter','latex')
ylabel('Linear Updraft Spacing (km)','Interpreter','latex')
legend({'$k$ = 120','$k$ = 24','$k$ = 6','$k$ = 4'},'Interpreter','latex')
exportgraphics(t,[filepath,'/../filtering_updraft_spacing.pdf'],'Resolution',300);
%%
t=figure;
plot(time_100,c_t1_k120_rot,'r-')
hold on
plot(time_100,c_t1_k24_rot,'g-')
hold on
plot(time_100,c_t1_k6_rot,'b-')
hold on
plot(time_100,c_t1_k4_rot,'k-')
xlim([time_100(1) time_100(end)])
ylim([2.8 3.1])
xlabel('Time (hrs)','Interpreter','latex')
ylabel('Linear Updraft Spacing (km)','Interpreter','latex')
legend({'$k$ = 120','$k$ = 12','$k$ = 6','$k$=4'},'Interpreter','latex')
exportgraphics(t,[filepath,'/../filtering_updraft_spacing.pdf'],'Resolution',300);
%%
t=figure;
plot(time_100,c_t1_k120,'k')
hold on
plot(time_100,c_t10_k120,'r')
hold on
plot(time_100,c_t100_k120,'g')
hold on
plot(time_1000,c_t1000_k120,'b')
hold on
plot(time_6000,c_t6000_k120,'c')
xlim([time_100(1) time_100(end)])
ylim([2.8 3.1])
xlabel('Time (hrs)','Interpreter','latex')
ylabel('Linear Updraft Spacing (km)','Interpreter','latex')
legend({'$\delta t$ = 1','$\delta t$ = 10','$\delta t$ = 100','$\delta t$ = 1000','$\delta t$ = 6000'},'Interpreter','latex')
exportgraphics(t,[filepath,'/../averaging_updraft_spacing.pdf'],'Resolution',300);
%%
% t=figure;
% plot(time_100,f_t1_k120,'k')
% hold on
% plot(time_100,f_t10_k120,'r')
% hold on
% plot(time_100,f_t100_k120,'g')
% hold on
% plot(time_1000,f_t1000_k120,'b')
% hold on
% plot(time_6000,f_t6000_k120,'c')
% xlim([time_100(1) time_100(end)])
% ylim([2.8 3.1])
% xlabel('Time (hrs)','Interpreter','latex')
% ylabel('Linear Updraft Spacing (km)','Interpreter','latex')
% legend({'$\delta t$ = 1','$\delta t$ = 10','$\delta t$ = 100','$\delta t$ = 1000','$\delta t$ = 6000'},'Interpreter','latex')
% exportgraphics(t,[filepath,'/../fourier_updraft_spacing.pdf'],'Resolution',300);