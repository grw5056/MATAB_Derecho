clc;
clear;
close all;
%%
addpath('/glade/work/gwarner/Functions')
load('w_inst_color.txt')
%%
%filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
filepath='/glade/derecho/scratch/gwarner/het_ug8_lx6_phi90/output/';
%filename='w_slice_output/w_z16_t0001_k120_01.nc';
filename='w_z32_snap_k006.nc';
ug=8;
Lx=6000;
Ly=6000;
nx=240;
ny=240;
avg_width=1;
filter=6;
lvl=32;
%%
w=ncread([filepath,filename],'w')*ug;
timesteps=ncread([filepath,filename],'time');
%%
xtick_input_w=[(Lx/nx/1000),2,4,6,8,10,12];
xticklabel_input_w={'0','2','4','6','8','10','12'};
ytick_input_w=[(Ly/ny/1000),2,4,6,8,10,12];
yticklabel_input_w={'0','2','4','6','8','10','12'};
%% Spacing Domain
x=Lx/nx:Lx/nx:Lx;
y=Ly/ny:Ly/ny:Ly;
[X,Y]=meshgrid(x,y);
%%
for k=1:length(timesteps)
    t=figure;
    %if k==1
    %    t=figure;
    %else
    %    t=figure('Visible','off');
    %end
    contourf(X/1000,Y/1000,squeeze(w(:,:,k))','Linecolor','none')
    pbaspect([1 1 1])
    xlabel('$x$ (km)','Interpreter','latex','FontSize',12)
    ylabel('$y$ (km)','Interpreter','latex','FontSize',12)
    xticks(xtick_input_w)
    xticklabels(xticklabel_input_w)
    yticks(ytick_input_w)
    yticklabels(yticklabel_input_w)
    clim([-10 10])
    colormap(w_inst_color)
    colorbar()
    text(12.6,12.4,'m s$^{-1}$','FontSize',12,'Interpreter','Latex')
    text(0.01,1.02,['timestep = ',num2str(timesteps(k))],'units','normalized')
    set(gca,'Fontsize',12)
    exportgraphics(t,[filepath,'/figures/w_snapshots/w_z',num2str(lvl),'_t',num2str(avg_width),'_k',num2str(filter),'_',num2str(timesteps(k)),'.png'],'Resolution',300);
end