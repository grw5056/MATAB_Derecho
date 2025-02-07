clc;
clear;
%close all;
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
%filepath='/glade/derecho/scratch/gwarner/het_ug8_lx6_phi90/output/';
filename='w_slice_output/r_z16_t0001_k120.nc';
%filename='r_z32_snap_k120.nc';
ug=15;
Lx=12000;
Ly=12000;
nx=240;
ny=240;
avg_width=1;
lvl=16;
filter=120;
neg_angles=1;
contour_levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
use_levels=[0.2,0.3,0.4,0.5,0.6];
ABL_depth=1100;
%%
r=ncread([filepath,filename],'r');
timesteps=ncread([filepath,filename],'time');
%%
dx=Lx/nx;
dy=Ly/ny;
x=Lx/nx:Lx/nx:Lx;
x=(x-(Lx/2))/1000;
y=Ly/ny:Ly/ny:Ly;
y=(y-(Ly/2))/1000;
%%
for k=1:length(timesteps)
    r_x=squeeze(r(120:end,120,k));
    integ_x=trapz(r_x)*dx;
    r_y=squeeze(r(120,120:end,k));
    integ_y=trapz(r_y)*dy;
    ratio(k)=integ_x/integ_y;
    integ(k)=sqrt((integ_x.^2)+(integ_y.^2));
end
%%
t=figure;
plot(double(timesteps(1:end))*0.1/60/60,integ)
xlim([3 4.5])
xlabel('Simulation Time (hrs)','Interpreter','latex','FontSize',12)
ylim([1000 1500])
ylabel('$L$ (m)','Interpreter','latex','FontSize',12)
exportgraphics(t,[filepath,'figures/integral_length_scale.pdf'],'Resolution',300);