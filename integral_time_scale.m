clc;
clear;
close all;
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
%filepath='/glade/derecho/scratch/gwarner/het_ug8_lx6_phi90/output/';
filename='w_slice_output/w_z16_t0001_k120.nc';
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
w_full=ncread([filepath,filename],'w');
timesteps=ncread([filepath,filename],'time');
%%
for n=1:length(timesteps)-1
    w_init=squeeze(w_full(:,:,n));
    w_init=w_init-mean(mean(w_init));
    clear corr
    for m=1:length(timesteps)-n
        w_shift=squeeze(w_full(:,:,n+m));
        w_shift=w_shift-mean(mean(w_shift));
        R=corrcoef(w_init(:),w_shift(:));
        corr(m)=R(1,2);
    end
    integ_time_scale(n)=trapz(corr)*10;
end
%%
t=figure;
plot(double(timesteps(1:end-1))*0.1/60/60,integ_time_scale)
xlim([3 4.5])
xlabel('Simulation Time (hrs)','Interpreter','latex','FontSize',12)
ylim([0 100])
ylabel('$T$ (s)','Interpreter','latex','FontSize',12)
exportgraphics(t,[filepath,'figures/integral_time_scale.pdf'],'Resolution',300);