clc;
clear;
close all;
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
%filepath='/glade/derecho/scratch/gwarner/het_ug8_lx6_phi90/output/';
filename='w_slice_output/w_z16_t0001_k004.nc';
%filename='r_z32_snap_k120.nc';
ug=15;
Lx=12000;
Ly=12000;
nx=240;
ny=240;
avg_width=1;
lvl=16;
filter=4;
neg_angles=1;
contour_levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
use_levels=[0.2,0.3,0.4,0.5,0.6];
ABL_depth=1100;
%%
w=ncread([filepath,filename],'w');
timesteps=ncread([filepath,filename],'time');
dx=Lx/nx;
dy=Ly/ny;
x=Lx/nx:Lx/nx:Lx;
x=(x-(Lx/2))/1000;
y=Ly/ny:Ly/ny:Ly;
y=(y-(Ly/2))/1000;
[X,Y]=meshgrid(x,y);
%
r = sqrt(X.^2 + Y.^2);  % Radial distance
theta = atan2d(Y, X);     % Azimuthal angle

% Step 4: Discretize the cylindrical coordinates
r_bins = linspace(0,max(max(r)),51);   % 50 radial bins
theta_bins = linspace(-180,180,51);    % 50 angular bins

% Step 5: Radial averaging (mean value within each r-bin)
radial_profile = zeros(1, length(r_bins)-1);
for i = 1:length(r_bins)-1
    for j = 1:length(theta_bins)-1
        mask = (r >= r_bins(i)) & (r < r_bins(i+1)) & (theta >= theta_bins(j)) & (theta < theta_bins(j+1));
        cyc(i,j) = mean(w(mask));
    end
end

