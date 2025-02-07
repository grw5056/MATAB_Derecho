clear;
clc;
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
tini=108000;
tend=160000;
Lz=4000;
dz=25;
ug=15;
T_scale=300;
wt_s=0.05;
g=9.81;
k=0.4;
u_heights=dz/2:dz:Lz+(dz/2);
w_heights=0:dz:Lz;
%%
theta=load([filepath,'theta_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
dTdz=load([filepath,'dTdz_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
u=load([filepath,'u_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
v=load([filepath,'v_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
w=load([filepath,'w_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
uw=load([filepath,'uw_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
txz=load([filepath,'txz_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
vw=load([filepath,'vw_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
tyz=load([filepath,'tyz_z_0',num2str(tini),'_0',num2str(tend),'.txt']);
%% ABL Depth
[max_val,max_ind]=max(dTdz);
zi=w_heights(max_ind);
%% friction velocity
u_temp=(u(2:end)+u(1:end-1))/2;
u_temp=cat(1,0,u_temp);
v_temp=(v(2:end)+v(1:end-1))/2;
v_temp=cat(1,0,v_temp);
uw_flux=uw-(u_temp.*w)+txz;
vw_flux=vw-(v_temp.*w)+tyz;
friction_velocity=(((uw_flux(1).^2)+(vw_flux(1).^2)).^(1/4))*ug;
%% Average ABL temperature
dTdz_temp=(dTdz(1:end-1)+dTdz(2:end))/2;
[max_val_temp,max_ind_temp]=max(dTdz_temp);
ABL_T=mean(theta(1:max_ind_temp))*T_scale;
%% Obukhov length
L=((friction_velocity.^3)*ABL_T)./(g*k*wt_s);
%% zi/L
stability=zi/L;