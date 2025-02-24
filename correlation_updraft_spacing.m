clc;
clear;
close all;
%%
filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
filename_r='w_slice_output/r_z16_t1000_k120.nc';
filename_w='w_slice_output/w_z16_t1000_k120.nc';
ug=15;
Lx=12000;
Ly=12000;
nx=240;
ny=240;
avg_width=100;
lvl=16;
filter=120;
contour_levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
use_levels=0.2;
ABL_depth=1100;
updraft_width=1.5; %Excepted updraft width normalized by zi
%% Colormap
colors = [
    0.8500, 0.3250, 0.0980;  % Red
    0.9290, 0.6940, 0.1250;  % Yellow (Gold)
    0.1960, 0.8030, 0.1960;  % Bright Green (Lighter Green)
    0.0940, 0.4110, 0.1030;  % Dark Green (Darker Green)
    0.3010, 0.7450, 0.9330;  % Light Blue
    0.6350, 0.0780, 0.1840;  % Dark Red
    0.4940, 0.1840, 0.5560;  % Purple
    0.9290, 0.1840, 0.6230;  % Magenta
    0.0780, 0.4780, 0.6980;  % Teal
];
%%
xtick_input_r=linspace(-(Lx/1000),(Lx/1000),7);
xticklabel_input_r=compose('%d',xtick_input_r);
xtick_input_r(1)=(-Lx)+(Lx/nx/1000);
ytick_input_r=linspace(-(Ly/1000),(Ly/1000),7);
yticklabel_input_r=compose('%d',ytick_input_r);
ytick_input_r(1)=(-Ly)+(Ly/ny/1000);
%%
r=ncread([filepath,filename_r],'r');
w=ncread([filepath,filename_w],'w');
timesteps=ncread([filepath,filename_r],'time');
%%
padding=Lx/2;
dx=Lx/nx;
dy=Ly/ny;
r=padarray(r,[padding/dx,padding/dy],'circular');
x=Lx/nx:Lx/nx:Lx;
x=(x-(Lx/2))/1000;
x=cat(2,x((end-(padding/dx))+1:end)-Lx/1000,x,x(1:(padding/dx))+Lx/1000);
y=Ly/ny:Ly/ny:Ly;
y=(y-(Ly/2))/1000;
y=cat(2,y((end-(padding/dy))+1:end)-Ly/1000,y,y(1:(padding/dx))+Ly/1000);
[X,Y]=meshgrid(x,y);
%%
corner_coordinates = [
    -(Lx/2000), -(Ly/2000);   % First corner (x1, y1)
    -(Lx/2000), (Ly/2000);   % Second corner (x2, y2)
    (Lx/2000), (Ly/2000);   % Third corner (x3, y3)
    (Lx/2000), -(Ly/2000)   % Fourth corner (x4, y4)
    -(Lx/2000), -(Ly/2000);
];
%%
for k=1:length(timesteps)
    %%
    t=figure('Visible','off');
    [~,hContourLines]=contour(x,y,squeeze(r(:,:,k))',contour_levels,'--','HandleVisibility','off');
    hold on
    plot(NaN,NaN,'color',colors(1,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(2,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(3,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(4,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(5,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(6,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(7,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(8,:),'linestyle','--')
    hold on
    plot(NaN,NaN,'color',colors(9,:),'linestyle','--')
    hold on
    plot(corner_coordinates(:,1), corner_coordinates(:,2), '-k', 'LineWidth', 2);
    colormap(colors);
    pbaspect([1 1 1])
    xlabel('$x$-shift (km)','Interpreter','latex','FontSize',12)
    ylabel('$y$-shift (km)','Interpreter','latex','FontSize',12)
    xticks(xtick_input_r)
    xticklabels(xticklabel_input_r)
    yticks(ytick_input_r)
    yticklabels(yticklabel_input_r)
    legend({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'},'location','eastoutside')
    text(0.01,1.02,['Timestep = ',num2str(timesteps(k))],'units','normalized','Interpreter','latex')
    set(gca,'Fontsize',12)
    contourData = hContourLines.ContourMatrix;
    cont_num=1;
    i=1;
    phi_num=1;
    temp_phi=0;
    while i <= size(contourData, 2)
        level = contourData(1, i); % The first element is the contour level
        numPoints = contourData(2, i); % Find the length of this contour line
        xCoords = contourData(1,i+1:i+numPoints);
        yCoords = contourData(2,i+1:i+numPoints);
        if xCoords(1)~=xCoords(end) || yCoords(1)~=yCoords(end) 
            i = i + numPoints + 1; % Move to the next contour line
        elseif length(xCoords)<=15
            i=i+numPoints+1;
        elseif all(round(xCoords,2)>-(Lx/2000)) && all(round(xCoords,2)<=(Lx/2000)) && all(round(yCoords,2)>-(Lx/2000)) && all(round(yCoords,2)<=(Lx/2000))
            if any(level==use_levels)
                [struc,temp_phi(phi_num)]=fit_ellipse(xCoords,yCoords);
                X0(cont_num)=struc.X0_in;
                Y0(cont_num)=struc.Y0_in;
                phi_num=phi_num+1;
                i = i + numPoints + 1; % Move to the next contour line
                cont_num=cont_num+1;
            else
                [struc]=fit_ellipse(xCoords,yCoords);
                X0(cont_num)=struc.X0_in;
                Y0(cont_num)=struc.Y0_in;
                i = i + numPoints + 1; % Move to the next contour line
                cont_num=cont_num+1;
            end
        elseif any(round(xCoords,2)>-(Lx/2000)) && any(round(xCoords,2)<=(Lx/2000)) && any(round(yCoords,2)>-(Lx/2000)) && any(round(yCoords,2)<=(Lx/2000))
            [struc]=fit_ellipse(xCoords,yCoords);
            if round(struc.X0_in,2)>-(Lx/2000) && round(struc.X0_in,2)<=(Lx/2000) && round(struc.Y0_in,2)>-(Lx/2000) && round(struc.Y0_in,2)<=(Lx/2000)
                X0(cont_num)=struc.X0_in;
                Y0(cont_num)=struc.Y0_in;
                i = i + numPoints + 1; % Move to the next contour line
                cont_num=cont_num+1;
            else
                i=i+numPoints+1;
            end
        else
            i=i+numPoints+1;
        end
    end
    HCR_angle(k)=HCR_orientation(squeeze(w(:,:,k)));
    phi(k)=mean(temp_phi);
    rotated_X0=(X0.*cosd(-HCR_angle(k)))-(Y0.*sind(-HCR_angle(k)));
    rotated_Y0=(X0.*sind(-HCR_angle(k)))+(Y0.*cosd(-HCR_angle(k)));
    scatter(X0,Y0,'filled','MarkerFaceColor','black','MarkerEdgeColor','black','HandleVisibility','off');
    hold on
    scatter(rotated_X0,rotated_Y0,'filled','MarkerFaceColor','red','MarkerEdgeColor','red','HandleVisibility','off');
    exportgraphics(t,[filepath,'figures/r_snapshots/r_z',num2str(lvl),'_t',num2str(avg_width),'_k',num2str(filter),'_',num2str(timesteps(k)),'.png'],'Resolution',300);
    %
    for n=1:length(rotated_Y0)
        if rotated_Y0(n)>(Lx/2000)
            rotated_Y0(n)=rotated_Y0(n)-(Lx/1000);
        elseif rotated_Y0(n)<=-(Lx/2000)
            rotated_Y0(n)=rotated_Y0(n)+(Lx/1000);
        end
    end
    group_averages=groupaverages(Y0,updraft_width*ABL_depth);
    if max(group_averages)>((Lx/2000)-1) && min(group_averages)<-((Lx/2000)-1)
        [minval,minloc]=min(group_averages);
        [maxval,maxloc]=max(group_averages);
        minval=minval+(Lx/1000);
        val=(minval+maxval)/2;
        if val>6
            val=val-(Lx/1000);
        end
        group_averages([minloc,maxloc])=[];
        group_averages=[group_averages,val];
        group_averages=sort(group_averages);
    end
    if length(group_averages)>1
        for n=1:length(group_averages)-1
            group_averages=sort(group_averages);
            diff(n)=group_averages(n+1)-group_averages(n);
        end
        spacing(k)=mean(diff);
    else
        spacing(k)=Ly;
    end
    group_averages_rot=groupaverages(rotated_Y0,1.5*ABL_depth);
    if max(group_averages_rot)>((Lx/2000)-1) && min(group_averages_rot)<-((Lx/2000)-1)
        [minval_rot,minloc_rot]=min(group_averages_rot);
        [maxval_rot,maxloc_rot]=max(group_averages_rot);
        minval_rot=minval_rot+(Lx/1000);
        val_rot=(minval_rot+maxval_rot)/2;
        if val_rot>6
            val_rot=val_rot-(Lx/1000);
        end
        group_averages_rot([minloc_rot,maxloc_rot])=[];
        group_averages_rot=[group_averages_rot,val_rot];
        group_averages_rot=sort(group_averages_rot);
    end
    if length(group_averages_rot)>1
        for n=1:length(group_averages_rot)-1
            group_averages_rot=sort(group_averages_rot);
            diff_rot(n)=group_averages_rot(n+1)-group_averages_rot(n);
        end
        spacing_rot(k)=mean(diff_rot);
    else
        spacing_rot(k)=Ly;
    end  
end
%%
save([filepath,'figures/r_snapshots/correlation_spacing_z',num2str(lvl),'_t',num2str(avg_width),'_k',num2str(filter),'.txt'],"spacing","-ascii")
save([filepath,'figures/r_snapshots/correlation_spacing_rot_z',num2str(lvl),'_t',num2str(avg_width),'_k',num2str(filter),'.txt'],"spacing_rot","-ascii")