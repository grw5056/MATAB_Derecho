clear
clc
close all
%%
%filepath='/glade/work/gwarner/HCR/hom_ug15_phi0_H005_large/output/';
filepath='/glade/derecho/scratch/gwarner/het_ug8_lx6_phi90/output/';
%filename='w_slice_output/w_z16_t1000_k120.nc';
filename='w_z32_snap_k006.nc';
dx=25;
dy=25;
ug=8;
Lx=6000;
Ly=6000;
lvl=32;
avg_width=1;
filter=6;
%%
x=dx:dx:Lx;
y=dy:dy:Ly;
[X,Y]=meshgrid(x,y);
%%
w=ncread([filepath,filename],'w')*ug;
timesteps=ncread([filepath,filename],'time');
time=double(timesteps)*.1/60/60;
%%
for n=1:length(timesteps)
    w_slice=w(:,:,n);
    W=fft2(w_slice);
    W_shifted = fftshift(W); % Shift zero frequency to the center
    magnitude = abs(W_shifted/(240*240));
    %
    [nx,ny] = size(w_slice);
    fx = (0:nx-1)/nx; % Frequency in x direction
    fy = (0:ny-1)/ny; % Frequency in x direction
    %
    fx = fx - (1/2); % Shift to center
    fy = fy - (1/2); % Shift to center
    %
    for m=1:3
        [peak_x,peak_y] = find(max(magnitude(:))==magnitude);
        %
        index_x(n,m)=peak_x(2);
        index_y(n,m)=peak_y(2);
        primary_frequencies_x = fx(peak_x);
        primary_frequencies_y = fy(peak_y);
        %
        wavelength_x(n,m) = (1 ./ primary_frequencies_x(2))*dx;
        wavelength_y(n,m) = (1 ./ primary_frequencies_y(2))*dy;

        amp(n,m)=max(magnitude(:));
        magnitude(peak_x(1),peak_y(1))=0;
        magnitude(peak_x(2),peak_y(2))=0;
    end
end
%%
indexs=cat(3,index_x,index_y);
% Initialize a cell array to store unique pairs
unique_pairs = zeros(1,2);
% Loop through all elements in the 1st and 2nd dimensions
for a = 1:size(indexs, 1)  % Loop over 100
    for b = 1:size(indexs, 2)  % Loop over 5
        % Extract the pair in the third dimension (size 2)
        pair = indexs(a, b, :);
        % Reshape pair to a row vector
        pair = reshape(pair, 1, 2);
        % Add the pair to the list of unique pairs (if it's not already there)
        if ~ismember(pair, unique_pairs, 'rows')
            unique_pairs = [unique_pairs; pair];
        end
    end
end
unique_pairs=unique_pairs(2:end,:);
%%
amp_sort=zeros(size(unique_pairs,1),size(w,3));
for n=1:size(w,3)
    w_slice=w(:,:,n);
    W=fft2(w_slice);
    W_shifted = fftshift(W); % Shift zero frequency to the center
    magnitude = abs(W_shifted/(240*240));
    for i=1:size(unique_pairs,1)
        amp_sort(i,n)=magnitude(unique_pairs(i,1),unique_pairs(i,2));
    end
end
%%
if size(unique_pairs,1)>6
    end_point=6;
else
    end_point=size(unique_pairs,1);
end
t=figure;
for n=1:end_point
    plot(time,amp_sort(n,:))
    if ((1 ./ fx(unique_pairs(n,1))*dx)/1000)==inf && ((1 ./ fy(unique_pairs(n,2))*dy)/1000)~=inf
        formattedValue{n} = ['$\lambda_{x}$ = $\infty$, $\lambda_{y}$ = ',num2str(((1 ./ fy(unique_pairs(n,2))*dy)/1000), '%.1f'),' km'];
    elseif ((1 ./ fx(unique_pairs(n,1))*dx)/1000)~=inf && ((1 ./ fy(unique_pairs(n,2))*dy)/1000)==inf
        formattedValue{n} = ['$\lambda_{x}$ = ',num2str(((1 ./ fx(unique_pairs(n,1))*dx)/1000), '%.1f'),' km, $\lambda_{y}$ = $\infty$'];
    elseif ((1 ./ fx(unique_pairs(n,1))*dx)/1000)==inf && ((1 ./ fy(unique_pairs(n,2))*dy)/1000)==inf
        formattedValue{n} = '$\lambda_{x}$ = $\infty$, $\lambda_{y}$ = $\infty$';
    else
        formattedValue{n} = ['$\lambda_{x}$ = ',num2str(((1 ./ fx(unique_pairs(n,1))*dx)/1000), '%.1f'),' km, $\lambda_{y}$ = ',num2str(((1 ./ fy(unique_pairs(n,2))*dy)/1000), '%.1f'),' km'];
    end    
    hold on
end
xlabel('Time (hr)','Interpreter','latex')
xlim([time(1) time(end)])
ylabel('Mode amplitude (m s$^{-1}$)','Interpreter','latex')
legend(formattedValue,'Interpreter','latex','Location','southoutside','NumColumns',2,'FontSize',12)
set(gca,'fontsize',12)
exportgraphics(t,[filepath,'figures/w_snapshots/fourier_modes_z',num2str(lvl),'_t',num2str(avg_width),'_k',num2str(filter),'.png'],'Resolution',300);
%%
angle=atan2d(wavelength_y(:,1),wavelength_x(:,1));
for n=1:size(w,3)
    w_slice=w(:,:,n);
    w_rot = imrotate(w_slice, angle(n), 'bilinear', 'crop');
    W=fft2(w_rot);
    W_shifted = fftshift(W); % Shift zero frequency to the center
    magnitude = abs(W_shifted/(240*240));
    [peak_x,peak_y] = find(max(magnitude(:))==magnitude);
    %
    index_y=peak_y(2);
    primary_frequencies_y = fy(peak_y);
    spacing(n) = (1 ./ primary_frequencies_y(2))*dy;
end
save([filepath,'figures/w_snapshots/fourier_spacing_z',num2str(lvl),'_t',num2str(avg_width),'_k',num2str(filter),'.txt'],"spacing","-ascii")