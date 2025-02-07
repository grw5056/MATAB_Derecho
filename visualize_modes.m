%%
% Parameters
x_wavelength =1500;  % Wavelength in the x direction (can be adjusted)
y_wavelength = 12000;  % Wavelength in the y direction (can be adjusted)
amplitude = 0.1;     % Amplitude of the wave (can be adjusted)
x_range = 12000;      % The range in the x direction
y_range = 12000;      % The range in the y direction
nx=240;
ny=240;

% Create a grid of x and y values
x = linspace(50, x_range, 240);  % x values from 0 to x_range, with 500 points
y = linspace(50, y_range, 240);  % y values from 0 to y_range, with 500 points
[X, Y] = meshgrid(x, y);        % Create 2D grid from x and y vectors

% Define the wave function: A*sin(kx*X + ky*Y)
k_x = 2*pi / x_wavelength;  % Wave number in the x direction
k_y = 2*pi / y_wavelength;  % Wave number in the y direction
for i=1:nx
    for j=1:ny
        Z(i,j) = amplitude * sin(k_x * x(i) + k_y * y(j));  % Sine wave function
    end
end
% Plot the 2D periodic wave
figure;
contourf(X,Y,Z)
