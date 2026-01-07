% =========================================================================
% Project: Free Space Optical Propagation (Angular Spectrum Method)
% Author:  sun-optica
% Date: 2026-01-07
% Description: 
%   This script simulates the diffraction of a plane wave through a 
%   square aperture using the Angular Spectrum Method (ASM).
% =========================================================================

clc; clear; close all;

%% 1. System Parameters (Physical Units)



wavelength = 632.8e-9;    % Wavelength (He-Ne Laser), unit: m
k = 2 * pi / wavelength;  % Wavenumber

L = 10e-3;                % Side length of the simulation window (10 mm)
N = 1024;                 % Number of sampling points (Power of 2 for FFT efficiency)
dx = L / N;               % Sampling interval
x = (-N/2 : N/2-1) * dx;  % Coordinate axis
[X, Y] = meshgrid(x, x);  % 2D Grid

z = 0.5;                  % Propagation distance (50 cm)

%% 2. Define Initial Field (Source Plane z=0)

% Creating a Square Aperture (2mm x 2mm)

aperture_width = 2e-3;
U0 = double(abs(X) < aperture_width/2 & abs(Y) < aperture_width/2);

% Visualization of Source


figure(1);
subplot(1, 2, 1);
imagesc(x*1e3, x*1e3, abs(U0).^2);
colormap('hot'); axis square; colorbar;
title('Source Intensity (z = 0)');
xlabel('x (mm)'); ylabel('y (mm)');

%% 3. Angular Spectrum Propagation (The Core Algorithm)
% Step 3.1: Fourier Transform of the initial field


U0_spectrum = fftshift(fft2(U0));

% Step 3.2: Define Frequency Coordinates

df = 1 / L;               % Frequency resolution
fx = (-N/2 : N/2-1) * df;
[FX, FY] = meshgrid(fx, fx);

% Step 3.3: Define Transfer Function H(fx, fy)
% H = exp(j * k * z * sqrt(1 - (lambda*fx)^2 - (lambda*fy)^2))

term_sqrt = 1 - (wavelength * FX).^2 - (wavelength * FY).^2;
term_sqrt(term_sqrt < 0) = 0; % Remove evanescent waves
H = exp(1j * k * z * sqrt(term_sqrt));

% Step 3.4: Propagation in Frequency Domain
Uz_spectrum = U0_spectrum .* H;

% Step 3.5: Inverse Fourier Transform to get field at z
Uz = ifft2(ifftshift(Uz_spectrum));

%% 4. Result Visualization
Intensity = abs(Uz).^2;

subplot(1, 2, 2);
imagesc(x*1e3, x*1e3, Intensity);
colormap('hot'); axis square; colorbar;
title(['Diffraction Pattern at z = ', num2str(z*100), ' cm']);
xlabel('x (mm)'); ylabel('y (mm)');

% Save the result for GitHub
% standardizing the output looks professional
saveas(gcf, '../results/diffraction_pattern.png');
disp('Simulation finished. Result saved to ../results/diffraction_pattern.png');
