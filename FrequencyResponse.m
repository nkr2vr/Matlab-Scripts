%% Created by Noah Rohrlich
%% June 2017
%% Designed to accompany Digital Signals Processing lab for instruction
%%  of ECE3430 at University of Virginia, led by Prof. Todd DeLong.
%% This program plots the theoretical frequency magnitude response for
%%  a 64-sample averaging filter and notch filter at Fs = 3840Hz. Also plots the 
%%  experimentally-collected data points (if any were entered in the vectors).

%% Define ranges and variables of importance
clear;clc;
figure_num = 1;
%Prepare the legend
legend_strings = {'Analytical'};

% These vectors are used to store the experimental data if available.  The
% first vector contains the frequencies for which experimental data was
% collected, and the second vector contains the magnitude of the output.
Exp_f_avg = [];
Exp_mag_avg = [];

Exp_f_notch = [];
Exp_mag_notch = [];

% Define the analytical frequency range
f = linspace(0,500,1000); %Hz
Fs = 3840;

% Define complex variable z
% Note: z = exp(j * 2*pi*f / Fs), a complex variable used often in DSP
z = cos(2*pi*f./Fs)+1j*sin(2*pi*f./Fs);
b = ones(length(f), 1)';    %Initialize numerator of transfer function

%% Calculate theoretical response

% Averaging Filter LCCDE:
% y[n] = 1/64 * (x[n] + x[n-1] + ... + x[n-63]) =>
%
% Transfer Function:
%        Y(z)   1   1 + z^-1 + z^-2 + ... + z^-63
% H(z) = --- = --- ------------------------------
%        X(z)   64                1
%
num_points = 64; %Number of samples being used in average
for i = 0:1:(num_points-1)
    b = b + z.^(-i);
end
a = num_points;      % Denominator of transfer function
h_avg = b./a;


% Notch Filter LCCDE:
% y[n] = 1/2 * (x[n] + x[n-d]) =>
%
% Transfer Function:
%        Y(z)  1 + z^(-d)
% H(z) = --- = ----------
%        X(z)      2

f_cutoff = 60; %Hz
d = Fs / (2*f_cutoff); % Magic index used to eliminate odd harmonics of a certain frequency
b = 1 + z.^(-d);
a = 2;      % Denominator of transfer function
h_notch = b./a;


%% Display results

%Averaging filter
figure(figure_num); figure_num = figure_num + 1;
plot(f,20*log10(abs(h_avg)), '-k');
grid on
hold on

% Read experimental data if its available
if (isempty(Exp_mag_avg) ~= 1)
    plot(Exp_f_avg,20*log10(Exp_mag_avg),'.b');
    legend_strings{2} = 'Experimental';
end

hold off
title('Frequency Response of Averaging Filter')
xlabel('Frequency f (Hz)')
ylabel('|H(z)| (dB)')
legend(legend_strings);
grid on

%Notch Filter
figure(figure_num); figure_num = figure_num + 1;
plot(f,20*log10(abs(h_notch)), '-k');
grid on
hold on

% Read experimental data if its available
if (isempty(Exp_mag_notch) ~= 1)
    plot(Exp_f_notch,20*log10(Exp_mag_notch),'.b');
    legend_strings{2} = 'Experimental';
elif (size(legend_strings) > 1)
    legend_strings(2) = []; %delete extra entry in legend
end

hold off
title('Frequency Response of Notch Filter')
xlabel('Frequency f (Hz)')
ylabel('|H(z)| (dB)')
legend(legend_strings);
grid on
