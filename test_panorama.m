%% Script for testing panoramaHRTF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;
%% Import HRIR from CIPIC database
load('CIPIC/subject_003/hrir_final','hrir_l', 'hrir_r');

%% Import input audio signal
[x, Fs] = audioread('audio/Piano FF C 3 off.wav');

%% Calculate input signal: White noise
Ts = 1/Fs; % sampling period
tn = 10; % signal duration in seconds
x = randn(tn*Fs,1); % generate White noise of duration tn

%% Calculate input signal: sine
Fs = 44100; % sampling frequency
t = 10; % signal duration
N = Fs*t; % number of samples
n = 0:N-1;
f = 333;
x = 0.5*sin(2*pi*f*n/Fs);

%% Calculate input signal: sine sweep
Fs = 1000;
t = 0:1/Fs:10;
x = chirp(t,0,20,1000,'linear', 45);

%% Sound source positions included in CIPIC database:
% phi = -80 -65 -55 -45:5:45 55 65 80  % azimuth
% theta = -45:5.625:230.625            % elevation
%% Test: moving source in one axis
% azimuth
phi_min = -80; % degrees
phi_max = 80; % degrees
step = 2; % degrees

phi0 = phi_min:step:phi_max;
theta0 = zeros(1,length(phi0));

% elevation
% theta_min = -40;
% theta_max = 80;
% step= 2;
% 
% theta0 = theta_min:step:theta_max;
% phi0 = zeros(1,length(theta0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Panorama 
y = panoramaHRTF(x, phi0, theta0, hrir_l, hrir_r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Play output
sound(y, Fs);
