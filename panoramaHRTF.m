function y = panoramaHRTF(x, phi0, theta0, hrir_l, hrir_r)
% This function pans moving sound source.
%  x - input audio signal (mono)
%  phi0 - vector of sound source positions: azimuth in degrees
%  theta0 - vector of sound source positions: elevation in degrees
%  hrir_l - HRIR for left ear (CIPIC databse)
%  hrir_r - HRIR for right ear (CIPIC databse)

%% Check input signal
% Transpose to column vector
if size(x,1)<size(x,2), x = x(1,:)'; end    
% Stereo to mono
if size(x,2) == 2 
    x=(x(:,1) + x(:,2))/2;
    disp('The input signal was mixed down to mono');
elseif size(x,2) > 2
    error('The input signal must be mono or stereo');
end

%% Check range of positions
if max(phi0) > 80 || min(phi0) < -80
    error('The azimuth angle must be between -80 ° and 80 °');
end
if max(theta0) > 230.625 || min(theta0) < -45
    error('The elevation angle must be in the range -45 ° to 230.625 °');
end
if length(phi0) ~= length(theta0)
    error('The number of entered azimuth angles does not match the number of elevation angles');
end

%% HRIR database variables
L = size(hrir_l,3); % impulse response lenght
% Sound sources positions included in CIPIC
phi = ([-80 -65 -55 -45:5:45 55 65 80]); % azimuth
theta = (-45:5.625:230.625); % elevation

%% Variables for segmentation and for FFT
Nx = length(x);             % number of samples of input signal
winlen = L;                 % window length
winover = winlen/2;         % window overlap
window = hamming(winlen,'periodic');
winshift = winlen-winover;  % window step
% [tf,~,maxDeviation] =  iscola(window,winover,'ola');
NFFT = 2^(ceil(log2(winlen+L-1))); % FFT lenght

%% Input signal segmentation
[xSeg, Nseg] = segmentation(x, window, winover); % xSeg(segments in columns) Nseg(number of segments)
clear x;
xSeg = [xSeg; zeros(NFFT-winlen, Nseg)]; % zero pad

%% Compute HRTF from HRIR for specified sound source positions
HRTF = interpolateHRTF(phi, theta, phi0, theta0, NFFT, hrir_l, hrir_r);

%% Calculate speed of moving source: period of direction change (change of HRTF)
T = (Nseg+1)/length(phi0);
n = 1:Nseg;
ind = floor(n*1/T)+1;

%% Fast convolution of input signal with HRIR. Overlap add
y = zeros(Nx + NFFT,2); % output audio signal (stereo)

for i = 1:Nseg
        % convolution
        ySeg = real(ifft(fft(xSeg(:,i),NFFT) .* HRTF(:,:,ind(i)), NFFT));
        % overlap add
        outindex = ((i-1)*winshift+1):((i-1)*winshift+NFFT);
        y(outindex,:) = y(outindex,:) + ySeg; 
end


end