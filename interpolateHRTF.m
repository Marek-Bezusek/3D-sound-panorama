function HRTF = interpolateHRTF(phi, theta, phi0, theta0,NFFT, hrir_l, hrir_r)
%  This function calculates the HRTF for the specified direction or directions. 
%  If HRIR direction is not included in the database, it performs linear interpolation (1D or 2D).
% 
%  phi - sound source positions from CIPIC: azimuth in degrees
%  theta - sound source positions from CIPIC: elevation in degrees
%  phi0 - sound source azimuth: value (static) or vector (moving)
%  theta0 - sound source elevation: value (static) or vector (moving)
%  NFFT - FFT length
%  hrir_l - HRIR for left ear (format: CIPIC databse)
%  hrir_r - HRIR for right ear (format: CIPIC databse)
%
% Returns:
%  HRTF - 3D matrix of interpolated HRTF. HRTF for right and left ear in columns (HRTF_l,HRTF_r).   
%       - the third dimension contains the HRTF pair for the k-th direction (phi0(k), theta0(k))
%  

%% Calculates HRTF for the specified direction
HRTF = zeros(NFFT,2,length(phi0));  % memory allocation
for k = 1:length(phi0)             
  
%% Check if interpolation is needed
[m, m_ind]=min(abs(phi-phi0(k)));        % find index of closest azimuth
[n, n_ind]=min(abs(theta-theta0(k)));    % find index of closest elevation

if m ~= 0 && n ~= 0 % check interpolation type (none, 2D, 1D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bilinear interpolation
% Interpolation point: H0(phi0,theta0)
% 
% azimuth axis ->
%             H11(phi_1,theta_1)----------------H21(phi_2,theta_1)
%              |                   |              |
%              |------------------H0(phi0,theta0) |
%              |                   |              |
%              |                   |              |
%              |                   |              |
%             H12(phi_1,theta_2)----------------H22(phi_2,theta_2)
% 
% Interpolation of azimuth:
%  Interpolate H11 and H21 on elevation theta_1:
%   H1(phi0,theta_1) = (phi_2-phi_0)/(phi_2-phi_1)*H11 + (phi_0-phi_1)/(phi_2-phi_1)*H21
%  Interpolate H12 and H22 on elevation theta_2:
%   H2(phi0,theta_2) = (phi_2-phi_0)/(phi_2-phi_1)*H12 + (phi_0-phi_1)/(phi_2-phi_1)*H22
% 
% Interpolation of elevation (H1 and H2):
%   H0(phi0,theta_0) = (theta_2-theta_0)/(theta_2-theta_1)*H1 + (theta_0-theta_1)/(theta_2-theta_1)*H2

    %% Finding the surrounding azimuth and elevation indices between which we will interpolate
    m_phi1 = find(phi <= phi0(k),1,'last'); % lower azimuth index
    m_phi2 = find(phi >= phi0(k),1,'first'); % upper azimuth index

    n_theta1 = find(theta <= theta0(k),1,'last'); % lower elevation index
    n_theta2 = find(theta >= theta0(k),1,'first'); % upper elevation index 

    %% Weighing coefficients
    % Azimuth weighing coefficients
    w_phi1 = (phi(m_phi2)-phi0(k))/(phi(m_phi2)-phi(m_phi1)); % weighing coefficient of lower azimuth (phi_1)
    w_phi2 = (phi0(k)-phi(m_phi1))/(phi(m_phi2)-phi(m_phi1)); % weighing coefficient of upper azimuth (phi_2)

    % Elevation weighing coefficients
    w_theta1 = (theta(n_theta2)-theta0(k))/(theta(n_theta2)-theta(n_theta1)); % weighing coefficient of lower elevation (theta_1)
    w_theta2 = (theta0(k)-theta(n_theta1))/(theta(n_theta2)-theta(n_theta1)); % weighing coefficient of upper elevation (theta_2)

    %% Calculate surrounding HRTF and interpolate
    H0 = zeros(NFFT,2); % memory allocation for interpolated HRTF pair
    L = size(hrir_l,3); % impulse response length
    window = hamming(L,'periodic');
    
    hrir = hrir_l;
    for i = 1:2        % loop for HRTF_l and HRTF_r calculation
    H = zeros(NFFT,4); % memory allocation for surrounding HRTF (impulse response length + zero pad)
    H(1:L,1) = squeeze(hrir(m_phi1,n_theta1,:)); % H11
    H(1:L,2) = squeeze(hrir(m_phi2,n_theta1,:)); % H21
    H(1:L,3) = squeeze(hrir(m_phi1,n_theta2,:)); % H12
    H(1:L,4) = squeeze(hrir(m_phi2,n_theta2,:)); % H22
    
    H(1:L,1:4) = H(1:L,1:4) .* window; % weight impulse response by window
    H = fft(H,NFFT,1); % surrounding HRTFs in columns

    % Interpolate azimuth
    H1 = w_phi1*H(:,1) + w_phi2*H(:,2); % H1(phi0,theta_1)
    H2 = w_phi1*H(:,3) + w_phi2*H(:,4); % H2(phi0,theta_2)

    % Interpolate elevation
    H0(:,i) =  w_theta1*H1 + w_theta2*H2;

    if i==2; break; end
    hrir = hrir_r;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D linear interpolation
% Interpolate only azimuth or elevation
% Interpolation point: H0(alpha0,beta0)
% 
% interpolation axis (alpha) ->
%             H1(alpha_1,beta0)-------H0---------H2(alpha_2,beta0)
% 
elseif m ~= 0 || n ~= 0     % check interpolation type (none, 2D, 1D)
    if m ~= 0               % variables for azimuth interpolation
        alpha = phi;
        alpha0 = phi0(k);
        n_beta0 = n_ind;    % elevation index 
    elseif n ~= 0           % variables for elevation interpolation
        alpha = theta;
        alpha0 = theta0(k);
        n_beta0 = m_ind;    % azimuth index
    end
    %% Finding the surrounding indices between which we will interpolate
    m_alpha1 = find(alpha <= alpha0,1,'last'); % lower index
    m_alpha2 = find(alpha >= alpha0,1,'first'); % upper index

    % Weighing coefficients
    w_alpha1 = (alpha(m_alpha2)-alpha0)/(alpha(m_alpha2)-alpha(m_alpha1)); % lower index weighing coefficient (alpha_1)
    w_alpha2 = (alpha0-alpha(m_alpha1))/(alpha(m_alpha2)-alpha(m_alpha1)); % upper index weighing coefficient (alpha_2)

    %% Calculate surrounding HRTF and interpolate
    H0 = zeros(NFFT,2); % memory allocation for interpolated HRTF pair
    L = size(hrir_l,3); % impulse response length
    window = hamming(L,'periodic');
    
    
    hrir = hrir_l;
    for i = 1:2        % loop for HRTF_l and HRTF_r calculation
    H = zeros(NFFT,2); % memory allocation for surrounding HRTF (impulse response length + zero pad)
    H(1:L,1) = squeeze(hrir(m_alpha1,n_beta0,:)); % H1 (alpha1,beta0)
    H(1:L,2) = squeeze(hrir(m_alpha2,n_beta0,:)); % H2 (alpha2,beta0)
    
    H(1:L,1:2) = H(1:L,1:2) .* window; % weight impulse response by window
    H = fft(H,NFFT,1); % surrounding HRTFs in columns

    % Interpolation (axis alpha)
    H0(:,i) = w_alpha1*H(:,1) + w_alpha2*H(:,2); % H0(alpha0,beta0)


    if i==2; break; end
    hrir = hrir_r;
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolation is not needed
% Calculate HRTF from HRIR included in CIPIC database
else
    H0 = zeros(NFFT,2);
    L = size(hrir_l,3); % impulse response lenght
    window = hamming(L,'periodic');
    
    H0(1:L,1) = squeeze(hrir_l(m_ind,n_ind,:)); 
    H0(1:L,2) = squeeze(hrir_r(m_ind,n_ind,:));
    
    H0(1:L,1:2) = H0(1:L,1:2) .* window; % weight impulse response by window
    H0 = fft(H0,NFFT,1); % HRTF in columns (HRFT_l, HRTF_r)
  
end % IF END - interpolation type (none, 2D, 1D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HRTF(:,:,k)=H0; % HRTF matrix (3D)
end % FOR LOOP END - calculation of HRTF for specified directions

end % FUNCTION END