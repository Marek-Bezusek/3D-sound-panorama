function [y, n_segments] = segmentation(x, winlen, winover)

% x         - input signal
% winlen    - window lenght or window function
% winover   - window overlap
% y         - segmented signal (segments in columns)
% n_segments - number of segments
%% Check the window length and overlap
wl = 0;
winover = floor(winover);

if(length(winlen) > 1)
    wl = length(winlen);
else
    wl = winlen;
end

%% Get number of segments
n_segments = ceil((length(x)-winover)/(wl-winover));     

%% Pad last segment with zeros if necessary
if(mod(length(x),wl) ~= 0)                          
    x(end+1:n_segments*wl) = 0;                         
end

%% Prepare matrix
y = zeros(wl, n_segments);                              

%% Indexes for segmentation process
idx_samples = (1:wl).';                                 % indexes of samples within a segment (segments in columns)
idx_step = 0:(wl-winover):(n_segments-1)*(wl-winover);  % indexes of segments within signal (segments in columns)
%% Segmentation
y(:) = x(idx_samples(:,ones(1,n_segments)) + idx_step(ones(1,wl),:)); 
%% Apply window
if(length(winlen) > 1)
    y = y.*winlen(:,ones(1,n_segments));                 
end