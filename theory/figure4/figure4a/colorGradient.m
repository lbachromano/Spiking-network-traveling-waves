function cmap = colorGradient(color1, color2, N)
    % color1, color2: 1x3 vectors (RGB), e.g. [1,0.5,0] for orange
    % N: number of distinct colors you want
    if nargin < 3
        N = 256;  % default to 256 if not specified
    end
    
    cmap = zeros(N, 3);
    for i = 1:N
        t = (i-1) / (N-1); 
        cmap(i, :) = (1 - t)*color1 + t*color2; 
    end
end
