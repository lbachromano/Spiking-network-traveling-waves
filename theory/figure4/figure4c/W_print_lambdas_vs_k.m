
% -------------------------------------------------------------------------
% Eigenvalue stability analysis across velocity parameter  
%
% This script visualizes the maximum real part of the
% eigenvalue spectrum as a function of spatial wavenumber k, for different
% values of propagation velocity. For each velocity:
%   - Load precomputed stability data (local eigenvalue minima per k)
%   - Extract representative dominant growth rate (max Re(lambda))
%   - Interpolate spectrum over k for visualization
%   - Identify most unstable wavenumber k*
%
% Outputs:
%   (1) Dispersion-like curves Re(lambda)(k) for each velocity
%   (2) Selected instability peak k*(v)
% -------------------------------------------------------------------------


List_velocities = [200 300 366.6667 433.3333 500];

% Parameters
orange = [1, 0.8, 0];
purple = [0.5, 0, 0.5];
N = length(List_velocities);
cmap_orange2purple = colorGradient(orange, purple, N);

pt_more = 1000;

g = 4.1;
nuext = 0.35;

figure; hold on;
xlabel('k'); ylabel('Re(\lambda)');
title('Real parts of eigenvalues vs k');

k_max = nan(1, N);

for array_id = 1:N

    filename = sprintf('Stability_peak_isotropy_nu_ext%.2f_g%.2f_v%.1f.mat', ...
                        nuext, g, List_velocities(array_id));

    if ~isfile(filename)
        warning('File not found: %s', filename);
        continue;
    end

    data = load(filename);
    list_k = data.List_ky;
    local_minima_all = data.local_minima_all;

    % Extract representative eigenvalue (max Re part)
    real_part_curve = nan(size(list_k));

    for i = 1:length(list_k)
        if ~isempty(local_minima_all{i})
            real_part_curve(i) = max(local_minima_all{i}(:,1));
        end
    end

    % Keep only valid points for interpolation
    valid = ~isnan(real_part_curve) & ~isnan(list_k);

    k_orig = list_k(valid);
    real_orig = real_part_curve(valid);

    % Require enough points for interpolation
    if numel(k_orig) < 2
        warning('Not enough valid points for v = %.1f', List_velocities(array_id));
        continue;
    end

    k_fine = linspace(k_orig(1), k_orig(end), pt_more);
    real_interp = interp1(k_orig, real_orig, k_fine, 'pchip');

    % Maximum location
    [max_val, max_idx] = max(real_interp);

    plot(k_fine, real_interp, 'LineWidth', 2, ...
         'Color', cmap_orange2purple(array_id,:));

    plot(k_fine(max_idx), max_val, '.', ...
         'Color', cmap_orange2purple(array_id,:), ...
         'MarkerSize', 30);

    k_max(array_id) = k_fine(max_idx);

end

axis square;
set(gca, 'FontSize', 18);

figure;
plot(List_velocities./1000, k_max, '.', 'MarkerSize', 30);
xlabel('v [mm/s]');
ylabel('k^* [mm]');
axis square;
set(gca, 'FontSize', 18);