%% Stability analysis: smoothed eigenvalue curves vs external drive
% Each curve plots the real value of the eigenvalue;
% different colors corresponds to different values of the wavenumber k 

clear
Set_parameters;

% list_varying_param was computed in units of nu_ext/param.nutheta;
% rescale to the population average of nu_{theta} (see eq 21 from the paper) for physical interpretation
ratio_nu_theta = param.nutheta / param.nt_pop_average;
 
%% Parameters
N = 51;
List_k = linspace(0.0, 1.5, N);

% Colormap (orange -> purple)
orange = [1, 0.8, 0];
purple = [0.5, 0, 0.5];
cmap = colorGradient(orange, purple, N);

%% Figure
figure; hold on;

for kx = 1:N-3
    % --- Load and rescale parameter axis ---
    ld = load(sprintf('Stability_g41_fwd_%d.mat', kx));
    list_param_rescaled = ld.list_varying_param * ratio_nu_theta;  % <-- rescale once

    % --- Aggregate minima across parameters ---
    all_param  = [];
    all_minima = [];
    for i = 1:length(ld.local_minima_all)
        min_h = ld.local_minima_all{i};
        if isempty(min_h), continue; end
        nmin = size(min_h, 1);
        all_param  = [all_param;  list_param_rescaled(i) * ones(nmin,1)];
        all_minima = [all_minima; min_h];
    end

    % --- Filter unstable branch ---
    mask       = all_minima(:,3) <= 1e-4;
    all_param  = all_param(mask);
    all_minima = all_minima(mask,:);
    if isempty(all_param), continue; end

    % --- Extract observable ---
    [xu, iu] = unique(all_param, 'stable');
    yu = all_minima(iu, 1);
    if numel(xu) < 5, continue; end

    % --- Interpolation + smoothing ---
    xq       = linspace(xu(1), xu(end), 1000);
    y_smooth = smoothdata(interp1(xu, yu, xq, 'makima'), 'sgolay', 100);

    % --- Map k -> color ---
    k_val = List_k(kx);
    idx   = round((k_val - min(List_k)) / (max(List_k) - min(List_k)) * (N-1)) + 1;
    idx   = max(1, min(N, idx));

    plot(xq, y_smooth, 'Color', cmap(idx,:), 'LineWidth', 2);
end

%% Reference lines and axes
xline(0.31  * ratio_nu_theta, '--k', 'LineWidth', 1);
xline(0.361 * ratio_nu_theta, '--k', 'LineWidth', 1);
yline(0,                      '--k', 'LineWidth', 1);

xlim([0.27, 0.7] * ratio_nu_theta);
xlabel('\nu_{ext}/\nu_{\theta}');
ylabel('Re \lambda');
axis square;
set(gca, 'FontSize', 18);

%% Colorbar
colormap(cmap);
caxis([min(List_k) max(List_k)]);
cb = colorbar;
cb.Label.String = 'k [mm]';
cb.FontSize = 16;
box on;