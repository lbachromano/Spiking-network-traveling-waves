 
% PLOT_ANISOTROPIC_DISPERSION
%
% Reconstructs Re(λ(kx, ky)) from upper-quadrant data assuming symmetry:
%   Re(λ(kx, ky)) = Re(λ(±kx, ±ky))
%
% Also extracts the location of maximal Re(λ) in a consistent way.
%
% INPUT FILES:
%   Stability_n_anisotropy_<rho>_<rho_e>_nu_ext0.35_g<g>_v300.0_kx<kx>.mat
%
% Each file must contain:
%   - local_minima_all (cell array over ky)
%   - each cell: [n_minima × ...], column 1 = Re(λ)
%
% KEY FEATURE:
%   Explicit aggregation over duplicate (kx, ky) points.

%% USER CHOICE: aggregation rule
agg_fun = @min;   % typical: @min (most unstable mode)

%% Parameters
rho   = 0.0;
rho_e = 0.7;
g     = 4.1;

Nx      = 21;
N_ky    = 21;

kx_list = linspace(0, 1.5, Nx);
ky_list = linspace(0, 1.5, N_ky);

file_template = ...
    'Stability_n_anisotropy_%.1f_%.1f_nu_ext0.35_g%.2f_v300.0_kx%.2f.mat';

%% Collect data (cell-based to avoid reallocations)
kx_all = cell(Nx,1);
ky_all = cell(Nx,1);
val_all = cell(Nx,1);

for ix = 1:Nx
    kx = kx_list(ix);

    fname = sprintf(file_template, rho, rho_e, g, kx);
    data  = load(fname);

    kx_tmp = [];
    ky_tmp = [];
    val_tmp = [];

    for iy = 1:N_ky
        minima = data.local_minima_all{iy};
        n_min  = size(minima,1);

        kx_tmp  = [kx_tmp; repmat(kx, n_min, 1)];
        ky_tmp  = [ky_tmp; repmat(ky_list(iy), n_min, 1)];
        val_tmp = [val_tmp; minima(:,1)];
    end

    kx_all{ix}  = kx_tmp;
    ky_all{ix}  = ky_tmp;
    val_all{ix} = val_tmp;
end

% Concatenate once
kx_all  = vertcat(kx_all{:});
ky_all  = vertcat(ky_all{:});
val_all = vertcat(val_all{:});

%% Apply symmetry (4 quadrants)
kx_full  = [ kx_all; -kx_all;  kx_all; -kx_all ];
ky_full  = [ ky_all;  ky_all; -ky_all; -ky_all ];
val_full = repmat(val_all, 4, 1);

%% Remove duplicates with controlled aggregation
pts = [kx_full, ky_full];

[pts_unique, ~, ic] = unique(pts, 'rows');
val_unique = accumarray(ic, val_full, [], agg_fun);

kx_clean = pts_unique(:,1);
ky_clean = pts_unique(:,2);

%% Interpolation
nGrid = 200;
xq = linspace(-2, 2, nGrid);
yq = linspace(-2, 2, nGrid);
[Xq, Yq] = meshgrid(xq, yq);

Zq = griddata(kx_clean, ky_clean, val_unique, Xq, Yq, 'cubic');
Zq(isnan(Zq)) = 0;

%% Plot
figure;
imagesc(xq, yq, Zq);
set(gca, 'YDir', 'normal');
axis square;

xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

% Color scaling (visible region only)
cols = xq >= -1.5 & xq <= 1.5;
rows = yq >= -1.5 & yq <= 1.5;
Zclip = Zq(rows, cols);

caxis([min(Zclip(:)), max(Zclip(:))]);

cb = colorbar;
cb.Label.String = 'Re \lambda';

xlabel('u (dir max spread)');
ylabel('v (dir min spread)');
set(gca,'FontSize', 16);

%% ---------------------------------------------------------
% MAXIMUM DETECTION (CONSISTENT WITH AGGREGATED FIELD)
%% ---------------------------------------------------------

tol_abs = 1e-3;

max_value = max(val_unique);

max_idx = find(abs(val_unique - max_value) <= tol_abs);

u_max = kx_clean(max_idx);
v_max = ky_clean(max_idx);

fprintf('Number of max points within tolerance: %d\n', numel(u_max));

disp('Locations of maxima (kx, ky):');
disp([u_max, v_max]);

 