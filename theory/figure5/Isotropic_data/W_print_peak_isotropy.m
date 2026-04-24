 
% PLOT_DISPERSION_FROM_QUADRANT
%
% Reconstructs Re(λ(kx, ky)) from upper-quadrant data using symmetry.
% Handles duplicate (kx, ky) points explicitly and allows controlled
% aggregation of multiple minima per point.
%
% REQUIREMENTS:
%   - Files: Stability_peak_isotropy_..._htm%d.mat
%   - Variable: local_minima_all (cell array)
%       each cell: [n_minima × ...], col 1 = Re(λ)
%
% KEY FEATURE:
%   Explicit aggregation over duplicate points using chosen operator.

%% USER CHOICE: aggregation method
% Options: @mean, @min, @max
agg_fun = @min;   % <-- CHANGE THIS depending on physics

%% Parameters
Nx      = 11;
N_ky    = 11;
kx_list = linspace(0, 1.5, Nx);
ky_list = linspace(0, 1.5, N_ky);

file_template = 'Stability_peak_isotropy_nu_ext0.35_g4.10_v300.0_htm%d.mat';

%% Collect raw data
kx_all = cell(Nx,1);
ky_all = cell(Nx,1);
val_all = cell(Nx,1);

for ix = 1:Nx
    kx = kx_list(ix);
    data = load(sprintf(file_template, ix));

    kx_tmp = [];
    ky_tmp = [];
    val_tmp = [];

    for iy = 1:N_ky
        minima = data.local_minima_all{iy};
        n_min  = size(minima, 1);

        kx_tmp  = [kx_tmp; repmat(kx, n_min, 1)];
        ky_tmp  = [ky_tmp; repmat(ky_list(iy), n_min, 1)];
        val_tmp = [val_tmp; minima(:,1)];
    end

    kx_all{ix}  = kx_tmp;
    ky_all{ix}  = ky_tmp;
    val_all{ix} = val_tmp;
end

% Concatenate
kx_all  = vertcat(kx_all{:});
ky_all  = vertcat(ky_all{:});
val_all = vertcat(val_all{:});

%% Apply symmetry (4 quadrants)
kx_full  = [ kx_all; -kx_all;  kx_all; -kx_all ];
ky_full  = [ ky_all;  ky_all; -ky_all; -ky_all ];
val_full = repmat(val_all, 4, 1);

%% --- REMOVE DUPLICATES WITH CONTROLLED AGGREGATION ---

pts = [kx_full, ky_full];

% Unique spatial locations
[pts_unique, ~, ic] = unique(pts, 'rows');

% Aggregate values per (kx, ky)
val_unique = accumarray(ic, val_full, [], agg_fun);

kx_clean = pts_unique(:,1);
ky_clean = pts_unique(:,2);

%% Interpolation grid
nGrid = 200;
xq = linspace(-1.5, 1.5, nGrid);
yq = linspace(-1.5, 1.5, nGrid);
[Xq, Yq] = meshgrid(xq, yq);

Zq = griddata(kx_clean, ky_clean, val_unique, Xq, Yq, 'cubic');

% Handle NaNs (outside convex hull)
Zq(isnan(Zq)) = 0;

%% Plot
figure;
imagesc(xq, yq, Zq);
set(gca, 'YDir', 'normal');
axis square;

xlim([-1 1]);
ylim([-1 1]);

% Color scaling restricted to visible region
cols = xq >= -1 & xq <= 1;
rows = yq >= -1 & yq <= 1;
Zclip = Zq(rows, cols);

caxis([min(Zclip(:)), max(Zclip(:))]);


% Labels
cb = colorbar;
cb.Label.String = 'Re \lambda';

xlabel('u (dir max spread)');
ylabel('v (dir min spread)');
set(gca, 'FontSize', 16);
 