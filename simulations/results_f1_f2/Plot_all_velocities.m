 
%% BOXPLOT

wv_data=[];
%Load data: velocity from Macaque recordings
ld = load(sprintf('Macaque_R_trials_100_DWfilter_3_thPGD_%.2f_the_20.mat',threshold_input));
Macaque_data = ld.Macaque_final;
wv_data  = Macaque_data.all_speeds;
 
%%

titlefile=sprintf('Waves_DWfilter_%d_thPGD_%.2f_the_%d.mat',width_filter,threshold_input,th_e);
ld = load(titlefile);
S_data_short = ld.Simulation_final;

wl_short = cell2mat({S_data_short.mean_wavelengths}');
wv_short =S_data_short.all_speeds;
%wv_short =S_data_short .mean_speeds;

 
 %%
  
edges = linspace(0, max([wv_data(:); wv_short(:)]), 101); % 100 bins

figure; hold on
histogram(wv_data,  edges, 'Normalization','pdf', ...
    'FaceColor','k','FaceAlpha',0.35);
histogram(wv_short, edges, 'Normalization','pdf', ...
    'FaceColor','b','FaceAlpha',0.5);
legend('Sata','Simulations'); box on
xlabel('Wave Velocity (cm/s)')
ylabel('PDF')
xlim([0 70])
 %%
% Combine velocities and create a cell?array grouping variable
vel_all   = [wv_data; wv_short];
group_vel = [ ...
    repmat({'Data'},   numel(wv_data),1); ...
    repmat({'Sim'}, numel(wv_short),1); ...
    ];
ax = gca;
ax.FontSize = 15;
% Plot box?plot for velocities
figure;
 
boxplot(vel_all, group_vel, ...
    'Notch',           'on', ...
    'Labels',          {'Data','Sim'}, ...
    'LabelOrientation','inline', ...
    'Positions',       [1, 1.3], ...    % pull the two boxes closer
    'Widths',          0.25);   


ylabel('Wave Velocity (cm/s)')
set(gca, 'FontSize', 15)
% Change group label font size
h = findobj(gca, 'Tag', 'Box');  % Forces the boxplot to fully render
drawnow                          % Ensures all graphics objects are available
xt = findall(gca, 'Type', 'text');  % Find x-axis text labels
set(xt, 'FontSize', 15)
xtickangle(45)
% Find all text objects in current axes
textHandles = findall(gca, 'Type', 'text');

% Filter by x-position to find x-axis group labels (usually at y=0 or near)
for i = 1:length(textHandles)
    if textHandles(i).Position(2) < min(ylim) + 0.05 * range(ylim)  % near x-axis
        textHandles(i).Rotation = 45;  % rotate 45 degrees
        textHandles(i).HorizontalAlignment = 'right';  % better alignment
        textHandles(i).FontSize = 15;
        
    end
end
grid on
ylim([0 70])
pbaspect([1 1 1])


