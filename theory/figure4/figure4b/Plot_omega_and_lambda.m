% -------------------------------------------------------------------------
% SUMMARY
% This script aggregates simulation outputs across multiple values of the
% external drive parameter (nu), automatically extracted from filenames.
% For each condition, it loads wavelength and oscillation frequency data,
% computes summary statistics (mean and standard deviation), and compares
% them against mean-field theoretical predictions.
%
% The script produces two figures:
% 1) Spatial frequency (k) vs external drive, including error bars from
%    simulated wavelengths and mean-field predictions, with instability
%    region boundaries indicated.
% 2) Oscillation frequency (omega) vs external drive, comparing simulated
%    dominant frequencies to mean-field predictions, again highlighting
%    instability bounds.
% -------------------------------------------------------------------------

clear

% 1. List all files beginning with 'wavelengths_nu' (any extension)
files = dir('wavelengths_nu*.*');
% 2. Preallocate vector
nFiles = numel(files);
vals   = NaN(nFiles,1);
% 3. Extract the three decimal number after 'wavelengths_nu'
for k = 1:nFiles
    name = files(k).name;
    tok = regexp(name, '^wavelengths_nu([0-9]+\.[0-9]{3})', 'tokens');
    if ~isempty(tok)
        vals(k) = str2double(tok{1}{1});
    end
end
% 4. Remove any that didn not match and sort ascending
vals = vals(~isnan(vals));
vals = sort(vals);
% 5. (Optional) display
disp('Found wavelengths_nu values:');
disp(vals);


%%

all_d={};
mean_w=zeros(length(vals),1);
std_w=zeros(length(vals),1);
mean_om=zeros(length(vals),1);
std_om=zeros(length(vals),1);


for i=1:length(vals)
   
    ld=load(sprintf('wavelengths_nu%.3f.mat',vals(i)));   
    all_w{i}=(ld.wave_wavelengths);   
    mean_w(i)=mean(all_w{i},"omitnan");
    std_w(i)=std(all_w{i},"omitnan");
    ld=load(sprintf('oscil_freq_nu%.3f.mat',vals(i)));
    mean_om0(i)=ld.f0_LFP;
    mean_om1(i)=ld.f0_ave_LFP;
    std_om1(i)=ld.sd_LFP;
    mean_om2(i)=ld.domFreq_rate;
    std_om2(i)=ld.sigma_rate;
end

%%

bounds=load('Instability_region.mat','lb','ub');
ld=load('Wavelength_from_mean_field.mat','ext_field','k_max_list');

figure
plot(ld.ext_field, ld.k_max_list , 'k','linewidth',3);
hold on
errorbar(vals,2.*pi./(10.*mean_w),std_w.*2.*pi./(10.*mean_w.^2),'.','markersize',40) %let's use mm
hold on
xline(bounds.lb)
hold on
xline(bounds.ub)
hold on
xlim([0.75 1.05])
ylabel('k [mm^{-1}]','FontSize',18);
set(gca,'fontsize',18)


ld=load('Oscill_freq_from_mean_field.mat','','maxOmegas');

figure
plot(ld.ext_field,ld.maxOmegas,'k','linewidth',2)
hold on
errorbar(vals,mean_om2,std_om2,'.','markersize',40)
hold on
xline(bounds.lb)
hold on
xline(bounds.ub)
xlim([0.75 1.05])
xlabel('\nu_{ext}/\nu_{\theta}','FontSize',14);
ylabel('\omega [Hz]','FontSize',18);
set(gca,'fontsize',18)
 