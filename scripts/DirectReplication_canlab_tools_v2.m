%% Attempt to validate W-G brain biomarkers for Social Anxiety Disorder treatment outcome on a novel dataset
% This script conducts a direct replication of the Whitfield-Gabrieli 
% amygdala connectivity analyses. 
% It depends on SPM12 and the Canlab_Tools matlab libraries. 
% It does a direct replication -- using the clusters and the parameters from
% the original report and also estimates new parameter estimates.


%% Load the data

clear all

repodir = '/Users/yoni/Repositories/WhitfieldGabrieli2015_replication';
maskdir = fullfile(repodir, 'masks');
resultsdir = fullfile(repodir, 'results');
datadir = fullfile(repodir, 'data');
figdir = fullfile(repodir, 'figs');

tabl = readtable(fullfile(datadir, 'behav_data_cleaned.csv'));

%% load masks and seeds

% load in Amygdala WFU seed
brain_mask = fullfile(maskdir, 'mask.volume.brainmask.nii');
amy = fmri_data(fullfile(maskdir, 'Amy_mask_WFU.nii'));

% load in the W-G Amy "target regions"
pos_mask = fmri_data(fullfile(maskdir, 'mask_AMYG_targetsPos.ROIs.img'));
neg_mask = fmri_data(fullfile(maskdir, 'mask_AMYG_targetsNeg.ROIs.img'));

%% load in Amy corr maps

%load(fullfile(datadir, 'Amy_corr_maps.mat'));

% these are NOT already fisher transformed
amy_corr_dat = fmri_data(fullfile(datadir, 'Amy_corr_maps.nii'));

%% direct replication: extract seed and target values, compute Zpos and Zneg

% Extract the BOLD signal of the amygdala and the target regions provided by W-G.
pos = amy_corr_dat.extract_roi_averages(pos_mask, 'contiguous_regions');  % IF WANT TO DO FOR EVERY VOXEL, CHECK OUT APPLY_MASK('CORRELATION') 
neg = amy_corr_dat.extract_roi_averages(neg_mask, 'contiguous_regions');

Zpos = fisherz(pos.dat);
Zneg = mean(fisherz(cat(2,neg.dat)), 2);

% compute AMYG = Zpos - Zneg;
AMYG = Zpos - Zneg;
AMYGz = zscore(AMYG);

%save(matfile_fname, 'AMYG', 'AMYGz', 'MVPA', 'MVPAz');



%% make deltaLSAS score -- doing baseline - LOCF, for all Ss. This way change is a positive number, like they did it
tabl.deltaLSAS = tabl.BaselineLSAS - tabl.LSASpostCBTall_LOCF;

% drop Ss for whom we don't have a delta LSAS. useless
todrop = isnan(tabl.deltaLSAS);
tabl(todrop,:) = [];
height(tabl)

% remove intercept, for simplicity and fairness
tabl.deltaLSAS_mc = tabl.deltaLSAS - nanmean(tabl.deltaLSAS);
tabl.BaselineLSAS_mc = tabl.BaselineLSAS - nanmean(tabl.BaselineLSAS);

% Zscore the computed connectivity values - one more time after dropping
% some Ss
tabl.AMYGz = zscore(AMYG(~todrop));
    
%% test for normality

[h, p, adstat]=adtest(tabl.deltaLSAS)
[h, p, adstat]=adtest(tabl.AMYGz)
[h, p, adstat]=adtest(tabl.BaselineLSAS)

%% compact model -- predicting deltaLSAS from initialLSAS only

clc

% estimated on my data
mdl_initialLSAS = fitglm(tabl.BaselineLSAS_mc, tabl.deltaLSAS_mc, 'VarNames', {'BaselineLSAS' 'DeltaLSAS'}, 'intercept', false)
disp(mdl_initialLSAS.Rsquared)
MSE_full = mdl_initialLSAS.SSE / 42

% their param for Initial LSAS
tabl.pred_deltaLSAS_baselineLSAS = 0.6194*tabl.BaselineLSAS_mc;
r2theirBaselineLSAS = corr(tabl.pred_deltaLSAS_baselineLSAS, tabl.deltaLSAS, 'rows', 'complete') ^ 2
MSE_compact = mean(( tabl.pred_deltaLSAS_baselineLSAS - tabl.deltaLSAS_mc) .^ 2)

predictR2_initialLSAS_OLS = 1 - (MSE_full / MSE_compact)

%% test Amy direct replication: their clusters, their param estimates
clc

% Reported by W-G: LSAS_change = 0.6194*LSAS_pre +  8.6290*AMYG - 9.9763   (R2=.33, T(35)=3.45, p=0.001)
tabl.pred_deltaLSAS_amy = 0.6194*tabl.BaselineLSAS_mc + 8.6290*tabl.AMYGz; % - 9.9763;

MSE_full = mean(( tabl.pred_deltaLSAS_amy - tabl.deltaLSAS_mc) .^ 2);
predictR2 = 1 - (MSE_full / MSE_compact)

% No need for adjusted Rsquared -- since we are not fitting
% additional params, don't want to penalize for that
model_r2theirPEs = corr(tabl.pred_deltaLSAS_amy, tabl.deltaLSAS, 'rows', 'complete') ^ 2;
modelR2 = model_r2theirPEs - r2theirBaselineLSAS

%% Amy data analytic replication #1: their clusters, estimate new param estimates. 
% THERE DOES NOT SEEM TO BE BETTER PARAM ESTIMATES. DROP FOR NOW.
% re-center AMYGz b/c dropped some Ss
mdl_newPEs = fitglm([tabl.BaselineLSAS_mc, tabl.AMYGz], tabl.deltaLSAS_mc, 'VarNames', {'baselineLSAS', 'AMYGz', 'deltaLSAS'}, 'intercept',false)
mdl_newPEs.Rsquared 


%% figure panel a

create_figure('fig 1', 1, 1)

fs = 88;
scatter(tabl.pred_deltaLSAS_amy , tabl.deltaLSAS_mc, 'b', 'SizeData', fs);
scatter(tabl.pred_deltaLSAS_baselineLSAS , tabl.deltaLSAS_mc, 'r', 'SizeData', fs);

h = lsline; set(h(1), 'color', 'b'); set(h(2), 'color', 'r');
set(gca, 'FontSize', 20)
ylabel('Observed treatment response')
xlabel('Predicted treatment response')
legend('full model', 'compact model', 'Location', 'NW', 'FontSize', 16)

saveas(gcf, fullfile(figdir, 'Fig1c.png'))

%% figure panel b
create_figure('PEs')

% error in prediction
fullmodel_error = tabl.deltaLSAS_mc - tabl.pred_deltaLSAS_amy;
compactmodel_error = tabl.deltaLSAS_mc - tabl.pred_deltaLSAS_baselineLSAS;
diffs = fullmodel_error - compactmodel_error;
%histogram(diffs , 10)%, 'kx', 'LineWidth', 1, 'SizeData', fs*2)
bar(sort(diffs))
ylabel('Change in prediction error')
xlabel('Subject #')
set(gca, 'FontSize', 20)

saveas(gcf, fullfile(figdir, 'Fig1d.png'))


%% for CBT immediate only, CBT-postWL only
% wh = strcmp(tabl.Group,'WL');
% corr(tabl.pred_deltaLSAS_amy(wh), tabl.deltaLSAS(wh), 'rows', 'complete') ^ 2
% corr(tabl.pred_deltaLSAS_amy(~wh), tabl.deltaLSAS(~wh), 'rows', 'complete') ^ 2


%% Permutation test on whether additional variance explained by Amy is sig

nPerm = 10000;
modelR2_perm = nan(nPerm,1);
predictR2_perm = nan(nPerm,1);

%tabl2 = tabl(~wh,:);

for i=1:nPerm
    
    % UNCOMMENT ONE OF THESE FOUR OPTIONS
    
    % FOR ALL Ss
    pred = 0.6194*tabl.BaselineLSAS_mc + 8.6290*(tabl.AMYGz(randperm(height(tabl))));     % shuffle Amy-conn values, then apply their model
    
    MSE_full = mean(( pred - tabl.deltaLSAS_mc) .^ 2);
    predictR2_perm(i) = 1 - (MSE_full / MSE_compact);
    
    modelR2_perm(i) = (corr(pred, tabl.deltaLSAS) ^ 2) - r2theirBaselineLSAS;
    
    % FOR CBT-IMMEDIATE ONLY
    %pred = 0.6194*tabl2.BaselineLSAS_mc + 8.6290*(tabl2.AMYGz(randperm(height(tabl2)))); 
    %var_exp_perm(i) = (corr(pred, tabl2.deltaLSAS)); % ^ 2;
    
end



%% plot perm test results

create_figure('perm test - corr', 1, 2);

% model r squared
histogram(modelR2_perm)
%vline(([r2theirBaselineLSAS r2theirPEs prctile((var_exp_perm .^ 2), 95)]), {'r' 'b' 'k'}, {'BaselineLSAS only' 'observed R^2' '95th prctile of null distrb'})
vline(prctile(modelR2_perm, 95), 'k', { '95th prctile of null distrb'})
vline(modelR2, 'b', { 'observed result'})
xlabel('model R^2'), ylabel('# of permutations')

% predict r squared
subplot(1,2,2)
histogram(predictR2_perm)
%vline(sqrt([r2theirBaselineLSAS r2theirPEs]), {'r' 'b'}, {'BaselineLSAS only' 'observed R^2'})
vline(prctile(predictR2_perm, 95), 'k', { '95th prctile of null distrb'})
vline(predictR2, 'b', { 'observed result'})
xlabel('predict R2'), ylabel('# of permutations')

%% p values

% our p value, model R2
sum(modelR2 < modelR2_perm) / nPerm

% our p value, predict R2
sum(predictR2 < predictR2_perm) / nPerm

