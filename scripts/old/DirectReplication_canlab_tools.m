%% Attempt to validate W-G brain biomarkers for Social Anxiety Disorder treatment outcome on a novel dataset
% This script conducts a direct replication of the Whitfield-Gabrieli analyses
% including both the Amygdala and the MVPA analyses. It depends on SPM12 and the Canlab_Tools matlab libraries. Does a direct replication -- using the clusters and the parameters from
% the original report also tests "Data-Analytic Replication #1": original
% clusters, but estimates new parameter estimates.


%% Load the data

clear all

repodir = '/Users/yoni/Repositories/WhitfieldGabrieli2015_replication';
maskdir = fullfile(repodir, 'masks');
resultsdir = fullfile(repodir, 'results');
datadir = fullfile(repodir, 'data');
figdir = fullfile(repodir, 'figs');

tabl = readtable(fullfile(datadir, 'behav_data_cleaned.csv'));

load(fullfile(datadir, 'Amy_corr_maps.mat'));

%% load masks and seeds

%{
files in the 'masks' directory from the original W-G paper include:

mask_AMYG_targetsNeg/Pos:  Show the regions depicted in Figure 1.
Note that caudate/sgACC region was negatively connected to Amy, while other
regions were positively connected.  

mask_MVPA_seeds.ROIs.img:  Figure 2., left side.  contains two regions:
PCC and cerebellum.

mask_MVPA_targets.ROIs.img:  Figure 2., right side.  
%}

% load in Amygdala WFU seed
brain_mask = fullfile(maskdir, 'mask.volume.brainmask.nii');
amy = fmri_data(fullfile(maskdir, 'Amy_mask_WFU.nii'));

% load in the W-G Amy "target regions"
pos_mask = fmri_data(fullfile(maskdir, 'mask_AMYG_targetsPos.ROIs.img'));
neg_mask = fmri_data(fullfile(maskdir, 'mask_AMYG_targetsNeg.ROIs.img'));

% load the W-G MVPA seeds and targets
% MVPA_seeds = fmri_data(fullfile(maskdir, 'mask_MVPA_seeds.ROIs.img'));
% MVPA_targets = fmri_data(fullfile(maskdir, 'mask_MVPA_targets.ROIs.img'));
% 
% % separate the two MVPA seed regions
% tmp = region(MVPA_seeds);
% MVPAseed1 = region2fmri_data(tmp(1), MVPA_seeds);
% MVPAseed2 = region2fmri_data(tmp(2), MVPA_seeds);

% uncomment to view the clusters
%orthviews_multiple_objs({amy, pos_mask, neg_mask, MVPA_seeds, MVPA_targets});


%% direct replication: extract seed and target values, compute Zpos and Zneg

matfile_fname = fullfile(resultsdir, 'direct_rep.mat');

% if have previously run this, can just load the matfile 
if exist(matfile_fname, 'file')
    load(matfile_fname);
else % if not, re-run
     
    % AMY
    for s=1:length(Ss)     
        
        % load in subject's preprocessed data, compute connectivity between
        % bilateral Amy and originally reported target regions
        dat = fmri_data(Ss{s});    
     
        % Extract the BOLD signal of the amygdala and the target regions provided by W-G.
        amy_timecourse{s} = dat.extract_roi_averages(amy);
        target_timecourse_pos{s} = dat.extract_roi_averages(pos_mask, 'contiguous_regions');  % IF WANT TO DO FOR EVERY VOXEL, CHECK OUT APPLY_MASK('CORRELATION') 
        target_timecourse_neg{s} = dat.extract_roi_averages(neg_mask, 'contiguous_regions');
        
        % compute Zneg = Average connectivity (Fisher transformed correlations) between Amyg seed and mask_AMYG_targetsNeg.ROIs.img
        Zneg(s) = mean(fisherz(corr(amy_timecourse{s}.dat, cat(2,target_timecourse_neg{s}.dat))));
        % compute Zpos = Average connectivity (Fisher transformed correlations) between Amyg seed and mask_AMYG_targetsPos.ROIs.img
        Zpos(s) = mean(fisherz(corr(amy_timecourse{s}.dat, cat(2,target_timecourse_pos{s}.dat))));% Extract the BOLD signal of the amygdala and the target regions provided by W-G.
        
        % compute AMYG = Zpos - Zneg;
        AMYG(s) = Zpos(s) - Zneg(s);
    
    
        % MVPA
        
        % Extract the BOLD signal of the MVPA generated seeds and targets provided by W-G.
%         mvpa_seed1_timecourse{s} = dat.extract_roi_averages(MVPAseed1, 'contiguous_regions');
%         mvpa_seed2_timecourse{s} = dat.extract_roi_averages(MVPAseed2, 'contiguous_regions');
%         mvpa_targets_timecourse{s} = dat.extract_roi_averages(MVPA_targets, 'contiguous_regions');
% 
%         % compute Zpos = Average connectivity (Fisher transformed correlations) between MVPA seed cluster #1 and mask_MVPA_targets.ROIs.img
%         % compute Zneg = Average connectivity (Fisher transformed correlations) between MVPA seed cluster #2 and mask_MVPA_targets.ROIs.img
%         Z1(s) = mean(fisherz(corr(mvpa_seed1_timecourse{s}.dat, cat(2,mvpa_targets_timecourse{s}.dat))));
%         Z2(s) = mean(fisherz(corr(mvpa_seed2_timecourse{s}.dat, cat(2,mvpa_targets_timecourse{s}.dat))));
% 
%         % compute MVPA = Zpos - Zneg;
%         MVPA(s) = Z1(s) - Z2(s);
    end

    
    MVPAz = (MVPA - nanmean(MVPA)) / nanstd(MVPA); 

    save(matfile_fname, 'AMYG', 'AMYGz', 'MVPA', 'MVPAz');
end


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
tabl.AMYGz = zscore(AMYGz(~todrop))';
    
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

predictR2 = 1 - (MSE_full / MSE_compact)

%% test Amy direct replication: their clusters, their param estimates

% Reported by W-G: LSAS_change = 0.6194*LSAS_pre +  8.6290*AMYG - 9.9763   (R2=.33, T(35)=3.45, p=0.001)
tabl.pred_deltaLSAS_amy = 0.6194*tabl.BaselineLSAS_mc + 8.6290*tabl.AMYGz; % - 9.9763;

MSE_full = mean(( tabl.pred_deltaLSAS_amy - tabl.deltaLSAS_mc) .^ 2)
predictR2 = 1 - (MSE_full / MSE_compact)

% No need for adjusted Rsquared -- since we are not fitting
% additional params, don't want to penalize for that
model_r2theirPEs = corr(tabl.pred_deltaLSAS_amy, tabl.deltaLSAS, 'rows', 'complete') ^ 2

% key results
clc

predictR2
modelR2 = model_r2theirPEs - r2theirBaselineLSAS





%% for CBT immediate only, CBT-postWL only
% wh = strcmp(tabl.Group,'WL');
% corr(tabl.pred_deltaLSAS_amy(wh), tabl.deltaLSAS(wh), 'rows', 'complete') ^ 2
% corr(tabl.pred_deltaLSAS_amy(~wh), tabl.deltaLSAS(~wh), 'rows', 'complete') ^ 2


%% Amy data analytic replication #1: their clusters, estimate new param estimates. Here we need to adjust R2
% THERE DOES NOT SEEM TO BE BETTER PARAM ESTIMATES. DROP FOR NOW.
% re-center AMYGz b/c dropped some Ss
mdl_newPEs = fitglm([tabl.BaselineLSAS_mc, tabl.AMYGz], tabl.deltaLSAS_mc, 'VarNames', {'baselineLSAS', 'AMYGz', 'deltaLSAS'}, 'intercept',false)
mdl_newPEs.Rsquared 

%%
% ypred = feval(mdl_newPEs, [tabl.BaselineLSAS_mc, tabl.AMYGz]);
% MSE_full = mean(( ypred - tabl.deltaLSAS_mc) .^ 2)
% predictR2 = 1 - (MSE_full / MSE_compact)

%  I THINK NEED TO CROSS-VAL


%tabl.pred_deltaLSAS_amy_new_pe = mdl_newPEs.Fitted.Response;
%tabl.pred_deltaLSAS_amy_new_pe_mc = scale(mdl_newPEs.Fitted.Response,1);


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


v = get(gca,'Position');
%set(gca,'Position',[v(1) v(2)*1.5 v(3:4)])

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

%% Permutation test on whether additional variance explained by Amy is sig

nPerm = 1000;
modelR2_perm = nan(nPerm,1);
predictR2_perm = nan(nPerm,1);

%tabl2 = tabl(~wh,:);

for i=1:nPerm
    
    % UNCOMMENT ONE OF THESE FOUR OPTIONS
    
    % FOR ALL Ss
    pred = 0.6194*tabl.BaselineLSAS_mc + 8.6290*(tabl.AMYGz(randperm(height(tabl))));     % shuffle Amy-conn values, then apply their model
    MSE_full = mean(( pred - tabl.deltaLSAS_mc) .^ 2);
    predictR2_perm(i) = 1 - (MSE_full / MSE_compact);
    modelR2_perm(i) = (corr(pred, tabl.deltaLSAS)) ^ 2;
    
    % FOR CBT-IMMEDIATE ONLY
    %pred = 0.6194*tabl2.BaselineLSAS_mc + 8.6290*(tabl2.AMYGz(randperm(height(tabl2)))); 
    %var_exp_perm(i) = (corr(pred, tabl2.deltaLSAS)); % ^ 2;
    
    % FOR AMYG MODEL ONLY, NO BASELINE LSAS, all Ss
    %pred = 8.6290*(tabl.AMYGz(randperm(height(tabl)))); 
    %var_exp_perm(i) = (corr(pred, tabl.deltaLSAS)); % ^ 2;
    
    % FOR AMYG MODEL ONLY, NO BASELINE LSAS, only CBT immediate Ss
    % pred = 8.6290*(tabl2.AMYGz(randperm(height(tabl2)))); 
    % var_exp_perm(i) = (corr(pred, tabl2.deltaLSAS)); % ^ 2;
    
end

%save(fullfile(resultsdir, 'direct_rep_permtest.mat'), 'var_exp_perm');


%% see results. GOOD FOR first and second code blocks above, see below for amyg model only

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




%% plot residuals
create_figure('resids', 1, 3); 
violinplot(resid_amy.^2 - resid_compact.^2); legend off, title('Prediction error differences in Amyg vs. compact models')
hline(0)
subplot(1,3,2)
scatter(tabl.pred_deltaLSAS_baselineLSAS, tabl.pred_deltaLSAS_amy), lsline
xlabel('Predictions from compact model'), ylabel('Predictions from amyg model')
subplot(1,3,3)
scatter(resid_compact.^2, resid_amy.^2)
xlabel('abs(errors) from compact model'), ylabel('abs(errors) from amyg model'), 
hold on, plot([0 4000], [0 4000], 'k')
lsline

%% characterize spread of difference in errors

sum(diffs < -5)
sum(diffs > 5)

sum(diffs < -6.56)
sum(diffs > 6.56)








%% is AMYGZ more related to pre, post, or change?

corr([tabl.AMYGz tabl.BaselineLSAS tabl.LSASpostCBTall_LOCF tabl.deltaLSAS])

%% MODEL COMPARISON - LRT
% L.R.T is appropriate due to Neyman-Pearson Lemma, I believe, since I did
% not estimate parameters on either model, both are just applied

% can submit the deviance values to lratiotest
% the deviance is the sum of the squared errors. checked in Matlab. #OK
dev1 = sum( (tabl.pred_deltaLSAS_baselineLSAS - tabl.deltaLSAS_mc) .^ 2)
dev2 = sum( (tabl.pred_deltaLSAS_amy - tabl.deltaLSAS_mc) .^ 2)

[~,p,stat,cvalue] = lratiotest(dev1, dev2, 50)

p = 1 - chi2cdf(dev1 -  dev2, 50)

%% MODEL COMPARISON -- COMBINED

PRE = (dev1 - dev2) / dev1

%% MODEL COMPARISON -- OLD
% Compare Model 2 and Model 1. Model 3 is no better than Model 2, and
% required estimation of new parameters

% first, mean center predicted improvement from both models
tabl.pred_deltaLSAS_amy_mc = tabl.pred_deltaLSAS_amy - nanmean(tabl.pred_deltaLSAS_amy);
tabl.pred_deltaLSAS_baselineLSAS_mc = tabl.pred_deltaLSAS_baselineLSAS - nanmean(tabl.pred_deltaLSAS_baselineLSAS);

SSEcompact = nansum((tabl.deltaLSAS_mc - tabl.pred_deltaLSAS_baselineLSAS_mc) .^ 2)
SSEaug = nansum((tabl.deltaLSAS_mc - tabl.pred_deltaLSAS_amy_mc) .^ 2)

PRE = (SSEcompact - SSEaug) / SSEcompact
%%
n = sum(~isnan(tabl.deltaLSAS));
F = PRE / ((1 - PRE) / (n-0))
p = 1 - fcdf(F, 1, 42)

% OK, so AMYGz helps, but increase is n.s. with p = .26


% now explain an additional 5% of the variance, p = .15

% FORMULA FOR ADJUSTED R2
% p = 0;
% r2adj = 1 - (1-r2)*((n-1) / (n-p-1))
% b = r2adj;


%% BAR PLOT: adjusted R2 values for all models

% for the model using their PEs, k = 0, so Adjusted R2 = R2

create_figure('R2');
%bar([mdl_initialLSAS.Rsquared.Adjusted r2theirPEs mdl_newPEs.Rsquared.Adjusted]')
bar([r2theirBaselineLSAS model_r2theirPEs]')
ylabel('R^2')
set(gca,'XTick',1:2,'XTickLabel',{'Baseline LSAS' 'Baseline LSAS + Amy'})
xtickangle(30)
saveas(gcf, fullfile(figdir, 'R2amy_models.svg'))%, ylim([.1 .27])
%% SCATTER PLOT

% R2 barplot
create_figure('amy direct'); plot_correlation_samefig(tabl.deltaLSAS, tabl.pred_deltaLSAS_amy,[],'ko',0,0); %lsline; title('Amy - Direct')
xlabel('Observed delta LSAS'), ylabel('Pred delta LSAS')
%[r,p] = corr(tabl.deltaLSAS, tabl.pred_deltaLSAS_amy, 'rows', 'complete') 


%% canlab variance decompisiton