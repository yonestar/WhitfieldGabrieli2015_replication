figdir = '/beegfs/pl-active/ics/data/projects/wagerlab/labdata/collab/SADTRP_reboot/figs';

% load masks
maskdir = '/beegfs/pl-active/ics/data/projects/wagerlab/labdata/collab/SADTRP_reboot/masks';
brain_mask = fullfile(maskdir,'mask.volume.brainmask.nii');
amy = fmri_data(fullfile(maskdir, 'Amy_mask_WFU.nii'));%, brain_mask);

% load preproc'd data
datadir = '/beegfs/pl-active/ics/data/projects/wagerlab/labdata/collab/SADTRP_reboot/data';
Ss = filenames(fullfile(datadir, 'preprocessed', 'niftiDATA*'));

% load in behav data
tabl = readtable(fullfile(datadir, 'behav_data_cleaned.csv'));


%% load in corr maps. if file is found, skip next cell

load(fullfile(datadir, 'Amy_corr_maps.mat'))


%%  GENERATE AMY CONNECTIVITY "R-MAPS"    only if above load command fails
for i=1:length(Ss)
     
    % load data
    subjdat = fmri_data(Ss{i}, brain_mask);
    
    amy_timecourse = subjdat.extract_roi_averages(amy); % extract Amy timecourse
    amycorrs = corr(amy_timecourse.dat, subjdat.dat'); % compute corr between Amy and all other voxels
    %amycorrs(isnan(amycorrs)) = 0; % SOME VOXELS ARE ALL ZEROS (HAVE NO VARIANCE) AND YIELD A NAN  WHEN CORRELATED W/ AMY TIMECOURSE.  ALL THOSE NANS MAY LATER BE MESSING IT UP
    
    if i == 1 % first iteration, save back into dat
        amy_corr_dat = subjdat;
        amy_corr_dat.dat = atanh(amycorrs)'; % Fisher Z transform
    else % append
        amy_corr_dat.dat(:,i) = atanh(amycorrs)'; % Fisher Z transform
    end
end

save(fullfile(datadir, 'Amy_corr_maps.mat'), 'amy_corr_dat')

return


%% drop Ss with no behav data

wh = ~isnan(tabl.LSASpostCBTall_LOCF);
tabl = tabl(wh,:);
amy_corr_dat.dat = amy_corr_dat.dat(:,wh);

%%  CV FIND CLUSTERS, PREDICT TX RESPONSE IN HELD OUT
K = 5;  % 5-fold
indices = crossvalind('Kfold', height(tabl),K);

amy_corr_dat.X = [(tabl.LSASpostCBTall_LOCF - tabl.BaselineLSAS) tabl.BaselineLSAS double(categorical(tabl.Group))]; % delta LSAS, Baseline LSAS -- as in original report

%%
for i=1:K
    
    % set up training/test data
    train = get_wh_image(amy_corr_dat, indices ~= i);
    test = get_wh_image(amy_corr_dat, indices == i);
    
    % fit brain model - estimate how predictive of deltaLSAS each voxel's
    % conn w/ amy is
    out = regress(train, .001, 'unc', 'nodisplay');
    
    % threshold for clusters
    cl_train = threshold(get_wh_image(out.t, 1), .001, 'unc', 'k', 10); % arb cluster thrsh.  w/ k=822 didn't find anything
    cl_vals_train = train.extract_roi_averages(cl_train);
    
    % fit clusters model
    X = cat(2,cl_vals_train.dat); % predictors
    mdl = fitglm([baselineLSAS(wh)' X], deltaLSAS(wh));  
    
    % apply brain model to test data
    cl_vals_test = test.extract_roi_averages(cl_train);
    
    % apply cluster model to test data
    predLSAS(i) = mdl.Coefficients.Estimate(1) + (mdl.Coefficients.Estimate(2)*baselineLSAS(i)) + (cat(2,cl_vals_test.dat) * mdl.Coefficients.Estimate(3:end));

end

%% correlate predicted and actual LSAS to get an R value; compare to compact model
