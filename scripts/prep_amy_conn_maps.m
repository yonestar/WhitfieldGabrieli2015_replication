clear all

repodir = '/Users/yoni/Repositories/WhitfieldGabrieli2015_replication';

Ss = filenames(fullfile(repodir, 'full_data', 'preprocessed', '*nii'));

% load Amy mask
amy = fmri_data(fullfile(repodir, 'masks', 'Amy_mask_WFU.nii'));

% make "template" for corr maps to save into
amy_corr_maps = fmri_data(Ss{1});
amy_corr_maps.dat = amy_corr_maps.dat(:,1);

%% for each S
for i=1:length(Ss)
    
    % load in the full timeseries data
    dat = fmri_data(Ss{i});
    
    % BPF .01 - .1 -- already done in preproc.

    % compute voxwelwise amy corr map. Don't Fisher Z -- that comes later.
    % Should then be able to run direct rep v2 and reproduce prev values
    amy_timecourse = dat.extract_roi_averages(amy);
    
    % voxelwise corr. then append to 4D corr maps
    amy_corr_maps.dat(:,i) = corr(amy_timecourse.dat, dat.dat');   
   
end

%% save to nifti
write(amy_corr_maps, 'fname', fullfile(repodir, 'data', 'Amy_corr_maps.nii'))