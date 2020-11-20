%% find preproc data

datadir = '/beegfs/pl-active/ics/data/projects/wagerlab/labdata/collab/SADTRP_reboot/data';
Ss = filenames(fullfile(datadir, 'preprocessed', '*nii'));

%% prep behav dat

tablxls = readtable(fullfile(datadir,'IndividualCBTvsWLstudy_LSASbaselineT2andPostCBTall_forYoni.xlsx'));

% drop healthies
tablxls(tablxls.ID > 2000,:) = [];

tablxls.Properties.VariableNames{4} = 'BaselineLSAS';
tablxls.Properties.VariableNames{7} = 'LSASpostCBTall_LOCF';

% LSAS_T2: LSAS score at 2nd timepoint, after tx or wl
% LSASpostCBTall: LSAS score after CBT, whether administered immediately or
% after wl period

head(tablxls)

%% find RS data, add column

rawdata_Ss = filenames(fullfile(datadir, 'raw', '*'), 'char');
rawdata_Ss = str2num(rawdata_Ss(:,end-3:end));

sum(ismember(rawdata_Ss, tablxls.ID))

% only keep Ss with brain data
tablxls = tablxls(ismember(tablxls.ID, rawdata_Ss),:);

head(tablxls)

sum(~strcmp(tablxls.Group,'WL'))

writetable(tablxls, fullfile(datadir, 'behav_data_cleaned.csv'))