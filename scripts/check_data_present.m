basedir = '~/Dropbox/Anx_Goldin_Gross_Anxiety_data/';

% load rs-fMRI scans
funcs = filenames(fullfile(basedir,'BaselineResting','1*nii*'));
func_subjs = char(funcs);
func_subjs = str2num(func_subjs(:,67:70));

% load anatomicals
anats = filenames(fullfile(basedir,'data','72_SP','*'));
anat_subjs = char(anats);
anat_subjs = str2num(anat_subjs(:,end-3:end));

%% find Ss with a functional but not structural

[a,b]=setdiff(func_subjs, anat_subjs);

% just one is missing -- drop the missing S
funcs(b) = [];
func_subjs(b) = [];

%% read in behavioral data
tabl = readtable(fullfile(basedir,'data','IndividualCBTvsWLstudy_LSASbaselineT2andPostCBTall_forYoni.csv'));

% keep patients only
tabl = tabl(tabl.ID < 2000,:);

% do all func Ss have behav data? Yes.
setdiff(func_subjs,tabl.ID)