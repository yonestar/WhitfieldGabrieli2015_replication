% copy resting state runs into anatomical directories

basedir = '/pl/active/ics/data/projects/wagerlab/labdata/collab/SADTRP_reboot';
cd(fullfile(basedir,'data','raw'));

% for each func run, find its folder and copy it
funcnames = filenames(fullfile('BaselineResting','*'));

for i=1:length(funcnames)
    s = funcnames{i};
    
    copyfile(s, s(17:20))
end

%% delete folders with just anatomicals, no functionals

funcnames = filenames('1*');

for i=1:length(funcnames)
    
    if ~exist(filenames(fullfile(funcnames{i},'*baseline*'),'char'))
        delete(filenames(fullfile(funcnames{i},'anat*'),'char')) % delete the anatomical
        rmdir(funcnames{i}) % delete the dir
    end
    
end