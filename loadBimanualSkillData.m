% loads raw data from all subjects into a series of .mat data files (one per subject)
%
% This example loads a single 'chunk' (block of 5 trials)

clear all
% rot_vs_mir
names{1} = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
names{2} = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
groups = {'rot','mir'};
path = 'Data/rot_vs_mir/';
blocks = {'no_rot1','rot1','rot2','rot3'};
START = repmat([0.8 0.3],[10 1]);

% 10P2P
% names{1} = {'subj37','subj38','subj39','subj40','subj41'};
% names{2} = {'subj42','subj43','subj44'};
% groups = {'rot','mir'};
% path = 'Data/10P2P/';
% blocks = {'no_rot1','rot1','rot2','rot3'};
% START = repmat([0.8 0.3],[10 1]);

% no_P2P_15deg
% names{1} = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10'};
% groups = {'rot'};
% path = 'Data/no_P2P_15deg/';
% blocks = {'baseline','early','late','after'};
% START = repmat([0.6 0.25],[9 1]);
% START = [START; 0.6 0.3];

% no_P2P_15deg experiment 2
% names{1} = {'subj11','subj12','subj13','subj14','subj15','subj17','subj18','subj19','subj20','subj21'};
% groups = {'rot'};
% path = 'Data/no_P2P_15deg/';
% blocks = {'baseline','rot1','rot2','rot3','after'};
% START = repmat([0.8 0.3],[10 1]);

for i = 1:length(names)
    subjnames = names{i};
    Nsubj = length(subjnames);
    for subj = 1:Nsubj
        clear data
        disp(subjnames{subj});
        disp('    Loading Subject Data...');
        data = loadSubjData([path,subjnames{subj}],blocks,START(subj,:)); % load one chunk - loadSubjData(Subjname, {blocknames}); load in chunks of 5 based on target jumps

        % process data (smooth etc, rotate, get RT, etc.)
        disp('    Processing Data...')
        data = processData(data);

        % split into jump types - save this into a different data structure d
    %     disp('    Splitting data by jump type...')
    %         d{subj} = splitDatabyJump(data);

        % save data from this subject in a separate file
    %     fname = fullfile(['BimanualSkillData_S',num2str(subj)]);
    %     save(fname,'data')
        d.(groups{i}){subj} = data;
    end
end

% d = edit_initDir(d);
% save P2P d
disp('All Done')