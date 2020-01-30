% loads raw data from all subjects into a series of .mat data files (one per subject)
%
% This example loads a single 'chunk' (block of 5 trials)

clear all
names{1} = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7'};
names{2} = {'subj13','subj18','subj19','subj21','subj22'};
names{3} = {'subj1'};
path = {'Data/denovo_2day/','Data/denovo_5day/','Data/denovo_10day/'};
blockNames = {'B5','B14','B29'};
START = repmat([0.6 0.25],[10 1]);

for i = 1:length(names)
    blocks = blockNames{i};
    subjnames = names{i};
    Nsubj = length(subjnames);
    for subj = 1:Nsubj
        clear data
        disp(subjnames{subj});
        disp('    Loading Subject Data...');
        data = loadSubjData([path{i},subjnames{subj}],blocks,START(subj,:)); % load one chunk - loadSubjData(Subjname, {blocknames}); load in chunks of 5 based on target jumps

        % process data (smooth etc, rotate, get RT, etc.)
        disp('    Processing Data...')
        data = processData(data);

        % split into jump types - save this into a different data structure d
    %     disp('    Splitting data by jump type...')
    %         d{subj} = splitDatabyJump(data);

        % save data from this subject in a separate file
    %     fname = fullfile(['BimanualSkillData_S',num2str(subj)]);
    %     save(fname,'data')
        switch i 
            case 1
                d.day2{subj} = data;
            case 2
                d.day5{subj} = data;
            case 3
                d.day10{subj} = data;
        end
    end
end

% d = edit_initDir(d);
% save P2P d
disp('All Done')