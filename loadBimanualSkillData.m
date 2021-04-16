% loads raw data from all subjects into a series of .mat data files (one per subject)
%
% This example loads a single 'chunk' (block of 5 trials)

clear all
names{1} = {'subj1','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10','subj11','subj12','subj13','subj14'};
names{2} = {'subj13','subj15','subj17','subj18','subj19','subj21','subj22','subj23','subj24','subj25','subj26','subj27','subj28','subj29'};
names{3} = {'subj1','subj2','subj3','subj4','subj5'};
path = {'Data/denovo_2day/','Data/denovo_5day/','Data/denovo_10day/'};
% blockNames{1} = {'B1_baseline','B2','B5','B6_habit'};
% blockNames{2} = {'B1_baseline','B2','B5','B14','B15_habit'};
% blockNames{3} = {'B1_baseline','B2','B5','B14','B29','B30_habit'};

% use for plot_tortuosity.m and plot_traj.m
blockNames{1} = {'B1_baseline','B2','B3','B4','B5','B6_habit'};
blockNames{2} = {'B1_baseline','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15_habit'};
blockNames{3} = {'B1_baseline','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30_habit'};
START = [0.6 0.25];

for i = 1:length(names)
    blocks = blockNames{i};
    subjnames = names{i};
    Nsubj = length(subjnames);
    for subj = 1:Nsubj
        clear data
        disp(subjnames{subj});
        disp('    Loading Subject Data...');
        data = loadSubjData([path{i},subjnames{subj}],blocks,START); % load one chunk - loadSubjData(Subjname, {blocknames}); load in chunks of 5 based on target jumps

        % process data (smooth etc, rotate, get RT, etc.)
        disp('    Processing Data...')
        data = processData(data);

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

% save('variables/d.mat','d') % can't save variables larger than 2 GB
disp('All Done')

%% plot data

plot_traj(d);
plot_kinematics(d);
plot_habit(d);
plot_direction(d);

% load('variables/vmFit')
% plot_direction(d,vmFit)