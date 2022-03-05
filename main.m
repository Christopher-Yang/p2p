% This file is associated with the article "Emergence of habitual control
% in a novel motor skill over multiple days of practice". Running
% main.m generates plots for Figures 2, 4, and S2.
% 
% Part 1 loads raw data from the point-to-point task and performs initial
% data analysis, storing the data in a structure called "data". This
% structure is organized as follows: data.(group){subject}. See processData
% for more details on subfields in "data". Part 2 performs follow-up data 
% analysis and plots figures.
% 
% All analyses can be performed by running main.m

%% PART 1: INITIAL DATA  ANALYSIS

clear all

% path to data for each group
path = {'Data/denovo_2day/','Data/denovo_5day/','Data/denovo_10day/'};

% select names of subjects to be analyzed; names{1}: 2-day group; names{2}:
% 5-day group; names{3}: 10-day group
names{1} = {'subj1','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10','subj11','subj12','subj13','subj14'};
names{2} = {'subj13','subj15','subj17','subj18','subj19','subj21','subj22','subj23','subj24','subj25','subj26','subj27','subj28','subj29'};
names{3} = {'subj1','subj2','subj3','subj4','subj5'};

% select blocks to be analyzed for each group
blockNames{1} = {'B1_baseline','B2','B3','B4','B5','B6_habit'};
blockNames{2} = {'B1_baseline','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15_habit'};
blockNames{3} = {'B1_baseline','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30_habit'};

% starting position of the target: [x y]
START = [0.6 0.25];

% loop to analyze data
for i = 1:length(names) % loop over groups
    
    blocks = blockNames{i}; % blocks in current group
    subjnames = names{i}; % subjects in current group
    Nsubj = length(subjnames); % number of subjects in current group
    for subj = 1:Nsubj % loop over subjects
        
        clear d
        disp(subjnames{subj});
        
        % load raw data from data files
        disp('    Loading Subject Data...');
        d = loadSubjData([path{i},subjnames{subj}],blocks); 

        % initial data analysis (smooth trajectories, rotate, get RT, etc.)
        disp('    Processing Data...')
        d = processData(d);
        
        % store data in "data"
        switch i 
            case 1
                data.day2{subj} = d;
            case 2
                data.day5{subj} = d;
            case 3
                data.day10{subj} = d;
        end
    end
end

disp('All Done')

%% PART 2: FOLLOW-UP DATA ANALYSIS AND PLOTTING FIGURES

% plot Figure 2A, 4A, and Supplementary Figure 4A
plot_traj(data)

% plot Figure 2B-C
plot_direction(data)

% plot Figure 4B and Supplementary Figure 3A
plot_heatmap(data)

% plot Figure 4C-E and Supplementary Figure 3B
plot_flip(data)

% plot Supplementary Figure 1
plot_kinematics(data)

% plot Supplementary Figure 3C
% 
% set loadAccuracy = 1 if you want to use precomputed accuracy matrix
% 
% to compute accuracy matrix from scratch, set loadAccuracy = 0, which will
% take about 15 mins to run
loadAccuracy = 1;
modelRecovery(loadAccuracy)

% plot Supplementary Figure 4B
plot_away(data)
