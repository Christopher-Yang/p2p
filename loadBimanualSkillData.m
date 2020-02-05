% This file is associated with the article "De novo learning versus
% adaptation of continuous control in a manual tracking task." Running
% loadBimanualSkillData.m generates Figure 2A.
% 
% Part 1 analyzes the point-to-point data while Part 2 plots the figure.
% Analysis and plotting code can be run independently of one another.
%
% To generate Figure 2A, run the entire script, which should take about 5 
% minutes. To avoid running the analysis multiple times, the variable "d"
% can be saved and loaded for plotting as needed.

%% PART 1: ANALYSIS
clear all

% set variables for analysis
names{1} = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'}; % rotation group
names{2} = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'}; % mirror-reversal group
groups = {'rot','mir'}; % names of groups
path = 'Data/rot_vs_mir/'; % path to the data
blocks = {'baseline','pert1','pert2','pert3'}; % names of the blocks
START = repmat([0.8 0.3],[10 1]); % initial position of the target in meters

for i = 1:length(names)
    subjnames = names{i};
    Nsubj = length(subjnames);
    for subj = 1:Nsubj
        clear data
        disp(subjnames{subj});
        
        % load data
        disp('    Loading Subject Data...');
        data = loadSubjData([path,subjnames{subj}],blocks,START(subj,:));

        % process data (smooth etc, rotate, get RT, etc.)
        disp('    Processing Data...')
        data = processData(data);

        d.(groups{i}){subj} = data; % store data into d
    end
end

% save P2P d % save data structure
disp('All Done')

%% PART 2: FIGURES
% load P2P % load data structure

plot_direction(d)