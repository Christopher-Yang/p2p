% This file is associated with the article "Emergence of habitual control
% in a novel motor skill over multiple days of practice". Running
% main.m generates plots for Figures 2 and 4 as well as Supplementary
% Figures 1, 3 and 4.
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
path = 'Data/CorsiMentalRotation_2day/';

% select names of subjects to be analyzed; names{1}: 2-day group; names{2}:
% 5-day group; names{3}: 10-day group
subj_name = {'subj03','subj04','subj05','subj06','subj07','subj08','subj09','subj10','subj11','subj12','subj13','subj14','subj15','subj16','subj17','subj18'};

% select blocks to be analyzed for each group
block_name = {'d1_prac_avg1', 'd1_prac_bi1', 'd1_prac_bi2', 'd1_prac_bi3', 'd1_prac_bi4', 'd1_prac_bi5', 'd1_prac_bi6', 'd2_prac_bi7'};

% starting position of the target: [x y]
START = [0.6 0.25];
Nsubj = length(subj_name);

for subj = 1:Nsubj
    disp(subj_name{subj})
    
    % load raw data from data files
    disp('    Loading Subject Data...');
    d = loadSubjData([path subj_name{subj}], block_name);
    
    % initial data analysis (smooth trajectories, rotate, get RT, etc.)
    disp('    Processing Data...')
    d = processData(d);
    
    data{subj} = d;
end

disp('All Done')

%% PART 2: FOLLOW-UP DATA ANALYSIS AND PLOTTING FIGURES

% plot Figure 2A, 4A, and Supplementary Figure 2A
plot_traj(data)

% plot Figure 2B-C
plot_direction(data)

% plot Figure 2D-G
plot_kinematics(data)

% plot Figure 4B and Supplementary Figure 1A
plot_heatmap(data)
