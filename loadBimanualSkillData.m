% This file is associated with the article "De novo learning versus
% adaptation of continuous control in a manual tracking task." Running
% loadBimanualSkillData.m generates Figure 2A.
% 
% Part 1 analyzes the point-to-point data while Part 2 plots the figure.
% Analysis and plotting code can be run independently of one another.
%
% To generate Figure 2A, run the entire script, which should take less than 
% 5 minutes. To avoid running the analysis multiple times, the variable "d"
% can be saved and loaded for plotting as needed.

%% PART 1: ANALYSIS
clear all

% set variables for analysis
path = 'Data/online/pilot_learning/'; % path to the data

d = loadSubjData(path);
data = processData(d);

% save P2P d % save data structure
disp('All Done')

%% PART 2: FIGURES
% load P2P % load data structure

plot_direction(d)