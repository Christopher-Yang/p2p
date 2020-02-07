function data = processData(data)
% Compute initial reach direction. The output structure, data, is organized 
% as d.[group]{[subject number]}. The cell array for each subject contains 
% the following fields:
%
%   R: right-hand trajectory for each trial
%   C: cursor trajectory for each trial
%   Ntrials: number of trials
%   state: state of the trial, represented as a number
%       7: hold period before target presentation
%       3: time between target presentation and cursor leaving start circle
%       4: reach epoch
%       6: inter-trial interval (1 second)
%   time: movement timestamps
%   targetAbs: absolute target position
%   targetRel: target position relative to starting position
%   Cr: cursor trajectory rotated into common reference frame
%   init: time of movement initiation
%   iDir: time when initial reach direction is calculated
%   initDir: initial reach direction

for i=1:data.Ntrials
    % smooth trajectories
    data.Cr{i} = savgolayFilt(data.Cr{i}',3,7)';
    
    % compute velocities
    vel = diff(data.Cr{i})./(diff(data.time{i})/1000);
    
    % time of movement initiation
    data.init(i) = min(find(data.state{i}==4));
    
    % time when initial reach direction is calculated (150 ms after 
    % movement initiation)
    data.iDir(i) = data.init(i)+20; 
    
    % initial reach direction
    data.initDir(i) = atan2(vel(data.iDir(i),2),vel(data.iDir(i),1));
end