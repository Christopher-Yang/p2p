% processData() analyzes the raw data from single subjects produced by 
% loadSubjData(). The output structure "data" contains cell arrays where
% each element of the array corresponds to data from one trial. The fields
% in "data" correspond to the following:
% 
%   Cr: cursor trajectories rotated into a common reference frame
%   L: left hand trajectory
%   R: right hand trajectory
%   C: unrotated cursor trajectory
%   Ntrials: total number of trials across all analyzed blocks
%   state: state of the trial, represented as a number
%       7: hold period before target presentation
%       3: time between target presentation and cursor leaving start circle
%       4: reach epoch
%       6: inter-trial interval (1 second)
%   time: time data was collected
%   targAng: relative angle between starting and ending targets
%   targetAbs: absolute position of the ending targets
%   targetRel: relative position between consecutive targets where starting
%              position is defined as (0,0)
%   start: absolute position of the starting targets
%   go: time of go cue
%   init: time of movement initiation (when cursor left starting target)
%   end: time of movement end (cursor stationary in end target)
%   movtime: movement time
%   tanVel: tangential velocity of the cursor
%   pkVel: peak velocity of the cursor
%   initVel: initial velocity of the cursor
%   RT: reaction time:
%   iDir_x: time at which initial horizontal reach direction was assessed
%   away: boolean for whether initial horizontal reach on a trial was aimed
%         away from the target
%   pathlength: path length
%   iDir: time at which overall initial reach direction was assessed
%   initDir: initial reach direction for rotated cursor trajectories
%   initDir_noRot: initial reach direction for unrotated cursor
%                  trajectories

function data = processData(data)
delay = 20; % time after reach initiation to compute initial reach direction
fs = 130; % sampling rate

for i=1:data.Ntrials % loop over trials
    
    % compute various time markers
    time = (data.time{i} - data.time{i}(1))/1000; % time at which data were collected
    data.time{i} = time;
    data.go(i) = find(data.state{i}==3,1); % time of go cue
    data.init(i) = find(data.state{i}==4,1); % time of movement initiation
    data.end(i) = find(data.state{i}==4,1,'last'); % time of movement end
    data.movtime(i) = time(data.end(i))-time(data.init(i)); % movement time

    % smooth trajectories
    data.Cr{i} = sgolayfilt(data.Cr{i},3,7); % rotated cursor trajectory
    data.C{i} = sgolayfilt(data.C{i},3,7); % unrotated cursor trajectory
    
    % resample trajectories to make velocity calculation more consistent
    C = resample(data.C{i},time);
    Cr = resample(data.Cr{i},time);
        
    % compute velocities
    vel_C = diff(C)*fs; % cursor velocity unrotated
    data.tanVel{i} = sgolayfilt(sqrt(sum(vel_C.^2,2)),3,7); % tangential velocity
    vel_C = sgolayfilt(vel_C,3,7); % filter velocity
    vel_Cr = sgolayfilt(diff(Cr)*fs,3,7); % cursor velocity rotated and filtered
    data.pkVel(i) = max(data.tanVel{i}); % peak tangential velocity
    data.initVel(i) = data.tanVel{i}(data.init(i)+delay); % initial velocity
    
    % calculate reaction time
    idx2 = data.tanVel{i} > 0.1; % find timepoints where velocity exceed 0.1 m/s
    idx2(1:data.go(i)) = 0; % ignore all timepoints before go cue
    idx = find(idx2 == 1,1); % find first timepoint after go cue where velocity exceeded 0.1 m/s
    if isempty(idx) % if velocity didn't exceed 0.1 m/s, then set reaciton time to NaN
        data.RT(i) = NaN;
    else % else, compute reaction time
        data.RT(i) = time(idx) - time(data.go(i)) - 0.1; % subtract off 100 ms to account for Kinereach delay
    end
    
    % find timepoints where cursor deviated horizontally from the starting
    % target by more than 1 cm 
    threshold_x = abs(data.C{i}(:,1)-data.start(i,1)) > 0.01;
    threshold_x(1:data.init(i)) = 0; % ignore timepoints before movement initiation
    
    % if cursor doesn't deviate horizontally from starting target, then
    % don't analyze data
    if sum(threshold_x) == 0
        data.away(i) = NaN;
        
    else
        
        % find first timepoint where cursor deviated from starting target
        init_x = find(threshold_x==1,1);
        
        % if timepoint is too late during the trial, then analysis cannot
        % be performed because the analysis time with delay comes after the
        % trial ends
        if init_x+delay > size(vel_C,1)
            data.away(i) = NaN;
            
        % else, perform analysis
        else
            
            % find initial horizontal velocity
            data.iDir_x(i) = init_x+delay; % analysis time
            initVel_x = vel_C(data.iDir_x(i),1);
            
            % determine whether initial horizontal velocity is aimed
            % towards or away from target
            if data.targetRel(i,1) < 0 && initVel_x > 0 % away
                data.away(i) = 1;
            elseif data.targetRel(i,1) > 0 && initVel_x < 0 % away
                data.away(i) = 1;
            else % toward
                data.away(i) = 0;
            end
        end
    end
        
    % compute path length
    dpath = diff(data.Cr{i}(1:data.end(i),:));
    dL = sqrt(sum(dpath.^2,2));
    data.pathlength(i) = sum(dL);
    
    % initial reach direction for rotated trajectories
    iDir = data.init(i)+delay; % 150 ms (20 time steps) after movement initiation
    initDir = atan2(vel_Cr(iDir,2),vel_Cr(iDir,1));
    while initDir >= pi % wrap directions of pi to -pi
        initDir = initDir-2*pi;
    end
    
    % initial reach direction for unrotated trajectories
    initDir_noRot = atan2(vel_C(iDir,2),vel_C(iDir,1));
    while initDir_noRot >= pi % wrap directions of pi to -pi
        initDir_noRot = initDir_noRot-2*pi;
    end
    
    % store reach directions
    data.iDir(i) = iDir; % analysis time
    data.initDir(i) = initDir; % initial reach direction for rotated trajectory
    data.initDir_noRot(i) = initDir_noRot; % initial reach direction for unrotated trajectory
end