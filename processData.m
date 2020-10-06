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

disp('Loading...');

Nsubj = length(data);
Ntrials = length(data{1}.frameDrops);
fields = {'time','C','Cr'};

for k = 1:Nsubj
    disp(['   Subj ' num2str(k)]);
    for i = 1:Ntrials
        for j = 1:length(fields)
            data{k}.(fields{j}) = interpolate_data(data{k}.(fields{j}));
        end
        
%         Cr = data{k}.Cr{i};
        
        % smooth trajectories
%         data{k}.Cr{i} = savgolayFilt(data{k}.Cr{i}',3,7)';

        % compute velocities
        vel = diff(data{k}.Cr{i})./(diff(data{k}.time{i})/1000);

        % time of movement initiation
        init = find(data{k}.begin{i} == 1, 1);

        % time when initial reach direction is calculated (150 ms after 
        % movement initiation)
        data{k}.iDir(i) = init+2;

        % initial reach direction
        data{k}.initDir(i) = atan2(vel(data{k}.iDir(i),2),vel(data{k}.iDir(i),1));
    end
end
end

% interpolates missing data (note that this method may interpolate negative
% times instead of 0 at the start of the trial)
function output = interpolate_data(data)
    Ntrials = size(data,2);
    for k = 1:Ntrials
        Ncols = size(data{k},2);
        for i = 1:Ncols
            d = data{k}(:,i);
            while sum(isnan(d))
                nans = find(isnan(d),1);
                numbers = find(~isnan(d));
                before = numbers(numbers < nans);
                after = numbers(numbers > nans);
                if isempty(before)
                    d(1) = d(2) - (d(3)-d(2));
                elseif isempty(after)
                    d(end) = d(end-1) - (d(end-2)-d(end-1));
                else
                    before = before(end);
                    after = after(1);
                    d(before:after) = linspace(d(before),d(after),after-before+1);
                end
            end
            data{k}(:,i) = d;
        end
    end
    output = data;
end