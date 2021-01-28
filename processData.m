function data = processData(data)
% analyze cursor data for path length, RT, etc smooth trajectories

for i=1:data.Ntrials
    % smooth trajectories
    data.Cr{i} = savgolayFilt(data.Cr{i}',3,7)';
    data.Nr{i} = savgolayFilt(data.Nr{i}',3,7)';
    data.C{i} = savgolayFilt(data.C{i}',3,7)';
    
    % compute velocities
    vel_Cr = diff(data.Cr{i})./(diff(data.time{i})/1000);
    data.tanVel{i} = sqrt(sum(vel_Cr.^2,2)); % tangential velocity
    data.pkVel(i) = max(data.tanVel{i}); % peak tangential velocity
    velFilt = sgolayfilt(data.tanVel{i},3,9); % smooth velocity trajectory
    
    % some of the differences in time are 0, leading to Inf velocities;
    % replace with NaN
    velFilt(velFilt == Inf) = NaN;
    velFilt(velFilt == -Inf) = NaN;
    data.velFilt{i} = velFilt;
    
    vel_C = diff(data.C{i})./(diff(data.time{i})/1000);
    
    data.go(i) = min(find(data.state{i}==3)); % time of go cue
    data.init(i) = min(find(data.state{i}==4)); % time of movement initiation
    data.end(i) = max(find(data.state{i}==4)); % time of movement end
    
    %% identify initial reach direction in x- and y-axes 
    if i == 1
        threshold_x = abs(data.C{i}(:,1)-0.6) > 0.01;
        threshold_y = abs(data.C{i}(:,2)-0.25) > 0.01;
    else
        threshold_x = abs(data.C{i}(:,1)-data.targetAbs(i-1,1)) > 0.01;
        threshold_y = abs(data.C{i}(:,2)-data.targetAbs(i-1,2)) > 0.01;
    end
    
    if sum(threshold_y) == 0
        data.init_y(i) = NaN;
        data.initDir_y(i) = NaN;
        data.incorrectReach_y(i) = NaN;
    else
        data.init_y(i) = min(find(threshold_y==1));
        if data.init_y(i)+20 > size(vel_C,1)
            data.init_y(i) = NaN;
            data.initDir_y(i) = NaN;
            data.incorrectReach_y(i) = NaN;
        else
            data.initDir_y(i) = vel_C(data.init_y(i)+20,2);
            
            % incorrectReach_y = 1 if there is a habit, incorrectReach_y = 0 if there's no
            % habit
            if data.targetRel(i,2) < 0 && data.initDir_y(i) > 0
                data.incorrectReach_y(i) = 1;
            elseif data.targetRel(i,2) > 0 && data.initDir_y(i) < 0
                data.incorrectReach_y(i) = 1;
            else
                data.incorrectReach_y(i) = 0;
            end
        end
    end
    
    if sum(threshold_x) == 0
        data.init_x(i) = NaN;
        data.initDir_x(i) = NaN;
        data.incorrectReach_x(i) = NaN;
    else
        data.init_x(i) = min(find(threshold_x==1));
        if data.init_x(i)+20 > size(vel_C,1)
            data.init_x(i) = NaN;
            data.initDir_x(i) = NaN;
            data.incorrectReach_x(i) = NaN;
        else
            data.initDir_x(i) = vel_C(data.init_x(i)+20,1);
            
            % incorrectReach_x = 1 if there is a habit, incorrectReach_x = 0 if there's no
            % habit
            if data.targetRel(i,1) < 0 && data.initDir_x(i) > 0
                data.incorrectReach_x(i) = 1;
            elseif data.targetRel(i,1) > 0 && data.initDir_x(i) < 0
                data.incorrectReach_x(i) = 1;
            else
                data.incorrectReach_x(i) = 0;
            end
        end
    end
    
    %%
    % should technically use data.init_x(i)+20 and similar for init_y, but
    % on some trials the initiation time is close to the end of the trial,
    % meaning the reach ends before adding on 150 ms 
    if ~isnan(data.init_x(i)) && ~isnan(data.init_y(i))
        data.lag(i) = data.time{i}(data.init_x(i))-data.time{i}(data.init_y(i));
    else
        data.lag(i) = NaN;
    end
    
    data.RT(i) = data.time{i}(data.init(i))-data.time{i}(data.go(i)); % reaction time
    data.movtime(i) = data.time{i}(data.end(i))-data.time{i}(data.init(i)); % movement time
    
    % compute path length
    dpath = diff(data.Cr{i}(1:data.end(i),:));
    dL = sqrt(sum(dpath.^2,2));
    data.pathlength(i) = sum(dL);
    
    % compute path length in null space
    dpath_null = diff(data.Nr{i}(1:data.end(i),:));
    dL = sqrt(sum(dpath_null.^2,2));
    data.pathlength_null(i) = sum(dL);
    
    % compute initial reach direction
    data.iDir(i) = data.init(i)+20; % 100 ms (13 time steps) after initiation
    data.initDir(i) = atan2(vel_Cr(data.iDir(i),2),vel_Cr(data.iDir(i),1));
    data.initDir(i) = data.initDir(i)-pi/2;
    while data.initDir(i) > pi
        data.initDir(i) = data.initDir(i)-2*pi;
    end
    while data.initDir(i) < -pi
        data.initDir(i) = data.initDir(i)+2*pi;
    end
    
    data.initDir_noRot(i) = atan2(vel_C(data.iDir(i),2),vel_C(data.iDir(i),1));
    data.initDir_noRot(i) = data.initDir_noRot(i);
    while data.initDir_noRot(i) > pi
        data.initDir_noRot(i) = data.initDir_noRot(i)-2*pi;
    end
    while data.initDir_noRot(i) < -pi
        data.initDir_noRot(i) = data.initDir_noRot(i)+2*pi;
    end
    
end