function data = processData(data)
% analyze cursor data for path length, RT, etc smooth trajectories

delay = 20; % time after reach initiation to compute initial reach direction

for i=1:data.Ntrials
    % smooth trajectories
    data.Cr{i} = sgolayfilt(data.Cr{i},3,7);
    data.C{i} = sgolayfilt(data.C{i},3,7);
    
    fs = 130;
    time = (data.time{i} - data.time{i}(1))/1000;
    C = resample(data.C{i},time);
    Cr = resample(data.Cr{i},time);
        
    % compute velocities
    vel_C = diff(C)*fs; % cursor velocity unrotated
    vel_Cr = diff(Cr)*fs; % cursor velocity rotated
    
    data.tanVel{i} = sqrt(sum(vel_Cr.^2,2)); % tangential velocity
    data.pkVel(i) = max(data.tanVel{i}); % peak tangential velocity
    velFilt = sgolayfilt(data.tanVel{i},3,7); % smooth velocity trajectory
    vel_C = sgolayfilt(vel_C,3,7);
    data.vel_C = vel_C;
    
    % some of the differences in time are 0, leading to Inf velocities;
    % replace with NaN
%     velFilt(velFilt == Inf) = NaN;
%     velFilt(velFilt == -Inf) = NaN;
    data.velFilt{i} = velFilt;
    
    data.go(i) = min(find(data.state{i}==3)); % time of go cue
    data.init(i) = min(find(data.state{i}==4)); % time of movement initiation
    data.end(i) = max(find(data.state{i}==4)); % time of movement end
    
    data.initVel(i) = data.tanVel{i}(data.init(i)+delay);
    data.initVel_filt(i) = velFilt(data.init(i)+delay);
    
    %% identify initial reach direction in x- and y-axes 
    if i == 1
        threshold_x = abs(data.C{i}(:,1)-0.6) > 0.01;
    else
        threshold_x = abs(data.C{i}(:,1)-data.targetAbs(i-1,1)) > 0.01;
    end
    
    threshold_x(1:data.init(i)) = 0;
    if sum(threshold_x) == 0
        data.init_x(i) = NaN;
        data.initVel_x(i) = NaN;
        data.incorrectReach_x(i) = NaN;
    else
        data.init_x(i) = find(threshold_x==1,1,'first');
        if data.init_x(i)+delay > size(vel_C,1)
            data.init_x(i) = NaN;
            data.initVel_x(i) = NaN;
            data.incorrectReach_x(i) = NaN;
        else
            data.initVel_x(i) = vel_C(data.init_x(i)+delay,1);
            
            % incorrectReach_x = 1 if there is a habit, incorrectReach_x = 0 if there's no
            % habit
            if data.targetRel(i,1) < 0 && data.initVel_x(i) > 0
                data.incorrectReach_x(i) = 1;
            elseif data.targetRel(i,1) > 0 && data.initVel_x(i) < 0
                data.incorrectReach_x(i) = 1;
            else
                data.incorrectReach_x(i) = 0;
            end
        end
    end
    
    data.RT(i) = data.time{i}(data.init(i))-data.time{i}(data.go(i)) - 100; % reaction time
    if isnan(data.init_x(i))
        data.RT_x(i) = NaN;
    else
        RT_x = data.time{i}(data.init_x(i))-data.time{i}(data.go(i)) - 100; % reaction time
        if RT_x < 0
            data.RT_x(i) = NaN;
        else
            data.RT_x(i) = RT_x;
        end
    end
        
    data.movtime(i) = data.time{i}(data.end(i))-data.time{i}(data.init(i)); % movement time
    
    % compute path length
    dpath = diff(data.Cr{i}(1:data.end(i),:));
    dL = sqrt(sum(dpath.^2,2));
    data.pathlength(i) = sum(dL);
    
    % compute initial reach direction
    iDir = data.init(i)+delay; % 100 ms (13 time steps) after initiation
    initDir = atan2(vel_Cr(iDir,2),vel_Cr(iDir,1));
    initDir = initDir-pi/2;
    while initDir >= pi
        initDir = initDir-2*pi;
    end
    while initDir < -pi
        initDir = initDir+2*pi;
    end
        
    initDir_noRot = atan2(vel_C(iDir,2),vel_C(iDir,1));
    while initDir_noRot >= pi
        initDir_noRot = initDir_noRot-2*pi;
    end
    while initDir_noRot < -pi
        initDir_noRot = initDir_noRot+2*pi;
    end
    
    data.iDir(i) = iDir;
    data.initDir(i) = initDir;
    data.initDir_noRot(i) = initDir_noRot;
end