function data = processData(data)
% analyze cursor data for path length, RT, etc smooth trajectories

for i=1:data.Ntrials
    % smooth trajectories
    data.Cr{i} = savgolayFilt(data.Cr{i}',3,7)';
    data.Nr{i} = savgolayFilt(data.Nr{i}',3,7)';
    
    % compute velocities
    vel = diff(data.Cr{i})./(diff(data.time{i})/1000);
    data.tanVel{i} = sqrt(sum(vel.^2,2));
    data.pkVel(i) = max(data.tanVel{i});
    
    data.go(i) = min(find(data.state{i}==3)); % time of go cue
    data.init(i) = min(find(data.state{i}==4)); % time of movement initiation
    data.end(i) = max(find(data.state{i}==4)); % time of movement end
    
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
    data.initDir(i) = atan2(vel(data.iDir(i),2),vel(data.iDir(i),1));
    data.initDir(i) = data.initDir(i)-pi/2;
    while data.initDir(i) > pi
        data.initDir(i) = data.initDir(i)-2*pi;
    end
    while data.initDir(i) < -pi
        data.initDir(i) = data.initDir(i)+2*pi;
    end
    
    % compute velocity and position at time of jump (could be useful to
    % know)
    data.vel{i} = vel;
    data.vel_at_jump(i,:) = vel(data.ipertonset(i),:)';
    data.pos_at_jump(i,:) = data.Cr{i}(data.ipertonset(i),:);
end