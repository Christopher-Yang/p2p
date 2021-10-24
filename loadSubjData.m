function [data] = loadSubjData(subjname,blocknames)
% load a single subject's timed response target jump data

START_X = 0.6;
START_Y = 0.25;

Nblocks = length(blocknames);
trial = 1;
tFileFull = [];
for blk=1:Nblocks
    path = [subjname,'/',blocknames{blk}];
    tFile = dlmread([path,'/tFile.tgt'],' ',0,0);
    fnames = dir(path);
    Ntrials = size(tFile,1);
    for j=1:Ntrials
        d = dlmread([path,'/',fnames(j+2).name],' ',6,0);

        L{trial} = d(:,1:2); % left hand X and Y
        R{trial} = d(:,3:4); % right hand X and Y
        C{trial} = d(:,5:6); % cursor X and Y

        % absolute target location
        targetAbs(trial,1) = tFile(j,2)+START_X;
        targetAbs(trial,2) = tFile(j,3)+START_Y;
        
        % determine relative target location
        if(j>1)
            start(trial,:) = targetAbs(trial-1,:);    
        else
            start(trial,:) = [START_X START_Y];
        end
        targetRel(trial,:) = targetAbs(trial,:)-start(trial,:);
        
        tAng = atan2(targetRel(trial,2),targetRel(trial,1));
        targAng(trial) = tAng;
        
        itarg = find(d(:,7)==3); % time of target movement
        if(isempty(itarg))
            itargonset = 1;
        else
            itargonset(trial) = min(itarg);
        end        
        
        imov = find(d(:,7)==4); % time of movement onset
        if(isempty(imov))
            imoveonset = 1;
        else
            imoveonset(trial) = min(imov);
        end        
        state{trial} = d(:,7); % trial 'state' at each time point
        time{trial} = d(:,9); % time during trial
        
        trial = trial+1;
    end
    tFileFull = [tFileFull; tFile(:,1:5)]; % copy of trial table
end

% compute target angle

% store all info in data structure 'data'
data.L = L;
data.R = R;
data.C = C;

data.Ntrials = size(targetRel,1);
data.tFile = tFileFull;

data.state = state;
data.time = time;
data.itargonset = itargonset;
data.imoveonset = imoveonset;

data.targAng = targAng;
data.targetAbs = targetAbs;
data.targetRel = targetRel;
data.start = start;

% rotate data into common coordinate frame - start at (0,0), target at
% (0,.12)
for j=1:data.Ntrials % iterate through all trials
    theta = atan2(data.targetRel(j,2),data.targetRel(j,1))-pi/2;
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    
    data.Cr{j} = (R*(data.C{j}'-repmat(start(j,:),size(data.C{j},1),1)'))';
end
