function data = loadSubjData(subjname,blocknames,START)
% Load a single subject's point-to-point data for a given block 
% (blocknames) and starting position of the target (START). The output 
% variable "data" contains cell arrays to be used for analysis in
% processData(). See that function for details on the meaning of each
% array.

START_X = START(1);
START_Y = START(2);

Nblocks = length(blocknames);
trial = 1;
tFileFull = [];
for blk=1:Nblocks % iterate over blocks
    path = [subjname,'/',blocknames{blk}];
    tFile = dlmread([path,'/tFile.tgt'],' ',0,0);
    fnames = dir(path);
    Ntrials = size(tFile,1);
    for j=1:Ntrials % iterate over trials
        d = dlmread([path,'/',fnames(j+2).name],' ',6,0);

        R{trial} = d(:,3:4); % right hand X and Y
        C{trial} = d(:,5:6);% cursor X and Y

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
        
        state{trial} = d(:,7); % trial 'state' at each time point
        time{trial} = d(:,9); % time during trial
        
        trial = trial+1;
    end
    tFileFull = [tFileFull; tFile(:,1:5)]; % copy of trial table
end

% store all info in data structure 'data'
data.R = R;
data.C = C;
data.Ntrials = size(targetRel,1);
data.state = state;
data.time = time;
data.targetAbs = targetAbs;
data.targetRel = targetRel;

% rotate data into common coordinate frame - start at (0,0), target at
% (0,.12)
for j=1:data.Ntrials % iterate through all trials
    theta(j) = atan2(data.targetRel(j,2),data.targetRel(j,1))-pi/2;
    R = [cos(theta(j)) sin(theta(j)); -sin(theta(j)) cos(theta(j))];
    
    % rotated cursor trajectory
    data.Cr{j} = (R*(data.C{j}'-repmat(start(j,:),size(data.C{j},1),1)'))';
end