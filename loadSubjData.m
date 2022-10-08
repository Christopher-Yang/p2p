% loadSubjData() loads raw data from data files
% 
%   subjname: names of subjects
%   blocknames: names of blocks
%
% The data structure "data" contains fields to be used for analysis in
% processData(). See that function for a full description of "data"

function data = loadSubjData(subjname,blocknames)

% starting position of the target on trial 1
START_X = 0.6;
START_Y = 0.25;

Nblocks = length(blocknames); % total number of blocks
trial = 1; % counter for storing data

for blk=1:Nblocks % loop over blocks
    path = [subjname,'/',blocknames{blk}]; % path to data
    tFile = dlmread([path,'/tFile.tgt'],' ',0,0); % loads file used to set target positions in experiment
    fnames = dir(path); % names of files in the file path
    Ntrials = size(tFile,1); % number of trials
    for j=1:Ntrials % loop over trials
        d = dlmread([path,'/',fnames(j+2).name],' ',8,0); % read data files

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
        targetRel(trial,:) = targetAbs(trial,:)-start(trial,:); % relative target position
        tAng = atan2(targetRel(trial,2),targetRel(trial,1));
        targAng(trial) = tAng; % relative target angle
        
        % target position bin relative to y-axis; bins numbered from 1
        % (closest)-4 (furthest)
        if (tAng>=3*pi/8 && tAng<5*pi/8) || (tAng>=-5*pi/8 && tAng<-3*pi/8)
            tBin = 1;
        elseif (tAng>=2*pi/8 && tAng<3*pi/8) || (tAng>=5*pi/8 && tAng<6*pi/8) || (tAng>=-3*pi/8 && tAng<-2*pi/8) || (tAng>=-6*pi/8 && tAng<-5*pi/8)
            tBin = 2;
        elseif (tAng>=pi/8 && tAng<2*pi/8) || (tAng>=6*pi/8 && tAng<7*pi/8) || (tAng>=-2*pi/8 && tAng<-pi/8) || (tAng>=-7*pi/8 && tAng<-6*pi/8)
            tBin = 3;
        elseif (tAng>=7*pi/8 && tAng<pi) || (tAng>=-pi && tAng<-7*pi/8) || (tAng>=-pi/8 && tAng<pi/8)
            tBin = 4;
        end
        targBin(trial) = tBin;
        
        state{trial} = d(:,7); % trial 'state' at each time point
        time{trial} = d(:,9); % time during trial
        
        trial = trial+1; % trial counter
    end
end

% rotate data into common coordinate frame - start at (0,0), target at
% (0.12,0)
Ntrials = size(targetRel,1);
for j=1:Ntrials % iterate through all trials
    theta = targAng(j); % angle to target
    rotMat = [cos(theta) sin(theta); -sin(theta) cos(theta)]; % make rotation matrix
    data.Cr{j} = (rotMat*(C{j}'-repmat(start(j,:),size(C{j},1),1)'))'; % rotate cursor trajectory
end

% store all info in data structure 'data'
data.L = L;
data.R = R;
data.C = C;
data.Ntrials = Ntrials;
data.state = state;
data.time = time;
data.targAng = targAng;
data.targetAbs = targetAbs;
data.targetRel = targetRel;
data.targBin = targBin;
data.start = start;
