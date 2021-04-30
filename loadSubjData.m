function [data] = loadSubjData(subjname,blocknames,START)
% load a single subject's timed response target jump data

START_X = START(1);
START_Y = START(2);

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
        N{trial} = [L{trial}(:,1) R{trial}(:,2)]; % null space movements

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
        pert(trial) = tFile(j,5);
        
        tAng = atan2(targetRel(trial,2),targetRel(trial,1));
        targAng(trial) = tAng;
        
        % target position bin relative to mirror axis; bins numbered from 1
        % (closest)-4 (furthest)
%         if (tAng>=pi/8 && tAng<3*pi/8) || (tAng<-5*pi/8 && tAng>=-7*pi/8)
%             tBin = 1;
%         elseif (tAng>=0 && tAng<pi/8) || (tAng>=3*pi/8 && tAng<pi/2) || (tAng<-pi/2 && tAng>=-5*pi/8) || (tAng<-7*pi/8 && tAng>=-pi)
%             tBin = 2;
%         elseif (tAng>=pi/2 && tAng<5*pi/8) || (tAng>=7*pi/8 && tAng<pi) || (tAng<0 && tAng>=-pi/8) || (tAng<-3*pi/8 && tAng >=-pi/2)
%             tBin = 3;
%         elseif (tAng>=5*pi/8 && tAng<7*pi/8) || (tAng<-pi/8 && tAng>=-3*pi/8)
%             tBin = 4;
%         end
        
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
% data.targAng = atan2(targetRel(:,2),targetRel(:,1));
data.targAng = targAng;
data.targIndex = data.targAng*(8/(2*pi));
data.targDist = sqrt(sum(targetRel(:,1:2)'.^2));
data.targBin = targBin;

% store all info in data structure 'data'
data.L = L;
data.R = R;
data.C = C;
data.N = N;

data.Ntrials = size(targetRel,1);
data.tFile = tFileFull;
data.pert = pert;

data.state = state;
data.time = time;
data.itargonset = itargonset;
data.imoveonset = imoveonset;

data.blocknames = blocknames;

% placeholders - these will be computed later
% data.RT = d0;
% data.iDir = d0;

data.targetAbs = targetAbs;
data.targetRel = targetRel;
data.start = start;

% rotate data into common coordinate frame - start at (0,0), target at
% (0,.12)
for j=1:data.Ntrials % iterate through all trials
    theta(j) = atan2(data.targetRel(j,2),data.targetRel(j,1))-pi/2;
    R = [cos(theta(j)) sin(theta(j)); -sin(theta(j)) cos(theta(j))];
    
    data.Cr{j} = (R*(data.C{j}'-repmat(start(j,:),size(data.C{j},1),1)'))';
    data.Nr{j} = (R*(data.N{j}'))';
    
    targetMir = [0 1; 1 0]* data.targetRel(j,:)';
    thetaMir(j) = atan2(targetMir(2),targetMir(1))-pi/2;
    R = [cos(thetaMir(j)) sin(thetaMir(j)); -sin(thetaMir(j)) cos(thetaMir(j))];
    
    data.Cr_mir{j} = (R*(data.C{j}'-repmat(start(j,:),size(data.C{j},1),1)'))';
    data.Nr_mir{j} = (R*(data.N{j}'))';
end
