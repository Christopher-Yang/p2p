function data = splitDatabyJump(data)
frameSize = 260; % max number of samples to keep around perturbation. 260 samples = 2s
% align all data to perturbation time
% figure out how many repetitions
Nblocks = data.Ntrials/150; % number of blocks (10 trials per block)

eliminate_bad_trials = 0; % eliminate bad trials (1 or 0)

% figure out size of data structure
post = data.end - data.ipertonset;
pre = data.ipertonset - data.init;
maxPost = max(max(post),frameSize);
maxPre = max(max(pre),frameSize);

% fill with data pre- and post-perturbation
for i=1:data.Ntrials
    % cursor
    CrX_post(i,:) = [data.Cr{i}(data.ipertonset(i)+1:data.end(i),1)' NaN*ones(1,maxPost-post(i))];   
    CrY_post(i,:) = [data.Cr{i}(data.ipertonset(i)+1:data.end(i),2)' NaN*ones(1,maxPost-post(i))];
    CrX_pre(i,:) = [NaN*ones(1,maxPre-pre(i)) data.Cr{i}(data.init(i):data.ipertonset(i),1)'];
    CrY_pre(i,:) = [NaN*ones(1,maxPre-pre(i)) data.Cr{i}(data.init(i):data.ipertonset(i),2)'];
    
    % null-space
    NrX_post(i,:) = [data.Nr{i}(data.ipertonset(i)+1:data.end(i),1)' NaN*ones(1,maxPost-post(i))];
    NrY_post(i,:) = [data.Nr{i}(data.ipertonset(i)+1:data.end(i),2)' NaN*ones(1,maxPost-post(i))];
    NrX_pre(i,:) = [NaN*ones(1,maxPre-pre(i)) data.Nr{i}(data.init(i):data.ipertonset(i),1)'];
    NrY_pre(i,:) = [NaN*ones(1,maxPre-pre(i)) data.Nr{i}(data.init(i):data.ipertonset(i),2)'];    
    
    pert_time(i,:) = data.time{i}(data.ipertonset(i));
end

% trim down to keep size manageable
CrX_post = CrX_post(:,1:frameSize);
CrY_post = CrY_post(:,1:frameSize);
CrX_pre = CrX_pre(:,end-frameSize+1:end);
CrY_pre = CrY_pre(:,end-frameSize+1:end);

NrX_post = NrX_post(:,1:frameSize);
NrY_post = NrY_post(:,1:frameSize);
NrX_pre = NrX_pre(:,end-frameSize+1:end);
NrY_pre = NrY_pre(:,end-frameSize+1:end);

% potentially rid of trials that were a long way off a straight line to the
% target
if(eliminate_bad_trials)
    ibad = find(abs(data.x_at_jump)>.02);
    CrX_post(ibad,:) = NaN;
    CrY_post(ibad,:) = NaN;
    CrX_pre(ibad,:) = NaN;
    CrY_pre(ibad,:) = NaN;
end

% split data into subsets organized by jump type
psize = [-.03 0 .03];

% figure out which targets NEVER jumped
%pertAll = reshape(data.tFile(:,4),96,Nblocks); % align all similar trials across blocks
%nopert = sum(abs(pertAll),2)==0; % index of trials which never jumped
%nopert = repmat(nopert,Nblocks,1); % figure out overall trial number for each no perturbation trial
%inopert = find(nopert); % trials that were not perturbed
%data.pert(inopert)=NaN; % NaN out trials that are never perturbed

for p=1:length(psize)
    
    if(~isnan(psize(p)))
        ip = find(data.pert==psize(p)); % trials with perturbation size i
    else
        ip = find(isnan(data.pert)); % trials with no perturbation
    end
    
    % get cursor trajectory on these trials
    dd{p}.CrX_pre = CrX_pre(ip,:);
    dd{p}.CrX_post = CrX_post(ip,:);
    dd{p}.CrY_pre = CrY_pre(ip,:);
    dd{p}.CrY_post = CrY_post(ip,:);
    
    % null-space trajectory
    dd{p}.NrX_pre = NrX_pre(ip,:);
    dd{p}.NrX_post = NrX_post(ip,:);
    dd{p}.NrY_pre = NrY_pre(ip,:);
    dd{p}.NrY_post = NrY_post(ip,:);
    
    % event times for alignment
    dd{p}.pert_time = pert_time(ip);
    
    % additional data
    dd{p}.target_index = data.targIndex(ip,:);
    dd{p}.RT = data.RT(ip);
    dd{p}.pathlength = data.pathlength(ip);
    dd{p}.pathlength_null = data.pathlength_null(ip);
    dd{p}.movtime = data.movtime(ip);
    dd{p}.pkVel = data.pkVel(ip);
    dd{p}.vel_at_jump = data.vel_at_jump(ip,:);
    dd{p}.x_at_jump = data.vel_at_jump(ip,:);
end

data = dd;