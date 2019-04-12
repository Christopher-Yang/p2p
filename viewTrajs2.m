function viewTrajs2(data,varargin)
% inspect trajectories and data from RT experiment
%
% inputs: data - a data structure from an RT experiment (created by
%                loadTRsubjDATa
%         [minRT maxRT] - range of RTs to allow (default [0 0.5]);
%
if(length(varargin)==0)
    igood = find(data.RT>0 & data.RT < .5);
else
    igood = find(data.RT>varargin{1}(1) & data.RT<varargin{1}(2));
end
[xxrt i_sort_rt] = sort(data.RT(igood));
rng = igood(i_sort_rt);
%igood = rng;%find(data.RT(rng)>0 & data.RT(rng) < .5);

figure(1); clf; hold on
subplot(3,2,1); hold on
plot(data.RT,data.reachDir,'.','color',.9*ones(1,3),'markersize',20)
plot(data.RT(rng),data.reachDir(rng),'.','color',.7*ones(1,3),'markersize',20)
xlabel('Reaction Time')
ylabel('Reach Error')
axis([0 .5 -180 180])

subplot(3,2,3); hold on
plot(data.RT,data.reachDir_absolute,'.','color',.9*ones(1,3),'markersize',20)
plot(data.RT(rng),data.reachDir_absolute(rng),'.','color',.7*ones(1,3),'markersize',20)
xlabel('Reaction Time')
ylabel('Absolute Reach Direction')
axis([0 .5 -180 180])

for i=1:length(rng)
    ii = rng(i);
    
    % plot reach direction error
    subplot(3,2,1); hold on
    if(i>1)
        plot(data.RT(rng(i-1)),data.reachDir(rng(i-1)),'k.','markersize',20)
    end
    plot(data.RT(ii),data.reachDir(ii),'r.','markersize',20)
    
    % plot absolute reach direction
    subplot(3,2,3); hold on
    if(i>1)
        plot(data.RT(rng(i-1)),data.reachDir_absolute(rng(i-1)),'k.','markersize',20)
    end
    plot(data.RT(ii),data.reachDir_absolute(ii),'r.','markersize',20)
    
    % plot velocity profile
    subplot(3,2,5); cla; hold on
    plot([1:length(data.tanVel{ii})]/130,data.tanVel{ii},'linewidth',2)
    plot(data.iInit(ii)/130,data.tanVel{ii}(data.iInit(ii)),'b.','linewidth',2,'markersize',20)
    xlabel('Time')
    ylabel('Velocity')
    
    % plot full trajectory
    subplot(3,2,[2 4 6]); cla; hold on
    %plot(data.handPos{ii}(1,data.iInit(ii):data.iEnd(ii)),data.handPos{ii}(2,data.iInit(ii):data.iEnd(ii)),'r')
    plot(data.handPos_rotated{ii}(1,data.iInit(ii):data.iEnd(ii)),data.handPos_rotated{ii}(2,data.iInit(ii):data.iEnd(ii)),'g','linewidth',2)
    %plot(data.targPos(ii,1),data.targPos(ii,2),'o','markersize',20)
    plot(0,.08,'o','markersize',20)
    %plot(data.handPos{ii}(1,data.iDir(ii)),data.handPos{ii}(2,data.iDir(ii)),'r.','markersize',16)
    plot(data.handPos_rotated{ii}(1,data.iDir(ii)),data.handPos_rotated{ii}(2,data.iDir(ii)),'g.','markersize',20)
    dirV = .015*[-sin(data.reachDir(ii)*pi/180) cos(data.reachDir(ii)*pi/180)];
    plot(data.handPos_rotated{ii}(1,data.iDir(ii))*[1 1] + [0 dirV(1)],data.handPos_rotated{ii}(2,data.iDir(ii))*[1 1] + [0 dirV(2)],'k','linewidth',2.5)
    axis equal
    axis([-.08 .08 -.1 .1])
    
    text(-.06,.02,['reachDir = ',num2str(data.reachDir(ii))])
    text(-.06,.015,['reachDirAbsolute = ',num2str(data.reachDir_absolute(ii))])
    text(-.06,.01,['targAng = ',num2str(data.targAng(ii))])
    text(-.06,.025,['Trial #: ',num2str(ii)]);
    
    %keyboard
    pause
end
