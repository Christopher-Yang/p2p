%% align velocity to target presentation

clear vel 
groups = {'day2','day5','day10'};
graph_names = {'2-day','5-day','10-day'};
Nsubj = [13 14 5];
Ngroup = 3;

% indices for dividing up the trials into blocks
trials{1} = 1:30;
for i = 1:29
    trials{i+1} = (i-1)*100 + 31:(i-1)*100 + 130;
end

% blocks for baseline, early, and late for each group
gblocks{1} = [1 2 3 5];
gblocks{2} = [1 2 3 5 9 12 14];
gblocks{3} = [1 2 3 5 9 12 15 18 21 24 27 29];

for k = 1:Ngroup
    Nmax = [];
    for j = 1:Nsubj(k)
        a = d.(groups{k}){j};
        Nmax(j) = max(cellfun('size',a.velFilt,1));
    end
    
    vel{k} = NaN(max(Nmax),length(gblocks{k})*100 - 70,Nsubj(k));
    for j = 1:Nsubj(k)
        a = d.(groups{k}){j};
        trialIdx = [trials{gblocks{k}}];
        for i = 1:length(trialIdx)
            v = a.velFilt{trialIdx(i)}(a.itargonset(trialIdx(i)):end);
            vel{k}(1:length(v),i,j) = v;
        end
    end
end

%% plot results

delt = 1/130; % sampling period

figure(1); clf
for i = 1:Ngroup % loop over groups
    
    xAxis = 0:delt:size(vel{i},1)*delt - delt; % times at which data were collected
    Nblock = length(gblocks{i}); % number of blocks to plot for each group
    
    % set colors to plot lines
    col1 = [255 0 0]./255;
    col2 = [255 215 0]./255;
    col = [linspace(col1(1),col2(1),Nblock-1)', linspace(col1(2),col2(2),Nblock-1)', linspace(col1(3),col2(3),Nblock-1)'];
    
    subplot(1,3,i); hold on
    
    % plot all data except baseline and last day of learning
    for j = 2:Nblock-1
        v = nanmean(vel{i}(:,trials{j},:),2)*100;
        s = shadedErrorBar(xAxis,nanmean(v,3),nanstd(v,[],3)/sqrt(Nsubj(i)));
        editErrorBar(s,col(j-1,:),2);
    end
    
    % plot baseline
    v = nanmean(vel{i}(:,1:30,:),2)*100;
    s = shadedErrorBar(xAxis,nanmean(v,3),nanstd(v,[],3)/sqrt(Nsubj(i)));
    editErrorBar(s,[0 0 0],2);
    
    % plot last day of learning
    v = nanmean(vel{i}(:,end-99:end,:),2)*100;
    s = shadedErrorBar(xAxis,nanmean(v,3),nanstd(v,[],3)/sqrt(Nsubj(i)));
    editErrorBar(s,col(end,:),2);
    
    axis([0 2 0 40])
    yticks(0:10:40)
    title(graph_names{i})
    xlabel('Time from target onset (s)')
    if i == 1
        ylabel('Tangential velocity (cm/s)')
    end
end

%%
figure(2); clf
for i = 1:3
    xAxis = 0:delt:size(vel{i},1)*delt - delt;
    
    subplot(2,3,i); hold on % plot baseline
    v = reshape(vel{i}(:,1:30,:),[length(xAxis) 30*Nsubj(i)]);
    plot(xAxis,v,'Color',[0 0 0 0.02])
    plot(xAxis,nanmean(nanmean(vel{i}(:,1:30,:),2),3),'r','LineWidth',2)
    axis([0 2 0 0.6])
    title(graph_names{i})
    if i == 1
        ylabel('Velocity')
    end
    
    subplot(2,3,i+3); hold on % plot late learning
    v = reshape(vel{i}(:,131:230,:),[length(xAxis) 100*Nsubj(i)]);
    plot(xAxis,v,'Color',[0 0 0 0.02])
    plot(xAxis,nanmean(nanmean(vel{i}(:,131:230,:),2),3),'r','LineWidth',2)
    axis([0 2 0 0.6])
    xlabel('Time (s)')
    if i == 1
        ylabel('Velocity')
    end
end

%% align velocity to movement onset
clear vel
groups = {'day2','day5','day10'};
Ngroup = 3;
Nsubj = [13 14 5];
Ntrials = [330 430 530];

for k = 1:Ngroup
    for j = 1:Nsubj(k)
        a = d.(groups{k}){j};
        
        before(j) = min(a.imoveonset);
        afterLength = NaN(Ntrials(k),1);
        
        for i = 1:Ntrials(k)
            afterLength(i) = length(a.velFilt{i}) - a.imoveonset(i);
        end
        
        after(j) = min(afterLength);
    end
    
    before_min(k) = min(before);
    after_min(k) = min(after);
    
    for j = 1:Nsubj(k)
        a = d.(groups{k}){j};
        
        for i = 1:Ntrials(k)
            idx = a.imoveonset(i) - before_min(k) + 1:a.imoveonset(i) + after_min(k);
            vel{k}(:,i,j) = a.velFilt{i}(idx);
        end
    end
end

%% plot results
gblock = [3 4 5];

figure(2); clf
for i = 1:Ngroup
    subplot(1,3,i); hold on
    plot(nanmean(nanmean(vel{i}(:,1:30,:),2),3),'k')
    plot(nanmean(nanmean(vel{i}(:,31:130,:),2),3))
    plot(nanmean(nanmean(vel{i}(:,trials{gblock(i)},:),2),3))
    ylim([0 0.5])
end
