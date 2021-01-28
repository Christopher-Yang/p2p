%% align velocity to target presentation

clear vel
blocks = {'Baseline','Early','Late'};
groups = {'day2','day5','day10'};
graph_names = {'2-day','5-day','10-day'};
Nsubj = [13 14 5];
Ntrials = [330 430 530];
Ngroup = 3;
trials{1} = 1:30;
trials{2} = 31:130;
trials{3} = 131:230;
trials{4} = 231:330;
trials{5} = 331:430;
trials{6} = 431:530;

for k = 1:Ngroup
    for j = 1:Nsubj(k)
        a = d.(groups{k}){j};
        Nmax(j) = max(cellfun('size',a.velFilt,1));
    end
    
    for j = 1:Nsubj(k)
        a = d.(groups{k}){j};
        for i = 1:Ntrials(k)
            v = a.velFilt{i}(a.itargonset(i):end);
            vel{k}(1:length(v),i,j) = v;
        end
    end
end

%% plot results
delt = 1/130;
gblock = [1 2 3 
          1 2 4 
          1 2 5];
col = [0 0 0
       0 191 255
       255 99 71]./255;

figure(1); clf
for i = 1:Ngroup
    xAxis = 0:delt:size(vel{i},1)*delt - delt;
    
    subplot(1,3,i); hold on
    for j = 1:3
        v = nanmean(vel{i}(:,trials{gblock(i,j)},:),2);
        s = shadedErrorBar(xAxis,nanmean(v,3),std(v,[],3)); % baseline
        editErrorBar(s,col(j,:),1);
    end
    axis([0 2 0 0.4])
    title(graph_names{i})
    xlabel('Time (s)')
    if i == 1
        ylabel('Velocity')
    elseif i == 3
        legend(blocks)
    end
end

%%
gblock = [3 4 5];

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
    v = reshape(vel{i}(:,trials{gblock(i)},:),[length(xAxis) 100*Nsubj(i)]);
    plot(xAxis,v,'Color',[0 0 0 0.02])
    plot(xAxis,nanmean(nanmean(vel{i}(:,trials{gblock(i)},:),2),3),'r','LineWidth',2)
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
