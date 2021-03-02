
clear pLength movtime RT
groups = {'day2','day5','day10'};
Ngroups = length(groups);
Nsubj = [13 14 5];
bin = 5;
dayStart = [1 47:60:527];
dayStartLabels = [1 231:300:2930];
habitBlocks = [431:530; 1331:1430; 2831:2930];

trials{1} = 1:6;
for i = 1:29
    trials{i+1} = (i-1)*20 + 7:(i-1)*20 + 26;
end

for i = 1:Ngroups
    pLength_incorrect{i} = [];
    movtime_incorrect{i} = [];
    RT_incorrect{i} = [];
    pLength_correct{i} = [];
    movtime_correct{i} = [];
    RT_correct{i} = [];
    
    for j = 1:Nsubj(i)
        a = d.(groups{i}){j};
        pLength{i}(:,j) = a.pathlength;
        movtime{i}(:,j) = a.movtime./1000;
        RT{i}(:,j) = a.RT;
        
        % need to use == 1/0 because some values are NaN
        incorrectIdx = a.incorrectReach_x(end-99:end) == 1;
        correctIdx = a.incorrectReach_x(end-99:end) == 0;

        last_tort = pLength{i}(habitBlocks(i,:),j);
        last_movtime = movtime{i}(habitBlocks(i,:),j);
        last_RT = RT{i}(habitBlocks(i,:),j);
        
        pLength_incorrect{i} = [pLength_incorrect{i}; last_tort(incorrectIdx)];
        movtime_incorrect{i} = [movtime_incorrect{i}; last_movtime(incorrectIdx)];
        RT_incorrect{i} = [RT_incorrect{i}; last_RT(incorrectIdx)];
        
        pLength_correct{i} = [pLength_correct{i}; last_tort(correctIdx)];
        movtime_correct{i} = [movtime_correct{i}; last_movtime(correctIdx)];
        RT_correct{i} = [RT_correct{i}; last_RT(correctIdx)];
    end
    
    Ntrials = size(RT{i},1);
    for j = 1:Ntrials/bin
        pLengthBin{i}(j,:) = mean(pLength{i}((j-1)*5+1:(j-1)*5+5,:),1);
        movtimeBin{i}(j,:) = mean(movtime{i}((j-1)*5+1:(j-1)*5+5,:),1);
        RTBin{i}(j,:) = mean(RT{i}((j-1)*5+1:(j-1)*5+5,:),1);
    end
    
    pLength_mu{i} = mean(pLengthBin{i},2);
    pLength_se{i} = std(pLengthBin{i},[],2)/sqrt(Nsubj(i));
    movtime_mu{i} = mean(movtimeBin{i},2);
    movtime_se{i} = std(movtimeBin{i},[],2)/sqrt(Nsubj(i));
    RT_mu{i} = mean(RTBin{i},2);
    RT_se{i} = std(RTBin{i},[],2)/sqrt(Nsubj(i));
end

col = [180 180 0
       0 191 255
       255 99 71]./255;
gblock = [3 2 1];
trialIdx = [5 14 29];
lw = 0.25;

figure(1); clf; hold on
for i = 1:3
    avg = mean(pLength_mu{i}(1:6));
    plot([trials{trialIdx(i)-1}(1) trials{trialIdx(i)}(end)],[avg avg],'Color',col(i,:),'LineWidth',4)
end

for i = 1:5
    for j = 1:3
        s = shadedErrorBar(trials{i},pLength_mu{gblock(j)}(trials{i}), pLength_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},pLength_mu{gblock(j)}(trials{i}), pLength_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 15:29
    s = shadedErrorBar(trials{i},pLength_mu{3}(trials{i}), pLength_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end

xticks(dayStart)
xticklabels(dayStartLabels)
yticks(0:0.25:1)
axis([1 566 0 1])
xlabel('Trial Number')
ylabel('pLength')
legend({'2-day','5-day','10-day'})
set(gca,'TickDir','out')

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/pLength','-dpdf','-painters')

figure(2); clf; hold on
for i = 1:3
    avg = mean(movtime_mu{i}(1:6));
    plot([trials{trialIdx(i)-1}(1) trials{trialIdx(i)}(end)],[avg avg],'Color',col(i,:),'LineWidth',4)
end

for i = 1:5
    for j = 1:3
        s = shadedErrorBar(trials{i},movtime_mu{gblock(j)}(trials{i}), movtime_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},movtime_mu{gblock(j)}(trials{i}), movtime_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 15:29
    s = shadedErrorBar(trials{i},movtime_mu{3}(trials{i}), movtime_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end

xticks(dayStart)
xticklabels(dayStartLabels)
yticks(0:2:8)
axis([1 566 0 8])
xlabel('Trial Number')
ylabel('Movement time (s)')
legend({'2-day','5-day','10-day'})
set(gca,'TickDir','out')

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/movtime','-dpdf','-painters')

figure(3); clf; hold on
for i = 1:3
    avg = mean(RT_mu{i}(1:6));
    plot([trials{trialIdx(i)-1}(1) trials{trialIdx(i)}(end)],[avg avg],'Color',col(i,:),'LineWidth',4)
end

for i = 1:5
    for j = 1:3
        s = shadedErrorBar(trials{i},RT_mu{gblock(j)}(trials{i}), RT_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},RT_mu{gblock(j)}(trials{i}), RT_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 15:29
    s = shadedErrorBar(trials{i},RT_mu{3}(trials{i}), RT_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end

xticks(dayStart)
xticklabels(dayStartLabels)
yticks(400:400:1800)
axis([1 566 400 1800])
xlabel('Trial Number')
ylabel('Reaction time (ms)')
legend({'2-day','5-day','10-day'})
set(gca,'TickDir','out')

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/RT','-dpdf','-painters')

%%
figure(4); clf
for i = 1:3
    subplot(3,3,(i-1)*3+1); hold on
    histogram(pLength_correct{i},0:0.25:8,'Normalization','pdf')
    histogram(pLength_incorrect{i},0:0.25:8,'Normalization','pdf')
    ylim([0 1.6])
    if i == 3
        xlabel('pLength')
    end
    
    subplot(3,3,(i-1)*3+2); hold on
    histogram(movtime_correct{i},0:0.25:10,'Normalization','pdf')
    histogram(movtime_incorrect{i},0:0.25:10,'Normalization','pdf')
    ylim([0 1])
    if i == 3
        xlabel('Movement time (s)')
    end
    
    subplot(3,3,(i-1)*3+3); hold on
    histogram(RT_correct{i},0:50:2000,'Normalization','pdf')
    histogram(RT_incorrect{i},0:50:2000,'Normalization','pdf')
    ylim([0 0.005])
    if i == 3
        xlabel('Reaction time (ms)')
    end
end

%% bin reaches by the target's location
for i = 1:Ngroup
    for j = 1:4 % loop over target direction bins
        trialIdx = habitBlocks(i,1);
        idx = find(d.(groups{i}){1}.targBin == j);
        idx = idx(idx >= trialIdx);
        pLength2{i}{j} = reshape(pLength{i}(idx,:), [length(idx)*Nsubj(i) 1]);
        movtime2{i}{j} = reshape(movtime{i}(idx,:), [length(idx)*Nsubj(i) 1]);
        RT2{i}{j} = reshape(RT{i}(idx,:), [length(idx)*Nsubj(i) 1]);
    end
    pLength2_mu{i} = cellfun(@mean,pLength2{i});
    movtime2_mu{i} = cellfun(@mean,movtime2{i});
    RT2_mu{i} = cellfun(@mean,RT2{i});

    pLength2_se{i} = cellfun(@std,pLength2{i})./sqrt(Nsubj(i));
    movtime2_se{i} = cellfun(@std,movtime2{i})./sqrt(Nsubj(i));
    RT2_se{i} = cellfun(@std,RT2{i})./sqrt(Nsubj(i));
end

figure(5); clf
subplot(2,3,4); hold on
for i = 1:Ngroup
    errorbar((i-1)*5 + (1:4),pLength2_mu{i},pLength2_se{i},'.','MarkerSize',20)
end
xticks([1 4 6 9 11 14])
xticklabels({'Close','Far','Close','Far','Close','Far'})
ylabel('pLength')
ylim([1 3])
title('Post')


subplot(2,3,5); hold on
for i = 1:Ngroup
    errorbar((i-1)*5 + (1:4),movtime2_mu{i},movtime2_se{i},'.','MarkerSize',20)
end
xticks([1 4 6 9 11 14])
xticklabels({'Close','Far','Close','Far','Close','Far'})
ylabel('Movement time (s)')
ylim([0.5 3.5])
title('Post')

subplot(2,3,6); hold on
for i = 1:Ngroup
    errorbar((i-1)*5 + (1:4),RT2_mu{i},RT2_se{i},'.','MarkerSize',20)
end
xticks([1 4 6 9 11 14])
xticklabels({'Close','Far','Close','Far','Close','Far'})
ylabel('Reaction Time (ms)')
ylim([400 900])
title('Post')

for i = 1:Ngroup
    for j = 1:4 % loop over target direction bins
        trialIdx = habitBlocks(i,1)-100;
        idx = find(d.(groups{i}){1}.targBin == j);
        idx = idx(idx >= trialIdx);
        pLength2{i}{j} = reshape(pLength{i}(idx,:), [length(idx)*Nsubj(i) 1]);
        movtime2{i}{j} = reshape(movtime{i}(idx,:), [length(idx)*Nsubj(i) 1]);
        RT2{i}{j} = reshape(RT{i}(idx,:), [length(idx)*Nsubj(i) 1]);
    end
    pLength2_mu{i} = cellfun(@mean,pLength2{i});
    movtime2_mu{i} = cellfun(@mean,movtime2{i});
    RT2_mu{i} = cellfun(@mean,RT2{i});

    pLength2_se{i} = cellfun(@std,pLength2{i})./sqrt(Nsubj(i));
    movtime2_se{i} = cellfun(@std,movtime2{i})./sqrt(Nsubj(i));
    RT2_se{i} = cellfun(@std,RT2{i})./sqrt(Nsubj(i));
end

subplot(2,3,1); hold on
for i = 1:Ngroup
    errorbar((i-1)*5 + (1:4),pLength2_mu{i},pLength2_se{i},'.','MarkerSize',20)
end
xticks([])
xticklabels({'Close','Far','Close','Far','Close','Far'})
ylabel('pLength')
ylim([1 3])
title('Late')

subplot(2,3,2); hold on
for i = 1:Ngroup
    errorbar((i-1)*5 + (1:4),movtime2_mu{i},movtime2_se{i},'.','MarkerSize',20)
end
xticks([])
ylabel('Movement time (s)')
ylim([0.5 3.5])
title('Late')

subplot(2,3,3); hold on
for i = 1:Ngroup
    errorbar((i-1)*5 + (1:4),RT2_mu{i},RT2_se{i},'.','MarkerSize',20)
end
xticks([])
ylabel('Reaction Time (ms)')
ylim([400 900])
title('Late')
legend(groups,'Location','northwest')
