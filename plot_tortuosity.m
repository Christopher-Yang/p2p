
clear tortuosity movtime RT
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
    tortuosity_incorrect{i} = [];
    movtime_incorrect{i} = [];
    RT_incorrect{i} = [];
    tortuosity_correct{i} = [];
    movtime_correct{i} = [];
    RT_correct{i} = [];
    
    for j = 1:Nsubj(i)
        a = d.(groups{i}){j};
        tortuosity{i}(:,j) = a.pathlength./0.12;
        movtime{i}(:,j) = a.movtime./1000;
        RT{i}(:,j) = a.RT;
        
        % need to use == 1/0 because some values are NaN
        incorrectIdx = a.incorrectReach_x(end-99:end) == 1;
        correctIdx = a.incorrectReach_x(end-99:end) == 0;

        last_tort = tortuosity{i}(habitBlocks(i,:),j);
        last_movtime = movtime{i}(habitBlocks(i,:),j);
        last_RT = RT{i}(habitBlocks(i,:),j);
        
        tortuosity_incorrect{i} = [tortuosity_incorrect{i}; last_tort(incorrectIdx)];
        movtime_incorrect{i} = [movtime_incorrect{i}; last_movtime(incorrectIdx)];
        RT_incorrect{i} = [RT_incorrect{i}; last_RT(incorrectIdx)];
        
        tortuosity_correct{i} = [tortuosity_correct{i}; last_tort(correctIdx)];
        movtime_correct{i} = [movtime_correct{i}; last_movtime(correctIdx)];
        RT_correct{i} = [RT_correct{i}; last_RT(correctIdx)];
    end
    
    Ntrials = size(RT{i},1);
    for j = 1:Ntrials/bin
        tortuosityBin{i}(j,:) = mean(tortuosity{i}((j-1)*5+1:(j-1)*5+5,:),1);
        movtimeBin{i}(j,:) = mean(movtime{i}((j-1)*5+1:(j-1)*5+5,:),1);
        RTBin{i}(j,:) = mean(RT{i}((j-1)*5+1:(j-1)*5+5,:),1);
    end
    
    tortuosity_mu{i} = mean(tortuosityBin{i},2);
    tortuosity_se{i} = std(tortuosityBin{i},[],2)/sqrt(Nsubj(i));
    movtime_mu{i} = mean(movtimeBin{i},2);
    movtime_se{i} = std(movtimeBin{i},[],2)/sqrt(Nsubj(i));
    RT_mu{i} = mean(RTBin{i},2);
    RT_se{i} = std(RTBin{i},[],2)/sqrt(Nsubj(i));
end

% col = [0 0 0
%        255 160 122
%        148 0 211]./255;
col = [180 180 0
       0 191 255
       255 99 71]./255;
gblock = [3 2 1];
trialIdx = [5 14 29];
lw = 0.25;

figure(1); clf; hold on
for i = 1:3
    avg = mean(tortuosity_mu{i}(1:6));
    plot([trials{trialIdx(i)-1}(1) trials{trialIdx(i)}(end)],[avg avg],'Color',col(i,:),'LineWidth',4)
end

for i = 1:5
    for j = 1:3
        s = shadedErrorBar(trials{i},tortuosity_mu{gblock(j)}(trials{i}), tortuosity_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},tortuosity_mu{gblock(j)}(trials{i}), tortuosity_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

for i = 15:29
    s = shadedErrorBar(trials{i},tortuosity_mu{3}(trials{i}), tortuosity_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end

xticks(dayStart)
xticklabels(dayStartLabels)
yticks(1:5)
axis([1 566 1 5])
xlabel('Trial Number')
ylabel('Tortuosity')
legend({'2-day','5-day','10-day'})
set(gca,'TickDir','out')

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/tortuosity','-dpdf','-painters')

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
    subplot(3,3,i); hold on
    histogram(tortuosity_correct{i},'Normalization','pdf')
    histogram(tortuosity_incorrect{i},'Normalization','pdf')
    
    subplot(3,3,i+3); hold on
    histogram(movtime_correct{i},'Normalization','pdf')
    histogram(movtime_incorrect{i},'Normalization','pdf')
    
    subplot(3,3,i+6); hold on
    histogram(RT_correct{i},'Normalization','pdf')
    histogram(RT_incorrect{i},'Normalization','pdf')
end