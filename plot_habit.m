%% organize data into easily plottable variables
groups = {'day2','day5','day10'};
Nsubj = [length(d.day2) length(d.day5) length(d.day10)];
Ngroup = length(groups);

trials{1} = 1:30;
for i = 1:29
    trials{i+1} = (i-1)*100 + 31:(i-1)*100 + 130;
end

gblocks = [1 2 5 6
          1 2 14 15
          1 2 29 30];
      
Nblock = size(gblocks,2);

for i = 1:Ngroup
    for j = 1:Nblock
        for k = 1:Nsubj(i)
            trialIdx = trials{gblocks(i,j)};
            
            a = d.(groups{i}){k}.incorrectReach_x(trialIdx);
            num.x.(groups{i})(k,j) = nansum(a);
            den.x.(groups{i})(k,j) = 100-sum(isnan(a));
            habit.x.(groups{i})(k,j) = 100*num.x.(groups{i})(k,j)/den.x.(groups{i})(k,j);
            
            a = d.(groups{i}){k}.incorrectReach_y(trialIdx);
            num.y.(groups{i})(k,j) = nansum(a);
            den.y.(groups{i})(k,j) = 100-sum(isnan(a));
            habit.y.(groups{i})(k,j) = 100*num.y.(groups{i})(k,j)/den.y.(groups{i})(k,j);
            
            a = d.(groups{i}){k};
            pathLength.(groups{i})(k,j) = mean(a.pathlength(trialIdx));
            movTime.(groups{i})(k,j) = mean(a.movtime(trialIdx))/1000;
            RT.(groups{i})(k,j) = mean(a.RT(trialIdx))/1000;
            
%             lag.(groups{i})(:,k,j) = d.(groups{i}){k}.lag((j-1)*100+1:(j-1)*100+100)/1000;
        end
    end
end

%% plot initial reach direction in x and y axes
group = 'day10';
subj = 1;
trial = 2831;

a = d.(group){subj};
t = (a.time{trial}-a.time{trial}(1))/1000;
init = a.init_x(trial)+20;
x = [t(init-20) t(init+20)];

y = a.C{trial}(init,1)-a.C{trial}(1,1);
m = a.initVel_x(trial);
b = y - m*t(init);

figure(1); clf; hold on
plot(t(end),a.targetAbs(trial,1)-a.C{trial}(1,1),'k.','MarkerSize',70,'HandleVisibility','off')
plot(t,a.C{trial}(:,1)-a.C{trial}(1,1),'k','LineWidth',1)
plot(t(init),a.C{trial}(init,1)-a.C{trial}(1,1),'.k','MarkerSize',30,'HandleVisibility','off')
% the commented out line plots the lag between left and right hand movement
% on a given trial
% plot([t(init) t(init)-a.lag(trial)/1000],[a.C{trial}(init,1)-a.C{trial}(1,1) a.C{trial}(init,1)-a.C{trial}(1,1)],'r','LineWidth',2)
plot(x,[m*x(1)+b m*x(2)+b],'k','LineWidth',3,'HandleVisibility','off')

init = a.init_y(trial)+20;
x = [t(init-20) t(init+20)];

y = a.C{trial}(init,2)-a.C{trial}(1,2);
m = a.initVel_y(trial);
b = y - m*t(init);

plot(t(end),a.targetAbs(trial,2)-a.C{trial}(1,2),'b.','MarkerSize',70,'HandleVisibility','off')
plot(t,a.C{trial}(:,2)-a.C{trial}(1,2),'b','LineWidth',1)
plot(t(init),a.C{trial}(init,2)-a.C{trial}(1,2),'.b','MarkerSize',30,'HandleVisibility','off')
plot(x,[m*x(1)+b m*x(2)+b],'b','LineWidth',3,'HandleVisibility','off')
xlabel('Time (s)')
ylabel('Position relative to start (m)')
ylim([-.15 .15])
legend({'Horizontal position','Vertical position'})

%% plot movement kinematics and incorrect-movement trials

% figure(2); clf
% for i = 1:Ngroup
%     subplot(3,3,i); hold on
%     plot(pathLength.(groups{i})','Color',[0 0 0 0.6])
%     plot(mean(pathLength.(groups{i}),1),'.r','MarkerSize',20)
%     xticks([])   
%     ylim([0 1])
%     if i == 1
%         ylabel('Path length( m)')
%         title('2-day')
%     elseif i == 2
%         title('5-day')
%     elseif i == 3
%         title('10-day')
%     end
%     
%     subplot(3,3,i+3); hold on
%     plot(movTime.(groups{i})','Color',[0 0 0 0.6])
%     plot(mean(movTime.(groups{i}),1),'.r','MarkerSize',20)
%     xticks([])
%     ylim([0 10])
%     if i == 1
%         ylabel('Movement time (s)')
%     end
%     
%     subplot(3,3,i+6); hold on
%     plot(RT.(groups{i})','Color',[0 0 0 0.6])
%     plot(mean(RT.(groups{i}),1),'.r','MarkerSize',20)
%     ylim([0 2.5])
%     xticks(1:4)
%     xticklabels({'Baseline','Early','Late','Flip'})
%     if i == 1
%         ylabel('Reaction time (s)')
%     end
% end

col = [180 180 0
       0 191 255
       255 99 71]./255;

% figure(3); clf
% for i = 1:Ngroup
%     subplot(2,3,i); hold on
%     plot(habit.x.(groups{i})','Color',[col(i,:) 0.5])
%     plot(mean(habit.x.(groups{i}),1),'.','Color',col(i,:),'MarkerSize',20)
%     xticks([])
%     ylim([0 60])
%     set(gca,'TickDir','out')
%     if i == 1
%         ylabel('Proportion of incorrect reaches: left hand')
%         title('2-day')
%     elseif i == 2
%         title('5-day')
%     elseif i == 3
%         title('10-day')
%     end
%     
%     subplot(2,3,i+3); hold on
%     plot(habit.y.(groups{i})','Color',[col(i,:) 0.7])
%     plot(mean(habit.y.(groups{i}),1),'.','Color',col(i,:),'MarkerSize',20)
%     xticks(1:4)
%     xticklabels({'Baseline','Early','Late','Flip'})
%     ylim([0 60])
%     set(gca,'TickDir','out')
%     if i == 1
%         ylabel('Proportion of incorrect reaches: right hand')
%     end
% end

rng(34);
ax = {'x','y'};
xAxis = [0.5 2.5 5 15
         1 3 8 15.5
         1.5 3.5 13 16];
figure(3); clf
for j = 1:2
    for i = 1:Ngroup
        subplot(1,2,j); hold on
        n = Nsubj(i);
%         plot(repmat([0 6 12 18],[n 1]) + (rand(n,4) - 0.5) + 1.5*i, habit.(ax{j}).(groups{i}), '.', 'MarkerSize', 20, 'Color', col(i,:))
%         plot([0 6 12 18] + 1.5*i, mean(habit.(ax{j}).(groups{i}),1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1)
        plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), habit.(ax{j}).(groups{i}), '.', 'MarkerSize', 20, 'Color', col(i,:))
        plot(xAxis(i,:), mean(habit.(ax{j}).(groups{i}),1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', col(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1)
        if i == 1
            xticks([1 3 5 8 13 15.5])
            xticklabels({'Baseline',1,2,5,10,'Flip'})
            xlabel('Day')
            set(gca,'TickDir','out')
            axis([0 16.5 0 60])
            if j == 1
                ylabel('Proportion of incorrect reaches')
                title('Left hand')
            else
                title('Right hand')
            end
        end
    end
end

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/habit','-dpdf','-painters')
%% plot proportion of incorrect reaches binned by distance from mirroring axis

for m = 1:4 % loop over target direction bins
    for i = 1:Ngroup
        for j = 1:Nblock
            trialIdx = trials{gblocks(i,j)};
            for k = 1:Nsubj(i)
                data = d.(groups{i}){k};
                idx = find(data.targBin(trialIdx) == m);
                overlap = data.incorrectReach_x(trialIdx(1)-1 + idx);
                habit2.(groups{i})(k,j,m) = nansum(overlap)/sum(~isnan(overlap));
            end
        end
    end
end

figure(4); clf
for i = 1:4
    for j = 1:Ngroup
        subplot(3,4,i+((j-1)*4)); hold on
        plot(habit2.(groups{j})(:,:,i)','Color',[0 0 0 0.6])
        plot(mean(habit2.(groups{j})(:,:,i),1),'.r','MarkerSize',20)
        ylim([0 0.7])
        xticks([])
        
        if j == 1
            if i == 1
                title('Closest')
                ylabel('2 day')
            elseif i == 2
                title('Close')
            elseif i == 3
                title('Far')
            else
                title('Farthest')
            end
        elseif j == 2 && i == 1
            ylabel('5 day')
        elseif j == 3
            xticks(1:3)
            xticklabels({'Baseline','Early','Late','Post-flip'})
            if i == 1
                ylabel('10 day')
            end
        end
    end
end

figure(5); clf; hold on
for j = 1:Ngroup
    plot((j-1)*5+1:(j-1)*5+4, permute(habit2.(groups{j})(:,4,:),[3 1 2]),'k')
    plot((j-1)*5+1:(j-1)*5+4, squeeze(mean(habit2.(groups{j})(:,4,:),1)),'.r','MarkerSize',20)
end
xticks([1 4 6 9 11 14])
xticklabels({'Closest','Farthest','Closest','Farthest','Closest','Farthest'})
ylabel('Propotion of incorrect reaches in flipped block')
title('left: 2-day; center: 5-day; right: 10-day')

%% plot initial velocity sorted by correct/incorrect reaches
for i = 1:Ngroup
    for k = 1:Nsubj(i)
        baseline = trials{gblocks(i,1)};
        early = trials{gblocks(i,2)};
        late = trials{gblocks(i,end-1)};
        flip = trials{gblocks(i,end)};
        
        trialIdx = {baseline, early, late, flip};
        for j = 1:length(trialIdx)
            z = trialIdx{j};
            incorrectIdx = d.(groups{i}){k}.incorrectReach_x(z) == 1;        
            correctIdx = d.(groups{i}){k}.incorrectReach_x(z) == 0;
            
            vel = d.(groups{i}){k}.initVel_filt(z);
            initVel.(groups{i})(j,k) = nanmean(vel);
            incorrectVel.(groups{i})(j,k) = nanmean(vel(incorrectIdx));
            correctVel.(groups{i})(j,k) = nanmean(vel(correctIdx));
        end
    end
end

col = lines;
col = col(1:7,:);

figure(6); clf
subplot(1,2,1); hold on
for i = 1:Ngroup
    for j = 1:4
        errorbar((j-1)*5+i,nanmean(initVel.(groups{i})(j,:)),nanstd(initVel.(groups{i})(j,:))/sqrt(Nsubj(i)),'.','LineWidth',2,'MarkerSize',20,'Color',col(i,:))
    end
end
xticks(2:5:17)
xticklabels({'Baseline','Early','Late','Flip'})
ylabel('Velocity (m/s)')
ylim([0 0.35])

subplot(1,2,2); hold on
for i = 1:Ngroup
    errorbar(i,nanmean(incorrectVel.(groups{i})(4,:)),nanstd(incorrectVel.(groups{i})(4,:))/sqrt(Nsubj(i)),'.','LineWidth',2,'MarkerSize',20,'Color',col(i,:))
    errorbar(i+4,nanmean(correctVel.(groups{i})(4,:)),nanstd(correctVel.(groups{i})(4,:))/sqrt(Nsubj(i)),'.','LineWidth',2,'MarkerSize',20,'Color',col(i,:),'HandleVisibility','off')
end
title('Data from just the flip block')
xticks([2 6])
xticklabels({'Incorrect','Correct'})
ylim([0 0.35])
legend({'2-day','5-day','10-day'})

%% plot the lag between left and right hand movement (doesn't work as is)
rng(1)
figure(7); clf
for i = 1:Ngroup
    subplot(1,3,i); hold on
    for j = 1:3
        plot(j+0.3*(rand(Ntrial,Nsubj(i))-0.5),lag.(groups{i})(:,:,j),'k.','MarkerSize',5)
        plot(j,nanmean(nanmean(lag.(groups{i})(:,:,j),1),2),'.r','MarkerSize',40)
    end
    if i == 1
        title('2 days of training')
        ylabel('LH start - RH start (s)')
    elseif i == 2
        title('5 days of training')
    elseif i == 3
        title('10 days of training')
    end
    xticks(1:3)
    xticklabels({'Early learning','Late learning','Post-flip'})
    ylim([-6 6])
end
