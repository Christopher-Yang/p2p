%% plot initial reach direction in x and y axes
group = 'day10';
subj = 2;
trial = 220;

a = d.(group){subj};
t = (a.time{trial}-a.time{trial}(1))/1000;
init = a.init_x(trial)+20;
x = [t(init-30) t(init+30)];

y = a.C{trial}(init,1)-a.C{trial}(1,1);
m = a.initDir_x(trial);
b = y - m*t(init);

figure(1); clf; hold on
plot(t(end),a.targetAbs(trial,1)-a.C{trial}(1,1),'k.','MarkerSize',70)
plot(t,a.C{trial}(:,1)-a.C{trial}(1,1),'k','LineWidth',1)
plot(t(init),a.C{trial}(init,1)-a.C{trial}(1,1),'.k','MarkerSize',30)
% the commented out line plots the lag between left and right hand movement
% on a given trial
% plot([t(init) t(init)-a.lag(trial)/1000],[a.C{trial}(init,1)-a.C{trial}(1,1) a.C{trial}(init,1)-a.C{trial}(1,1)],'r','LineWidth',2)
plot(x,[m*x(1)+b m*x(2)+b],'k','LineWidth',3)

init = a.init_y(trial)+20;
x = [t(init-30) t(init+30)];

y = a.C{trial}(init,2)-a.C{trial}(1,2);
m = a.initDir_y(trial);
b = y - m*t(init);

plot(t(end),a.targetAbs(trial,2)-a.C{trial}(1,2),'b.','MarkerSize',70)
plot(t,a.C{trial}(:,2)-a.C{trial}(1,2),'b','LineWidth',1)
plot(t(init),a.C{trial}(init,2)-a.C{trial}(1,2),'.b','MarkerSize',30)
plot(x,[m*x(1)+b m*x(2)+b],'b','LineWidth',3)
xlabel('Time (s)')
ylabel('Horizontal position (m)')
ylim([-.15 .15])
%% organize data into easily plottable variables
groups = {'day2','day5','day10'};
Nsubj = [12 8 4];
Nblock = 3;
Ngroup = length(groups);
Ntrial = 100;

for i = 1:Ngroup
    for j = 1:Nblock
        for k = 1:Nsubj(i)
            a = d.(groups{i}){k}.incorrectReach_x((j-1)*100+1:(j-1)*100+100);
            num.x.(groups{i})(k,j) = nansum(a);
            den.x.(groups{i})(k,j) = 100-sum(isnan(a));
            habit.x.(groups{i})(k,j) = 100*num.x.(groups{i})(k,j)/den.x.(groups{i})(k,j);
            
            a = d.(groups{i}){k}.incorrectReach_y((j-1)*100+1:(j-1)*100+100);
            num.y.(groups{i})(k,j) = nansum(a);
            den.y.(groups{i})(k,j) = 100-sum(isnan(a));
            habit.y.(groups{i})(k,j) = 100*num.y.(groups{i})(k,j)/den.y.(groups{i})(k,j);
            
            pathLength.(groups{i})(k,j) = mean(d.(groups{i}){k}.pathlength((j-1)*100+1:(j-1)*100+100));
            movTime.(groups{i})(k,j) = mean(d.(groups{i}){k}.movtime((j-1)*100+1:(j-1)*100+100))/1000;
            RT.(groups{i})(k,j) = mean(d.(groups{i}){k}.RT((j-1)*100+1:(j-1)*100+100))/1000;
            
            lag.(groups{i})(:,k,j) = d.(groups{i}){k}.lag((j-1)*100+1:(j-1)*100+100)/1000;
        end
    end
end

%% plot movement kinematics and incorrect-movement trials
figure(2); clf
for i = 1:Ngroup
    subplot(4,3,i); hold on
    plot(habit.x.(groups{i})','Color',[0 0 0 0.6])
    plot(mean(habit.x.(groups{i}),1),'.k','MarkerSize',20)
    xticks([])
    ylim([0 60])
    if i == 1
        ylabel('Incorrect reach (%)')
        title('2 days of training')
    elseif i == 2
        title('5 days of training')
    elseif i == 3
        title('10 days of training')
    end
    
    subplot(4,3,i+3); hold on
    plot(pathLength.(groups{i})','Color',[0 0 0 0.6])
    plot(mean(pathLength.(groups{i}),1),'.k','MarkerSize',20)
    xticks([])
    ylim([0 0.6])
    if i == 1
        ylabel('Path length (m)')
    end
    
    subplot(4,3,i+6); hold on
    plot(movTime.(groups{i})','Color',[0 0 0 0.6])
    plot(mean(movTime.(groups{i}),1),'.k','MarkerSize',20)
    xticks([])
    ylim([0 6])
    if i == 1
        ylabel('Movement time (s)')
    end
    
    subplot(4,3,i+9); hold on
    plot(RT.(groups{i})','Color',[0 0 0 0.6])
    plot(mean(RT.(groups{i}),1),'.k','MarkerSize',20)
    xticks(1:3)
    xticklabels({'Early learning','Late learning','Post-flip'})
    ylim([0 1.5])
    if i == 1
        ylabel('Reaction time (s)')
    end
end

% figure(3); clf
% for i = 1:Ngroup
%     subplot(4,3,i); hold on
%     plot(habit.x.(groups{i})(:,[1 3])','Color',[0 0 0 0.6])
%     plot(mean(habit.x.(groups{i})(:,[1 3]),1),'.k','MarkerSize',20)
%     ylim([0 60])
%     xticks([])
%     if i == 1
%         ylabel('Incorrect reach (%)')
%         title('2 days')
%     elseif i == 2
%         title('5 days')
%     elseif i == 3
%         title('10 days')
%     end
%     
%     subplot(4,3,i+3); hold on
%     plot(pathLength.(groups{i})(:,[1 3])','Color',[0 0 0 0.6])
%     plot(mean(pathLength.(groups{i})(:,[1 3]),1),'.k','MarkerSize',20)
%     xticks([])   
%     ylim([0 0.6])
%     if i == 1
%         ylabel('Path length (m)')
%     end
%     
%     subplot(4,3,i+6); hold on
%     plot(movTime.(groups{i})(:,[1 3])','Color',[0 0 0 0.6])
%     plot(mean(movTime.(groups{i})(:,[1 3]),1),'.k','MarkerSize',20)
%     ylim([0 6])
%     if i == 1
%         ylabel('Movement time (s)')
%     end
%     
%     subplot(4,3,i+9); hold on
%     plot(RT.(groups{i})(:,[1 3])','Color',[0 0 0 0.6])
%     plot(mean(RT.(groups{i})(:,[1 3]),1),'.k','MarkerSize',20)
%     ylim([0 1.5])
%     xticks(1:2)
%     xticklabels({'Early learning','Post-flip'})
%     if i == 1
%         ylabel('Reaction time (s)')
%     end
% end

figure(3); clf
for i = 1:Ngroup
    subplot(2,3,i); hold on
    plot(habit.x.(groups{i})','Color',[0 0 0 0.6])
    plot(mean(habit.x.(groups{i}),1),'.k','MarkerSize',20)
    xticks([])
    ylim([0 60])
    if i == 1
        ylabel('Incorrect reach: left hand (%)')
        title('2 days')
    elseif i == 2
        title('5 days')
    elseif i == 3
        title('10 days')
    end
    
    subplot(2,3,i+3); hold on
    plot(habit.y.(groups{i})','Color',[0 0 0 0.6])
    plot(mean(habit.y.(groups{i}),1),'.k','MarkerSize',20)
    xticks([])
    ylim([0 60])
    ylabel('Incorrect reach: right hand (%)')
end

%% plot the lag between left and right hand movement
rng(1)
figure(4); clf
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

