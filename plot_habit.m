function plot_habit(d)

% organize data into easily plottable variables
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
        end
    end
end

x = [reshape(habit.x.day2(:,3:end), [numel(habit.x.day2(:,3:end)) 1]); ...
    reshape(habit.x.day5(:,3:end), [numel(habit.x.day5(:,3:end)) 1]); ...
    reshape(habit.x.day10(:,3:end), [numel(habit.x.day10(:,3:end)) 1])];

groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
blockNames([1:13 27:40 55:59],1) = "late";
blockNames([14:26 41:54 60:64],1) = "flip";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, blockNames, subject, x, 'VariableNames', {'group','block','subject','away'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/away.csv')

%% plot movement kinematics and incorrect-movement trials

col = [180 180 0
       0 191 255
       255 99 71]./255;

rng(34);
ax = {'x','y'};
f = figure(1); clf
set(f,'Position',[200 200 270 150]);

for j = 1:2
    for i = 1:Ngroup
        subplot(1,2,j); hold on
        n = Nsubj(i);
        plot(repmat((i-1) + [1 5],[n 1]) + 0.5*(rand(n,2) - 0.5), habit.(ax{j}).(groups{i})(:,3:4), '.', 'MarkerSize', 12, 'Color', col(i,:))
        plot((i-1) + [1 5], mean(habit.(ax{j}).(groups{i})(:,3:4),1), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'LineWidth', 1)
        if i == 1
            xticks([2 6])
            xticklabels({'Late','Flip'})
            xlabel('Block')
            set(gca,'TickDir','out')
            axis([0 8 0 60])
            if j == 1
                ylabel('Aimed away (%)')
                title('Left hand')
            else
                title('Right hand')
            end
        end
    end
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/away','-dpdf','-painters')
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

figure(2); clf
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
            xticks(1:4)
            xticklabels({'Baseline','Early','Late','Post-flip'})
            if i == 1
                ylabel('10 day')
            end
        end
    end
end

figure(3); clf; hold on
for j = 1:Ngroup
    plot((j-1)*5+1:(j-1)*5+4, permute(habit2.(groups{j})(:,4,:),[3 1 2]),'k')
    plot((j-1)*5+1:(j-1)*5+4, squeeze(mean(habit2.(groups{j})(:,4,:),1)),'.r','MarkerSize',20)
end
xticks([1 4 6 9 11 14])
xticklabels({'Closest','Farthest','Closest','Farthest','Closest','Farthest'})
ylabel('Propotion of incorrect reaches in flipped block')
title('left: 2-day; center: 5-day; right: 10-day')

%% plot the lag between left and right hand movement (doesn't work as is)
% rng(1)
% figure(7); clf
% for i = 1:Ngroup
%     subplot(1,3,i); hold on
%     for j = 1:3
%         plot(j+0.3*(rand(Ntrial,Nsubj(i))-0.5),lag.(groups{i})(:,:,j),'k.','MarkerSize',5)
%         plot(j,nanmean(nanmean(lag.(groups{i})(:,:,j),1),2),'.r','MarkerSize',40)
%     end
%     if i == 1
%         title('2 days of training')
%         ylabel('LH start - RH start (s)')
%     elseif i == 2
%         title('5 days of training')
%     elseif i == 3
%         title('10 days of training')
%     end
%     xticks(1:3)
%     xticklabels({'Early learning','Late learning','Post-flip'})
%     ylim([-6 6])
% end

%% plot initial reach direction in x and y axes
% group = 'day10';
% subj = 1;
% trial = 2863;
% 
% a = d.(group){subj};
% t = (a.time{trial}-a.time{trial}(1))/1000;
% init = a.init_x(trial)+20;
% x = [t(init-20) t(init+20)];
% 
% y = a.C{trial}(init,1)-a.C{trial}(1,1);
% m = a.initVel_x(trial);
% b = y - m*t(init);
% 
% figure(8); clf; hold on
% plot(t(end),a.targetAbs(trial,1)-a.C{trial}(1,1),'k.','MarkerSize',70,'HandleVisibility','off')
% plot(t,a.C{trial}(:,1)-a.C{trial}(1,1),'k','LineWidth',1)
% plot(t(init),a.C{trial}(init,1)-a.C{trial}(1,1),'.k','MarkerSize',30,'HandleVisibility','off')
% % the commented out line plots the lag between left and right hand movement
% % on a given trial
% % plot([t(init) t(init)-a.lag(trial)/1000],[a.C{trial}(init,1)-a.C{trial}(1,1) a.C{trial}(init,1)-a.C{trial}(1,1)],'r','LineWidth',2)
% plot(x,[m*x(1)+b m*x(2)+b],'k','LineWidth',3,'HandleVisibility','off')
% 
% init = a.init_y(trial)+20;
% x = [t(init-20) t(init+20)];
% 
% y = a.C{trial}(init,2)-a.C{trial}(1,2);
% m = a.initVel_y(trial);
% b = y - m*t(init);
% 
% plot(t(end),a.targetAbs(trial,2)-a.C{trial}(1,2),'b.','MarkerSize',70,'HandleVisibility','off')
% plot(t,a.C{trial}(:,2)-a.C{trial}(1,2),'b','LineWidth',1)
% plot(t(init),a.C{trial}(init,2)-a.C{trial}(1,2),'.b','MarkerSize',30,'HandleVisibility','off')
% plot(x,[m*x(1)+b m*x(2)+b],'b','LineWidth',3,'HandleVisibility','off')
% xlabel('Time (s)')
% ylabel('Position relative to start (m)')
% ylim([-.15 .15])
% legend({'Horizontal position','Vertical position'})
end