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
            num.(groups{i})(k,j) = nansum(a);
            den.(groups{i})(k,j) = 100-sum(isnan(a));
            habit.(groups{i})(k,j) = 100*num.(groups{i})(k,j)/den.(groups{i})(k,j);
        end
    end
end

x = [reshape(habit.day2(:,3:end), [numel(habit.day2(:,3:end)) 1]); ...
    reshape(habit.day5(:,3:end), [numel(habit.day5(:,3:end)) 1]); ...
    reshape(habit.day10(:,3:end), [numel(habit.day10(:,3:end)) 1])];

groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
blockNames([1:13 27:40 55:59],1) = "late";
blockNames([14:26 41:54 60:64],1) = "flip";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, blockNames, subject, x, 'VariableNames', {'group','block','subject','away'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/away.csv')

%% Figure 4G

col = [180 180 0
       0 191 255
       255 99 71]./255;

rng(34);
f = figure(1); clf
set(f,'Position',[200 200 140 150]); hold on

for i = 1:Ngroup
    n = Nsubj(i);
    plot(repmat((i-1) + [1 5],[n 1]) + 0.5*(rand(n,2) - 0.5), habit.(groups{i})(:,3:4), '.', 'MarkerSize', 12, 'Color', col(i,:))
    plot((i-1) + [1 5], mean(habit.(groups{i})(:,3:4),1), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'LineWidth', 1)
    if i == 1
        xticks([2 6])
        xticklabels({'Late','Flip'})
        xlabel('Block')
        set(gca,'TickDir','out')
        axis([0 8 0 60])
        ylabel('Aimed away (%)')
    end
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/away','-dpdf','-painters')

end