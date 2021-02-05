
col = [198 156 109
       0 146 69
       191 126 255]./255;
subj = [2 6 5];
trials{1} = [1:10; 31:40; 331:340];
trials{2} = [1:10; 31:40; 1251:1260]; 
trials{3} = [1:10; 31:40; 2811:2820];
groups = {'day2','day5','day10'};
titles = {'Baseline','Early','Late'};
groupNames = {'2-day','5-day','10-day'};

figure(1); clf
for k = 1:3
    s = subj(k);
    for j = 1:3
        subplot(3,3,(k-1)*3+j); hold on
        a = d.(groups{k}){s};
        for i = trials{k}(j,:)
            scatter(a.targetAbs(i,1),a.targetAbs(i,2),12,'filled','MarkerFaceColor',[179 179 179]./255);
            plot(a.L{i}(:,1),a.L{i}(:,2),'Color',col(1,:))
            plot(a.R{i}(:,1),a.R{i}(:,2),'Color',col(2,:))
            plot(a.C{i}(:,1),a.C{i}(:,2),'Color',col(3,:))
        end
        axis([0.1 0.95 -0.2 0.65])
        axis square
%         axis equal
        if k == 1
            title(titles{j})
        end
        if j == 1
            ylabel(groupNames{k})
        end
    end
end
% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/traj','-dpdf','-painters')

%%
group = 'day10';
subj = 3;
trial = 2920;

a = d.(group){subj};

init = a.init_x(trial)+20;
vel = diff(a.C{trial});
vel = vel(init,:);
angle = atan2(vel(2),vel(1));
vector = 0.03*[cos(angle) sin(angle)];

figure(2); clf
subplot(1,2,1); hold on
plot(a.C{trial}(:,1), a.C{trial}(:,2),'Color',col(3,:),'LineWidth',1)
scatter(a.targetAbs(trial,1), a.targetAbs(trial,2), 200, 'k', 'filled', 'MarkerFaceAlpha',0.4)
scatter(a.targetAbs(trial-1,1), a.targetAbs(trial-1,2), 200, 'r', 'filled', 'MarkerFaceAlpha',0.4)
plot([a.C{trial}(init,1) a.C{trial}(init,1)+vector(1)], [a.C{trial}(init,2) a.C{trial}(init,2)+vector(2)],'r')
axis([0.45 0.7 0.2 0.45])
axis square
title('Incorrect reach')

trial = 2922;
init = a.init_x(trial)+20;
vel = diff(a.C{trial});
vel = vel(init,:);
angle = atan2(vel(2),vel(1));
vector = 0.03*[cos(angle) sin(angle)];

subplot(1,2,2); hold on
plot(a.C{trial}(:,1), a.C{trial}(:,2),'Color',col(3,:),'LineWidth',1)
scatter(a.targetAbs(trial,1), a.targetAbs(trial,2), 200, 'k', 'filled', 'MarkerFaceAlpha',0.4)
scatter(a.targetAbs(trial-1,1), a.targetAbs(trial-1,2), 200, 'r', 'filled', 'MarkerFaceAlpha',0.4)
plot([a.C{trial}(init,1) a.C{trial}(init,1)+vector(1)], [a.C{trial}(init,2) a.C{trial}(init,2)+vector(2)],'r')
axis([0.40 0.65 0.10 0.35])
axis square
title('Correct reach')

print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/traj_habit','-dpdf','-painters')