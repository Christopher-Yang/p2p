function plot_traj(d)

% plot 10 reaches from baseline, early, and late learning

% set variables for plotting
subj = [2 6 1]; % select which participants to plot from each group
trials{1} = [1:10; 31:40; 331:340]; % trials to plot for 2-day group
trials{2} = [1:10; 31:40; 1251:1260]; % trials to plot for 5-day group
trials{3} = [1:10; 31:40; 2811:2820]; % trials to plot for 10-day group
groups = {'day2','day5','day10'};
titles = {'Baseline','Early','Late'};
groupNames = {'2-day','5-day','10-day'};
Ngroup = length(groupNames);
Nblock = length(titles);
col = [198 156 109
       0 146 69
       100 149 237
       251 176 59]./255; % colors for plotting

% generate plot
figure(1); clf
for k = 1:Ngroup % loop over groups
    s = subj(k);
    
    for j = 1:Nblock % loop over blocks
        subplot(3,3,(k-1)*3+j); hold on
        a = d.(groups{k}){s};
        
        for i = trials{k}(j,:) % loop over trials
            plot(a.targetAbs(i,1),a.targetAbs(i,2), '.', 'Color', [1 0.4 0.4], 'MarkerSize', 13); % plot targets
            plot(a.L{i}(:,1),a.L{i}(:,2),'Color',col(1,:)) % plot left hand
            plot(a.R{i}(:,1),a.R{i}(:,2),'Color',col(2,:)) % plot right hand
            plot(a.C{i}(:,1),a.C{i}(:,2),'k') % plot cursor
        end
        axis([0.1 1 -0.2 0.65])
        axis square
        if k == 1
            title(titles{j})
            if j == 1
                plot([0.4 0.52],[0.5 0.5],'k','LineWidth',4)
            end
        end
        if j == 1
            ylabel(groupNames{k})
        end
    end
end
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/traj','-dpdf','-painters')

%% plot reaches from flip block

% set variables for plotting
delay = 20; % time after movement initiation when velocity is computed
trialIdx = [2922 467]; % select trials to plot
axisLims = [0.43 0.73 0.1 0.4;
            0.4 0.7 0.15 0.45]; % set axes for plots

% generate plot
f = figure(2); clf
set(f,'Position',[200 200 350 200]);
for i = 1:2 % loop to generate plot for toward or away trial
    
    if i == 1
        a = d.day10{3};
    else
        a = d.day2{1};
    end
    
    % set variables
    trial = trialIdx(i); % trial to be plotted
    targ = a.targetAbs(trial-1:trial,:); % target trajectory
    curs = a.C{trial}; % cursor trajectory
    
    % calculate direction of instantaneous velocity vector
    init = a.init_x(trial)+delay; % index of cursor position 150 ms after movement initiation
    vel = diff(curs); % compute velocity (not divided by time because only the direction is needed)
    vel = vel(init,:);
    angle = atan2(vel(2),vel(1)); % find direction of velocity vector
    vector = 0.03*[cos(angle) sin(angle)]; % scale vector
    
    subplot(1,2,i); hold on
    plot([targ(1,1) targ(1,1)], [targ(1,2)-0.15 targ(1,2)+0.15], 'k--', 'LineWidth', 1) % mirror axis
    plot([curs(init,1) curs(init,1)+vector(1)], [curs(init,2) curs(init,2)+vector(2)], 'Color', col(3,:), 'LineWidth',4) % instantaneous velocity vector
    plot(targ(1,1), targ(1,2), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 35) % starting target
    plot(targ(2,1), targ(2,2), '.', 'Color', [1 0.4 0.4], 'MarkerSize', 35) % ending target
    plot(targ(1,1) - (targ(2,1) - targ(1,1)), targ(2,2), 'o', 'Color', [1 0.4 0.4], 'MarkerSize', 10) % mirrored ending target
    plot(curs(:,1), curs(:,2),'k','LineWidth',1) % cursor trajectory
    axis(axisLims(i,:))
    axis square
    xticks([])
    yticks([])

    if i == 2
        plot([0.46 0.58],[0.2 0.2],'k','LineWidth',4)
    end
end
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/traj_habit','-dpdf','-painters')

%%
trial = 2896;

a = d.day10{2};
curs = a.C{trial}; % cursor trajectory

vel = diff(curs); % compute velocity (not divided by time because only the direction is needed)

init = a.init(trial) + delay;
initVel = vel(init,:);
angle = atan2(initVel(2),initVel(1)); % find direction of velocity vector
vector = 0.03*[cos(angle) sin(angle)]; % scale vector

init_x = a.init_x(trial) + delay;
initVel_x = vel(init_x,:);
angle = atan2(initVel_x(2),initVel_x(1)); % find direction of velocity vector
vector_x = 0.03*[cos(angle) sin(angle)]; % scale vector

targ = a.targetAbs(trial-1:trial,:); % target trajectory

f = figure(3); clf; hold on
set(f,'Position',[200 200 100 100]);
plot([targ(1,1) targ(1,1)], [targ(1,2)-0.15 targ(1,2)+0.15], 'k--', 'LineWidth', 1) % mirror axis
plot([curs(init,1) curs(init,1)+vector(1)], [curs(init,2) curs(init,2)+vector(2)], 'Color', col(3,:), 'LineWidth',4) % instantaneous velocity vector
plot([curs(init_x,1) curs(init_x,1)+vector_x(1)], [curs(init_x,2) curs(init_x,2)+vector_x(2)], 'Color', col(4,:), 'LineWidth',4) % instantaneous horizontal velocity
plot(targ(1,1), targ(1,2), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 35) % starting target
plot(targ(2,1), targ(2,2), '.', 'Color', [1 0.4 0.4], 'MarkerSize', 35) % ending target    
plot(a.C{trial}(:,1),a.C{trial}(:,2),'k') % plot cursor
axis([0.55 0.7 0.2 0.35])
axis square
xticks([])
yticks([])

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/traj_habit2','-dpdf','-painters')

end
