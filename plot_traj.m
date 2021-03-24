function plot_traj(d)

% plot 10 reaches from baseline, early, and late learning

% set variables for plotting
subj = [2 6 5]; % select which participants to plot from each group
trials{1} = [1:10; 31:40; 331:340]; % trials to plot for 2-day group
trials{2} = [1:10; 31:40; 1251:1260]; % trials to plot for 5-day group
trials{3} = [1:10; 31:40; 2811:2820]; % trials to plot for 10-day group
groups = {'day2','day5','day10'};
titles = {'Baseline','Early','Late'};
groupNames = {'2-day','5-day','10-day'};
Ngroup = length(groupNames);
Nblock = length(titles);
col = [198 156 109
       0 146 69]./255; % colors for plotting

% generate plot
figure(1); clf
for k = 1:Ngroup % loop over groups
    s = subj(k);
    
    for j = 1:Nblock % loop over blocks
        subplot(3,3,(k-1)*3+j); hold on
        a = d.(groups{k}){s};
        
        for i = trials{k}(j,:) % loop over trials
            plot(a.targetAbs(i,1),a.targetAbs(i,2), '.', 'Color', [1 0.4 0.4], 'MarkerSize', 15); % plot targets
            plot(a.L{i}(:,1),a.L{i}(:,2),'Color',col(1,:)) % plot left hand
            plot(a.R{i}(:,1),a.R{i}(:,2),'Color',col(2,:)) % plot right hand
            plot(a.C{i}(:,1),a.C{i}(:,2),'k') % plot cursor
        end
        axis([0.1 0.95 -0.2 0.65])
        axis square
        if k == 1
            title(titles{j})
        end
        if j == 1
            ylabel(groupNames{k})
        end
    end
end
% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/traj','-dpdf','-painters')

%% plot reaches from flip block

% set variables for plotting
subj = [1 1 3]; % select which participants to analyze from each group
trialIdx = [474 467 1362 1356 2922 2920]; % select trials to plot
axisLims = [0.4 0.7 0.15 0.45;
            0.4 0.7 0.15 0.45
            0.35 0.65 0.1 0.4;
            0.45 0.75 0.05 0.35
            0.43 0.73 0.1 0.4;
            0.45 0.75 0.15 0.45]; % set axes for plots

% generate plot
figure(2); clf
for j = 1:Ngroup % loop over groups
    
    a = d.(groups{j}){subj(j)};
    for i = 1:2 % loop to generate plot for toward or away trial
        
        % set variables
        plotIdx = (j-1)*2 + i; % indices for plots
        trial = trialIdx(plotIdx); % trial to be plotted
        targ = a.targetAbs(trial-1:trial,:); % target trajectory
        curs = a.C{trial}; % cursor trajectory
        
        % calculate direction of instantaneous velocity vector
        init = a.init_x(trial)+20; % index of cursor position 150 ms after movement initiation
        vel = diff(curs); % compute velocity (not divided by time because only the direction is needed)
        vel = vel(init,:);
        angle = atan2(vel(2),vel(1)); % find direction of velocity vector 
        vector = 0.03*[cos(angle) sin(angle)]; % scale vector 
        
        subplot(2,3,(i-1)*3+j); hold on
        plot([targ(1,1) targ(1,1)], [targ(1,2)-0.15 targ(1,2)+0.15], 'k--', 'LineWidth', 1) % mirror axis
        plot(targ(1,1), targ(1,2), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 35) % starting target
        plot(targ(2,1), targ(2,2), '.', 'Color', [1 0.4 0.4], 'MarkerSize', 35) % ending target
        plot(targ(1,1) - (targ(2,1) - targ(1,1)), targ(2,2), 'o', 'Color', [1 0.4 0.4], 'MarkerSize', 10) % mirrored ending target
        plot(curs(:,1), curs(:,2),'k','LineWidth',1) % cursor trajectory
        plot([curs(init,1) curs(init,1)+vector(1)], [curs(init,2) curs(init,2)+vector(2)], 'Color', [0 0 1 0.5], 'LineWidth',4) % instantaneous velocity vector
        axis(axisLims(plotIdx,:))
        axis square
        if i == 1
            title(groupNames{j})
            if j == 1
                ylabel('Toward trial')
            end
        elseif i == 2 && j == 1
            ylabel('Away trial')
        end
    end
end

end
% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/traj_habit','-dpdf','-painters')