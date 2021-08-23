function plot_kinematics(d)

% set variables for analysis
rng(8);
groups = {'day2','day5','day10'}; % groups 
Ngroups = length(groups); % total number of groups
Nsubj = [length(d.day2) length(d.day5) length(d.day10)]; % number of subjects in each group
bin = 5; % number of trials to bin together for averaging
habitBlocks = [431:530; 1331:1430; 2831:2930]; % trial indices of habit blocks for each group

% colors for plotting
col = [180 180 0
       0 191 255
       255 99 71]./255;

for i = 1:Ngroups
    RT_away{i} = [];
    RT_toward{i} = [];
    RTX_away{i} = [];
    RTX_toward{i} = [];
    RTY_away{i} = [];
    RTY_toward{i} = [];
    vel_away{i} = [];
    vel_toward{i} = [];
    velX_away{i} = [];
    velX_toward{i} = [];
    velY_away{i} = [];
    velY_toward{i} = [];
    
    for j = 1:Nsubj(i)
        a = d.(groups{i}){j};
        pLength{i}(:,j) = a.pathlength;
        movtime{i}(:,j) = a.movtime./1000;
        RT{i}(:,j) = a.RT;
        RTX{i}(:,j) = a.RT_x;
        RTY{i}(:,j) = a.RT_y;
        vel{i}(:,j) = a.initVel_filt;
        pkVel{i}(:,j) = a.pkVel;
        velX{i}(:,j) = abs(a.initVel_x);
        velY{i}(:,j) = abs(a.initVel_y);
        
        % need to use == 1/0 because some values are NaN
        awayIdx = a.incorrectReach_x(end-99:end) == 1;
        towardIdx = a.incorrectReach_x(end-99:end) == 0;

        last_RT = RT{i}(habitBlocks(i,:),j);
        last_RTX = RTX{i}(habitBlocks(i,:),j);
        last_RTY = RTY{i}(habitBlocks(i,:),j);
        last_vel = vel{i}(habitBlocks(i,:),j);
        last_velX = velX{i}(habitBlocks(i,:),j);
        last_velY = velY{i}(habitBlocks(i,:),j);
        
        RT_away{i}(j) = mean(last_RT(awayIdx));
        RT_toward{i}(j) = mean(last_RT(towardIdx));
        RTX_away{i}(j) = nanmean(last_RTX(awayIdx));
        RTX_toward{i}(j) = nanmean(last_RTX(towardIdx));
        RTY_away{i}(j) = nanmean(last_RTY(awayIdx));
        RTY_toward{i}(j) = nanmean(last_RTY(towardIdx));
        vel_away{i}(j) = nanmean(last_vel(awayIdx));
        vel_toward{i}(j) = nanmean(last_vel(towardIdx));
        velX_away{i}(j) = nanmean(last_velX(awayIdx));
        velX_toward{i}(j) = nanmean(last_velX(towardIdx));
        velY_away{i}(j) = nanmean(last_velY(awayIdx));
        velY_toward{i}(j) = nanmean(last_velY(towardIdx));
        
    end
    
    Ntrials = size(RT{i},1);
    for j = 1:Ntrials/bin
        pLengthBin{i}(j,:) = mean(pLength{i}((j-1)*5+1:(j-1)*5+5,:),1);
        movtimeBin{i}(j,:) = mean(movtime{i}((j-1)*5+1:(j-1)*5+5,:),1);
        RTBin{i}(j,:) = mean(RT{i}((j-1)*5+1:(j-1)*5+5,:),1);
        velBin{i}(j,:) = nanmean(vel{i}((j-1)*5+1:(j-1)*5+5,:),1);
        pkVelBin{i}(j,:) = nanmean(pkVel{i}((j-1)*5+1:(j-1)*5+5,:),1);
    end
    
    pLength_mu{i} = mean(pLengthBin{i},2);
    pLength_se{i} = std(pLengthBin{i},[],2)/sqrt(Nsubj(i));
    movtime_mu{i} = mean(movtimeBin{i},2);
    movtime_se{i} = std(movtimeBin{i},[],2)/sqrt(Nsubj(i));
    RT_mu{i} = mean(RTBin{i},2);
    RT_se{i} = std(RTBin{i},[],2)/sqrt(Nsubj(i));
    vel_mu{i} = nanmean(velBin{i},2);
    vel_se{i} = nanstd(velBin{i},[],2)./sqrt(sum(~isnan(velBin{i}),2));
    pkVel_mu{i} = nanmean(pkVelBin{i},2);
    pkVel_se{i} = nanstd(pkVelBin{i},[],2)./sqrt(sum(~isnan(pkVelBin{i}),2));
end

%% plot path length

% set variables for plotting
gblock = [3 2 1]; % set order in which to plot groups (10-day, 5-day, then 2-day)
trialIdx = [5 14 29]; % select which trials to plot from variable "trials"
lw = 0.25; % line width for plots
% dayStart = [7 47:60:527]; % xticks
% dayStartLabels = [31 231:300:2930]; % xticklabels
dayStart = [7 47 227 527];
dayStartLabels = {'1','2','5','10'};

% x-axis for binned trials
trials{1} = 1:6;
for i = 1:29
    trials{i+1} = (i-1)*20 + 7:(i-1)*20 + 26;
end

f = figure(1); clf
set(f,'Position',[200 200 170 120]);

% plot baseline data
subplot(1,8,1:2); hold on
for j = 1:3
    s = shadedErrorBar(trials{1}, pLength_mu{gblock(j)}(trials{1})*100, pLength_se{gblock(j)}(trials{1})*100);
    editErrorBar(s,col(gblock(j),:),lw);
end
ylabel('Path length (cm)')
xticks(1)
xticklabels('Baseline')
yticks(10:10:40)
axis([1 6 10 40])
set(gca,'TickDir','out')

% plot data from days 1-2
subplot(1,8,3:8); hold on

for i = 2:5
    for j = 1:3
        s = shadedErrorBar(trials{i},pLength_mu{gblock(j)}(trials{i})*100, pLength_se{gblock(j)}(trials{i})*100);
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 3-5
for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i}, pLength_mu{gblock(j)}(trials{i})*100, pLength_se{gblock(j)}(trials{i})*100);
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 6-10
for i = 15:29
    s = shadedErrorBar(trials{i}, pLength_mu{3}(trials{i})*100, pLength_se{3}(trials{i})*100);
    editErrorBar(s,col(3,:),lw);
end
xlabel('Trial Number')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(10:10:40)
axis([7 566 10 40])
set(gca,'TickDir','out')

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/path_length','-dpdf','-painters')

y = [mean(pLength{1}(1:30,:))'; mean(pLength{2}(1:30,:))'; mean(pLength{3}(1:30,:))'; mean(pLength{1}(end-199:end-100,:))'; mean(pLength{2}(end-199:end-100,:))'; mean(pLength{3}(end-199:end-100,:))'];

groupNames([1:13 33:45],1) = "2-day";
groupNames([14:27 46:59],1) = "5-day";
groupNames([28:32 60:64],1) = "10-day";
blockNames(1:32,1) = "baseline";
blockNames(33:64,1) = "late";
subject = [1:32 1:32]';
T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','plength'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/path_length.csv')

%% plot movement time
f = figure(2); clf
set(f,'Position',[200 200 170 120]);

% plot baseline data
subplot(1,8,1:2); hold on
for j = 1:3
    s = shadedErrorBar(trials{1},movtime_mu{gblock(j)}(trials{1}), movtime_se{gblock(j)}(trials{1}));
    editErrorBar(s,col(gblock(j),:),lw);
end
xticks(1)
xticklabels([])
yticks(0:2:8)
axis([1 6 0 6])
ylabel('Movement time (s)')
set(gca,'TickDir','out')

% plot data from days 1-2
subplot(1,8,3:8); hold on

for i = 2:5
    for j = 1:3
        s = shadedErrorBar(trials{i},movtime_mu{gblock(j)}(trials{i}), movtime_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 3-5
for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},movtime_mu{gblock(j)}(trials{i}), movtime_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 6-10
for i = 15:29
    s = shadedErrorBar(trials{i},movtime_mu{3}(trials{i}), movtime_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end
xlabel('Trial Number')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(0:2:8)
axis([7 566 0 6])
set(gca,'TickDir','out')

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/move_time','-dpdf','-painters')

y = [mean(movtime{1}(1:30,:))'; mean(movtime{2}(1:30,:))'; mean(movtime{3}(1:30,:))'; mean(movtime{1}(end-199:end-100,:))'; mean(movtime{2}(end-199:end-100,:))'; mean(movtime{3}(end-199:end-100,:))'];

groupNames([1:13 33:45],1) = "2-day";
groupNames([14:27 46:59],1) = "5-day";
groupNames([28:32 60:64],1) = "10-day";
blockNames(1:32,1) = "baseline";
blockNames(33:64,1) = "late";
subject = [1:32 1:32]';
T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','movtime'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/movtime.csv')

%% plot reaction time
f = figure(3); clf
set(f,'Position',[200 200 170 120]);

% plot baseline data
subplot(1,8,1:2); hold on
for j = 1:3
    s = shadedErrorBar(trials{1},RT_mu{gblock(j)}(trials{1}), RT_se{gblock(j)}(trials{1}));
    editErrorBar(s,col(gblock(j),:),lw);
end
xticks(1)
xticklabels([])
yticks(300:400:1700)
axis([1 6 300 1100])
ylabel('Reaction time (ms)')
set(gca,'TickDir','out')

% plot data from days 1-2
subplot(1,8,3:8); hold on
for i = 2:5
    for j = 1:3
        s = shadedErrorBar(trials{i},RT_mu{gblock(j)}(trials{i}), RT_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 3-5
for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},RT_mu{gblock(j)}(trials{i}), RT_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 6-10
for i = 15:29
    s = shadedErrorBar(trials{i},RT_mu{3}(trials{i}), RT_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end

xlabel('Trial Number')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(300:400:1700)
axis([7 566 300 1100])
set(gca,'TickDir','out')

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/rt','-dpdf','-painters')

y = [mean(RT{1}(1:30,:))'; mean(RT{2}(1:30,:))'; mean(RT{3}(1:30,:))'; mean(RT{1}(end-199:end-100,:))'; mean(RT{2}(end-199:end-100,:))'; mean(RT{3}(end-199:end-100,:))'];

groupNames([1:13 33:45],1) = "2-day";
groupNames([14:27 46:59],1) = "5-day";
groupNames([28:32 60:64],1) = "10-day";
blockNames(1:32,1) = "baseline";
blockNames(33:64,1) = "late";
subject = [1:32 1:32]';
T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','RT'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/RT.csv')

%% plot peak reach velocity
f = figure(4); clf
set(f,'Position',[200 200 170 120]);

% plot baseline data
subplot(1,8,1:2); hold on
for j = 1:3
    s = shadedErrorBar(trials{1},pkVel_mu{gblock(j)}(trials{1}), pkVel_se{gblock(j)}(trials{1}));
    editErrorBar(s,col(gblock(j),:),lw);
end
xticks(1)
xticklabels([])
yticks(0.1:0.2:0.5)
axis([1 6 0.1 0.5])
ylabel('Peak velocity (m/s)')
set(gca,'TickDir','out')

% plot data from days 1-2
subplot(1,8,3:8); hold on

for i = 2:5
    for j = 1:3
        s = shadedErrorBar(trials{i},pkVel_mu{gblock(j)}(trials{i}), pkVel_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 3-5
for i = 6:14
    for j = 1:2
        s = shadedErrorBar(trials{i},pkVel_mu{gblock(j)}(trials{i}), pkVel_se{gblock(j)}(trials{i}));
        editErrorBar(s,col(gblock(j),:),lw);
    end
end

% plot data from days 6-10
for i = 15:29
    s = shadedErrorBar(trials{i},pkVel_mu{3}(trials{i}), pkVel_se{3}(trials{i}));
    editErrorBar(s,col(3,:),lw);
end

xticks(dayStart)
xticklabels(dayStartLabels)
yticks(0.1:0.2:0.5)
axis([7 566 0.1 0.5])
set(gca,'TickDir','out')

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/peak_vel','-dpdf','-painters')

y = [mean(pkVel{1}(1:30,:))'; mean(pkVel{2}(1:30,:))'; mean(pkVel{3}(1:30,:))'; mean(pkVel{1}(end-199:end-100,:))'; mean(pkVel{2}(end-199:end-100,:))'; mean(pkVel{3}(end-199:end-100,:))'];

groupNames([1:13 33:45],1) = "2-day";
groupNames([14:27 46:59],1) = "5-day";
groupNames([28:32 60:64],1) = "10-day";
blockNames(1:32,1) = "baseline";
blockNames(33:64,1) = "late";
subject = [1:32 1:32]';
T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','pkVel'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/pkVel.csv')

%% compare path length, movement time, reaction time from flip block

for i = 1:3
    pLength3{i} = mean(pLength{i}(habitBlocks(i,:),:),1);
    movtime3{i} = mean(movtime{i}(habitBlocks(i,:),:),1);
    RT3{i} = mean(RT{i}(habitBlocks(i,:),:),1);
end

figure(5); clf
subplot(1,3,1); hold on
for i = 1:3
    plot(i + 0.5*(rand(1,Nsubj(i)) - 0.5), pLength3{i}*100, '.', 'Color', col(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(i, mean(pLength3{i})*100, 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10)
end
axis([0.5 3.5 10 45])
xticks([])
yticks(10:10:40)
ylabel('Path length (cm)')
set(gca, 'TickDir', 'out')
legend({'2-day','5-day','10-day'})

subplot(1,3,2); hold on
for i = 1:3
    plot(i + 0.5*(rand(1,Nsubj(i)) - 0.5), movtime3{i}, '.', 'Color', col(i,:), 'MarkerSize', 20)
    plot(i, mean(movtime3{i}), 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10)
end
axis([0.5 3.5 1 4.5])
xticks([])
yticks(1:4)
title('Data from flip block')
ylabel('Movement time (s)')
set(gca, 'TickDir', 'out')

subplot(1,3,3); hold on
for i = 1:3
    plot(i + 0.5*(rand(1,Nsubj(i)) - 0.5), RT3{i}, '.', 'Color', col(i,:), 'MarkerSize', 20)
    plot(i, mean(RT3{i}), 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10)
end
axis([0.5 3.5 350 1050])
xticks([])
yticks(400:200:1000)
ylabel('Reaction time (ms)')
set(gca, 'TickDir', 'out')

%% plot kinematics across blocks with subjects connected by lines

for i = 1:3
    pLength3{i}(1,:) = mean(pLength{i}(1:30,:),1);
    pLength3{i}(2,:) = mean(pLength{i}(habitBlocks(i,:)-100,:),1);
    pLength3{i}(3,:) = mean(pLength{i}(habitBlocks(i,:),:),1);

    movtime3{i}(1,:) = mean(movtime{i}(1:30,:),1);
    movtime3{i}(2,:) = mean(movtime{i}(habitBlocks(i,:)-100,:),1);
    movtime3{i}(3,:) = mean(movtime{i}(habitBlocks(i,:),:),1);

    RT3{i}(1,:) = mean(RT{i}(1:30,:),1);
    RT3{i}(2,:) = mean(RT{i}(habitBlocks(i,:)-100,:),1);
    RT3{i}(3,:) = mean(RT{i}(habitBlocks(i,:),:),1);

    vel3{i}(1,:) = nanmean(vel{i}(1:30,:),1);
    vel3{i}(2,:) = nanmean(vel{i}(habitBlocks(i,:)-100,:),1);
    vel3{i}(3,:) = nanmean(vel{i}(habitBlocks(i,:),:),1);
end

subject = [repelem({'2-day'},13) repelem({'5-day'},14) repelem({'10-day'},5)]';
meas(:,1) = [vel3{1}(2,:) vel3{2}(2,:) vel3{3}(2,:)];
meas(:,2) = [vel3{1}(3,:) vel3{2}(3,:) vel3{3}(3,:)];

t = table(subject, meas(:,1), meas(:,2), 'VariableNames', {'Subject','meas1','meas2'});
Meas = table([1 2]', 'VariableNames', {'Block'});

rm = fitrm(t, 'meas1-meas2~Subject', 'WithinDesign', Meas);


% velocity = [reshape(vel3{1}, [numel(vel3{1}) 1]); reshape(vel3{2}, [numel(vel3{2}) 1]); ...
%     reshape(vel3{3}, [numel(vel3{3}) 1])];
% dlmwrite('C:/Users/Chris/Documents/R/habit/data/velocity.csv', velocity)

figure(6); clf
subplot(2,2,1); hold on
for i = 1:3
    plot(3*(i-1) + (1:3), pLength3{i}*100, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(3*(i-1) + (1:3), mean(pLength3{i},2)*100, 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
end
xticks(1:9)
xlim([0.5 9.5])
xticklabels({'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip'})
ylabel('Path length (cm)')
set(gca,'TickDir','out')
legend({'2-day','5-day','10-day'})

subplot(2,2,2); hold on
for i = 1:3
    plot(3*(i-1) + (1:3), movtime3{i}, 'Color', col(i,:))
    plot(3*(i-1) + (1:3), mean(movtime3{i},2), 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
end
xticks(1:9)
xlim([0.5 9.5])
xticklabels({'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip'})
ylabel('Movement time (s)')
set(gca,'TickDir','out')

subplot(2,2,3); hold on
for i = 1:3
    plot(3*(i-1) + (1:3), RT3{i}, 'Color', col(i,:))
    plot(3*(i-1) + (1:3), mean(RT3{i},2), 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
end
xticks(1:9)
xlim([0.5 9.5])
xticklabels({'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip'})
ylabel('Reaction time (ms)')
set(gca,'TickDir','out')

subplot(2,2,4); hold on
for i = 1:3
    plot(3*(i-1) + (1:3), vel3{i}, 'Color', col(i,:))
    plot(3*(i-1) + (1:3), mean(vel3{i},2), 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
end
xticks(1:9)
xlim([0.5 9.5])
xticklabels({'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip', 'Baseline', 'Late', 'Flip'})
ylabel('Initial reach velocity (m/s)')
set(gca,'TickDir','out')

%% plot reaction time and initial velocity of away vs toward trials

figure(7); clf
subplot(1,2,1); hold on
for i = 1:Ngroups
    n = length(RT_toward{i});
    
    % plot away trials
    plot(i + 0.5*(rand(1,n) - 0.5), RT_away{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility','off')
    plot(i, mean(RT_away{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:))
    
    % plot toward trials
    plot(i+4 + 0.5*(rand(1,n) - 0.5),RT_toward{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(i+4, mean(RT_toward{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:), 'HandleVisibility', 'off')
end    
xticks([2 6])
xticklabels({'Away','Toward'})
ylim([300 1100])
yticks(300:200:1100)
ylabel('Reaction time (ms)')
set(gca, 'TickDir', 'out')
legend({'2-day','5-day','10-day'},'Location','northwest')

subplot(1,2,2); hold on
for i = 1:Ngroups
    n = length(vel_toward{i});
    
    % plot away trials
    plot(i + 0.5*(rand(1,n) - 0.5), vel_away{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility','off')
    plot(i, mean(vel_away{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:))
    
    % plot toward trials
    plot(i+4 + 0.5*(rand(1,n) - 0.5),vel_toward{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(i+4, mean(vel_toward{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:), 'HandleVisibility', 'off')
end
xticks([2 6])
xticklabels({'Away','Toward'})
ylim([0.05 0.35])
ylabel('Tangential initial velocity (m/s)')
set(gca, 'TickDir', 'out')


f = figure(15); clf
set(f,'Position',[200 200 140 120]); hold on
for i = 1:Ngroups
    n = length(RTX_toward{i});
    
    % plot away trials
    plot(i + 0.5*(rand(1,n) - 0.5), RTX_away{i}, '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility','off')
    plot(i, mean(RTX_away{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', col(i,:))
    
    % plot toward trials
    plot(i+4 + 0.5*(rand(1,n) - 0.5),RTX_toward{i}, '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot(i+4, mean(RTX_toward{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'HandleVisibility', 'off')
end    
xticks([2 6])
xticklabels({'Away','Toward'})
yticks(200:400:1400)
axis([0.5 7.5 200 1400])
ylabel('Horizontal reaction time (ms)')
set(gca, 'TickDir', 'out')

% subplot(1,4,2); hold on
% for i = 1:Ngroups
%     n = length(RTY_toward{i});
%     
%     % plot away trials
%     plot(i + 0.5*(rand(1,n) - 0.5), RTY_away{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility','off')
%     plot(i, mean(RTY_away{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:))
%     
%     % plot toward trials
%     plot(i+4 + 0.5*(rand(1,n) - 0.5),RTY_toward{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility', 'off')
%     plot(i+4, mean(RTY_toward{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:), 'HandleVisibility', 'off')
% end    
% xticks([2 6])
% xticklabels({'Away','Toward'})
% ylim([200 1800])
% ylabel('Vertical reaction time (ms)')
% set(gca, 'TickDir', 'out')

% subplot(1,2,2); hold on
% for i = 1:Ngroups
%     n = length(vel_toward{i});
%     
%     % plot away trials
%     plot(i + 0.5*(rand(1,n) - 0.5), velX_away{i}, '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility','off')
%     plot(i, mean(velX_away{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', col(i,:))
%     
%     % plot toward trials
%     plot(i+4 + 0.5*(rand(1,n) - 0.5),velX_toward{i}, '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility', 'off')
%     plot(i+4, mean(velX_toward{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'HandleVisibility', 'off')
% end
% xticks([2 6])
% xticklabels({'Away','Toward'})
% yticks(0:0.1:0.3)
% axis([0.5 7.5 0 0.3])
% ylabel('Horizontal initial velocity (m/s)')
% set(gca, 'TickDir', 'out')

% subplot(1,4,4); hold on
% for i = 1:Ngroups
%     n = length(vel_toward{i});
%     
%     % plot away trials
%     plot(i + 0.5*(rand(1,n) - 0.5), velY_away{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility','off')
%     plot(i, mean(velY_away{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:))
%     
%     % plot toward trials
%     plot(i+4 + 0.5*(rand(1,n) - 0.5),velY_toward{i}, '.', 'MarkerSize', 20, 'Color', col(i,:), 'HandleVisibility', 'off')
%     plot(i+4, mean(velY_toward{i}), 'ko', 'LineWidth', 1, 'MarkerSize', 10, 'MarkerFaceColor', col(i,:), 'HandleVisibility', 'off')
% end
% xticks([2 6])
% xticklabels({'Away','Toward'})
% yticks(0:0.1:0.3)
% ylim([0 0.3])
% ylabel('Vertical initial velocity (m/s)')
% set(gca, 'TickDir', 'out')

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/horizontal_kinematics','-dpdf','-painters')

%% correlate reaction time and initial velocity to proportion of away trials

habit = NaN(14,Ngroups);
% store proportion of away trials in "habit"
for i = 1:Ngroups
    n = length(d.(groups{i}));
    
    for j = 1:n
        a = d.(groups{i}){j}.incorrectReach_x(habitBlocks(i,:));
        num = nansum(a);
        den = 100-sum(isnan(a));
        habit(j,i) = 100*num/den;
    end
end

% average reaction times within subjects
for i = 1:3
    pLength3{i} = pLength{i}(habitBlocks(i,:)-100,:);
    pLength3{i} = mean(pLength3{i},1);
    
    movtime3{i} = movtime{i}(habitBlocks(i,:)-100,:);
    movtime3{i} = mean(movtime3{i},1);
    
    RT3{i} = RT{i}(habitBlocks(i,:)-100,:); 
    RT3{i} = mean(RT3{i},1);
    
    vel3{i} = vel{i}(habitBlocks(i,:)-100,:); 
    vel3{i} = nanmean(vel3{i},1);
end

% plot data
figure(8); clf
subplot(1,4,1); hold on
for i = 1:3
    plot(RT3{i}, habit(~isnan(habit(:,i)),i), '.', 'Color', col(i,:), 'MarkerSize', 20)
end

% plot best-fit lines from least-squares regression
idx = [200 700];
for i = 1:3
    p = polyfit(RT3{i}', habit(~isnan(habit(:,i)),i), 1);
    plot(idx, p(1)*idx + p(2), 'Color', col(i,:))
end

% xlim([300 1100])
ylim([10 60])
xlabel('Reaction time (ms)')
ylabel('Proportion of away trials (%)')
set(gca, 'TickDir', 'out')

subplot(1,4,2); hold on
for i = 1:3
    plot(vel3{i}, habit(~isnan(habit(:,i)),i), '.', 'Color', col(i,:), 'MarkerSize', 20)
end

% plot best-fit lines from least-squares regression
idx = [0 0.6];
for i = 1:3
    p = polyfit(vel3{i}', habit(~isnan(habit(:,i)),i), 1);
    plot(idx, p(1)*idx + p(2), 'Color', col(i,:))
end

% xlim([0.05 0.35])
ylim([10 60])
xlabel('Velocity (m/s)')
ylabel('Proportion of away trials (%)')
set(gca, 'TickDir', 'out')
legend({'2-day','5-day','10-day'}, 'Location', 'northeast')

subplot(1,4,3); hold on
for i = 1:3
    plot(pLength3{i}*100, habit(~isnan(habit(:,i)),i), '.', 'Color', col(i,:), 'MarkerSize', 20)
end

% plot best-fit lines from least-squares regression
idx = [10 40];
for i = 1:3
    p = polyfit(pLength3{i}'*100, habit(~isnan(habit(:,i)),i), 1);
    plot(idx, p(1)*idx + p(2), 'Color', col(i,:))
end

ylim([10 60])
xlabel('Path length (cm)')
ylabel('Proportion of away trials (%)')
set(gca, 'TickDir', 'out')

subplot(1,4,4); hold on
for i = 1:3
    plot(movtime3{i}, habit(~isnan(habit(:,i)),i), '.', 'Color', col(i,:), 'MarkerSize', 20)
end

% plot best-fit lines from least-squares regression
idx = [0.5 2.5];
for i = 1:3
    p = polyfit(movtime3{i}', habit(~isnan(habit(:,i)),i), 1);
    plot(idx, p(1)*idx + p(2), 'Color', col(i,:))
end

ylim([10 60])
xlabel('Movement time (s)')
ylabel('Proportion of away trials (%)')
set(gca, 'TickDir', 'out')
legend({'2-day','5-day','10-day'}, 'Location', 'northeast')
%% bin reaches by the target's location

% store data in variables
for k = 1:2 % either plots late learning (k=1) or post-learning (k=2)
    for i = 1:Ngroups % loop over groups
        for j = 1:4 % loop over target direction bins
            
            % for group i, get index of first trial for a block
            if k == 1
                trialIdx = habitBlocks(i,1)-100; % indices for late learning
            else
                trialIdx = habitBlocks(i,1); % indices for post-learning
            end
            
            % extract trials in direction bin j
            idx = find(d.(groups{i}){1}.targBin == j);
            idx = idx(idx >= trialIdx);
            
            % average data within subjects
            pLength2{k}{i}(:,j) = mean(pLength{i}(idx,:),1);
            movtime2{k}{i}(:,j) = mean(movtime{i}(idx,:),1);
            RT2{k}{i}(:,j) = mean(RT{i}(idx,:),1);
        end
        
        % average data across subjects
        pLength2_mu{k}{i} = mean(pLength2{k}{i},1);
        movtime2_mu{k}{i} = mean(movtime2{k}{i},1);
        RT2_mu{k}{i} = mean(RT2{k}{i},1);
    end
end

figure(9); clf
for k = 1:2
    
    % plot path length data
    subplot(2,3,3*(k-1)+1); hold on
    for i = 1:Ngroups
        for j = 1:4
            plot(6*(i-1)+j + 0.3*(rand(Nsubj(i),1)-0.5),pLength2{k}{i}(:,j)*100,'.','Color',col(i,:),'MarkerSize',20)
        end
        plot(6*(i-1)+(1:4), pLength2_mu{k}{i}*100, 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
    end
    if k == 1
        xticks([])
        title('Late')
    else
        xticks([1 4 7 10 13 16])
        xticklabels({'Close','Far','Close','Far','Close','Far'})
        title('Post')
    end
    ylabel('Path length (cm)')
    axis([0 17 10 55])
    yticks(10:15:55)
    set(gca,'TickDir','out')
    
    % plot movement time data
    subplot(2,3,3*(k-1)+2); hold on
    for i = 1:Ngroups
        for j = 1:4
            plot(6*(i-1)+j + 0.3*(rand(Nsubj(i),1)-0.5),movtime2{k}{i}(:,j),'.','Color',col(i,:),'MarkerSize',20)
        end
        plot(6*(i-1)+(1:4), movtime2_mu{k}{i}, 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
    end
    if k == 1
        xticks([])
        title('Late')
    else
        xticks([1 4 7 10 13 16])
        xticklabels({'Close','Far','Close','Far','Close','Far'})
        title('Post')
    end
    ylabel('Movement time (s)')
    axis([0 17 0.5 5])
    yticks(1:5)
    set(gca,'TickDir','out')
    
    % plot reaction time data
    subplot(2,3,3*(k-1)+3); hold on
    for i = 1:Ngroups
        for j = 1:4
            plot(6*(i-1)+j + 0.3*(rand(Nsubj(i),1)-0.5),RT2{k}{i}(:,j),'.','Color',col(i,:),'MarkerSize',20,'HandleVisibility','off')
        end
        plot(6*(i-1)+(1:4), RT2_mu{k}{i}, 'ko', 'MarkerFaceColor', col(i,:), 'MarkerSize', 10, 'LineWidth', 1)
    end
    if k == 1
        xticks([])
        title('Late')
        legend(groups,'Location','northwest')
    else
        xticks([1 4 7 10 13 16])
        xticklabels({'Close','Far','Close','Far','Close','Far'})
        title('Post')
    end
    ylabel('Reaction Time (ms)')
    axis([0 17 300 1300])
    yticks(300:300:1200)
    set(gca,'TickDir','out')
end
