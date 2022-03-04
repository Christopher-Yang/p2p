% plots various kinematic measures as a function of trials

function plot_kinematics(data)

% set variables for analysis
rng(8);
groups = {'day2','day5','day10'}; % groups 
Ngroups = length(groups); % total number of groups
Nsubj = [length(data.day2) length(data.day5) length(data.day10)]; % number of subjects in each group
bin = 5; % number of trials to bin together for averaging

% colors for plotting
col = [180 180 0
       0 191 255
       255 99 71]./255;

for i = 1:Ngroups % loop over groups
    
    % store data from all subjects
    for j = 1:Nsubj(i)
        a = data.(groups{i}){j};
        pLength{i}(:,j) = a.pathlength;
        movtime{i}(:,j) = a.movtime;
        RT{i}(:,j) = a.RT;
        pkVel{i}(:,j) = a.pkVel;
    end
    
    % group data into bins of 5 trials
    Ntrials = size(RT{i},1);
    for j = 1:Ntrials/bin
        pLengthBin{i}(j,:) = mean(pLength{i}((j-1)*bin+1:(j-1)*bin+bin,:),1);
        movtimeBin{i}(j,:) = mean(movtime{i}((j-1)*bin+1:(j-1)*bin+bin,:),1);
        RTBin{i}(j,:) = mean(RT{i}((j-1)*bin+1:(j-1)*bin+bin,:),1,'omitnan');
        pkVelBin{i}(j,:) = mean(pkVel{i}((j-1)*bin+1:(j-1)*bin+bin,:),1);
    end
    
    % average binned data and compute standard error
    pLength_mu{i} = mean(pLengthBin{i},2);
    pLength_se{i} = std(pLengthBin{i},[],2)/sqrt(Nsubj(i));
    movtime_mu{i} = mean(movtimeBin{i},2);
    movtime_se{i} = std(movtimeBin{i},[],2)/sqrt(Nsubj(i));
    RT_mu{i} = mean(RTBin{i},2);
    RT_se{i} = std(RTBin{i},[],2)/sqrt(Nsubj(i));
    pkVel_mu{i} = mean(pkVelBin{i},2);
    pkVel_se{i} = std(pkVelBin{i},[],2)./sqrt(sum(~isnan(pkVelBin{i}),2));
end

% set variables for plotting
gblock = [3 2 1]; % set order in which to plot groups (10-day, 5-day, then 2-day)
trialIdx = [5 14 29]; % select which trials to plot from variable "trials"
lw = 0.25; % line width for plots
dayStart = [7 47 227 527]; % index of day 1, 2, 5, and 10
dayStartLabels = {'1','2','5','10'};
dayStart2 = 47:60:527; % index for all days

% x-axis for binned trials
trials{1} = 1:6;
for i = 1:29
    trials{i+1} = (i-1)*20 + 7:(i-1)*20 + 26;
end

%% Supplementary Figure 1A
f = figure(4); clf; hold on
set(f,'Position',[200 200 140 140]);

% vertical lines delineating days
plot([dayStart2; dayStart2],[0.1 0.5],'Color',[0.8 0.8 0.8])

% plot baseline data
for j = 1:3
    avg = mean(pkVel_mu{j}(trials{1}));
    plot([7 dayStart(j+1)+40], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
end

% plot data from days 1-2
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
xlabel('Day')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(0.1:0.2:0.5)
ylabel('Peak velocity (m/s)')
axis([7 566 0.1 0.5])
set(gca,'TickDir','out')

% prints figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/peak_vel','-dpdf','-painters')

% write data table for analysis in R
% 
% % vector of data for statistical analysis in R
% y = [mean(pkVel{1}(31:130,:))'; mean(pkVel{2}(331:430,:))'; mean(pkVel{3}(1231:1330,:))'; mean(pkVel{1}(end-199:end-100,:))'; mean(pkVel{2}(end-199:end-100,:))'; mean(pkVel{3}(end-199:end-100,:))'];
% 
% % labels data points in y to create R data frame
% groupNames([1:13 33:45],1) = "2-day";
% groupNames([14:27 46:59],1) = "5-day";
% groupNames([28:32 60:64],1) = "10-day";
% blockNames(1:32,1) = "before";
% blockNames(33:64,1) = "after";
% subject = [1:32 1:32]';
% T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','pkVel'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/pkVel.csv')

%% Supplementary Figure 1B

f = figure(5); clf; hold on
set(f,'Position',[200 200 140 140]);

% vertical lines delineating days
plot([dayStart2; dayStart2],[10 40],'Color',[0.8 0.8 0.8])

% plot baseline data
for j = 1:3
    avg = 100*mean(pLength_mu{j}(trials{1}));
    plot([7 dayStart(j+1)+40], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
end

% plot data from days 1-2
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
xlabel('Day')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(10:10:40)
ylabel('Path length (cm)')
axis([7 566 10 40])
set(gca,'TickDir','out')

% prints figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/path_length','-dpdf','-painters')

% write data table for analysis in R
% y = [mean(pLength{1}(31:130,:))'; mean(pLength{2}(331:430,:))'; mean(pLength{3}(1231:1330,:))'; mean(pLength{1}(end-199:end-100,:))'; mean(pLength{2}(end-199:end-100,:))'; mean(pLength{3}(end-199:end-100,:))'];
% T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','plength'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/path_length.csv')

%% Supplementary Figure 1C
f = figure(6); clf; hold on
set(f,'Position',[200 200 140 140]);

% vertical lines delineating days
plot([dayStart2; dayStart2],[0 6],'Color',[0.8 0.8 0.8])

% plot baseline data
for j = 1:3
    avg = mean(movtime_mu{j}(trials{1}));
    plot([7 dayStart(j+1)+40], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
end

% plot data from days 1-2
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
xlabel('Day')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(0:2:8)
ylabel('Movement time (s)')
axis([7 566 0 6])
set(gca,'TickDir','out')

% prints figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/move_time','-dpdf','-painters')

% write data table for analysis in R
% y = [mean(movtime{1}(31:130,:))'; mean(movtime{2}(331:430,:))'; mean(movtime{3}(1231:1330,:))'; mean(movtime{1}(end-199:end-100,:))'; mean(movtime{2}(end-199:end-100,:))'; mean(movtime{3}(end-199:end-100,:))'];
% T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','movtime'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/movtime.csv')

%% Supplementary Figure 1D
f = figure(7); clf; hold on 
set(f,'Position',[200 200 140 140]);

% vertical lines delineating days
plot([dayStart2; dayStart2],[0 2],'Color',[0.8 0.8 0.8])

% plot baseline data
for j = 1:3
    avg = mean(RT_mu{j}(trials{1}));
    plot([7 dayStart(j+1)+40], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
end


% plot data from days 1-2
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
xlabel('Day')
xticks(dayStart)
xticklabels(dayStartLabels)
yticks(.5:.5:2)
ylabel('Reaction time (s)')
axis([7 566 .2 2])
set(gca,'TickDir','out')

% prints figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/rt','-dpdf','-painters')

% write data table for analysis in R
% y = [mean(RT{1}(31:130,:))'; mean(RT{2}(331:430,:))'; mean(RT{3}(1231:1330,:))'; mean(RT{1}(end-199:end-100,:))'; mean(RT{2}(end-199:end-100,:))'; mean(RT{3}(end-199:end-100,:))'];
% T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','RT'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/RT.csv')

end
