% plots various kinematic measures as a function of trials

function plot_kinematics(data)

% set variables for analysis
rng(8);
Nsubj = length(data); % number of subjects in each group
bin = 5; % number of trials to bin together for averaging

% colors for plotting
col = [180 180 0
       0 191 255
       255 99 71]./255;

% store data from all subjects
for j = 1:Nsubj
    a = data{j};
    pLength(:,j) = a.pathlength;
    movtime(:,j) = a.movtime;
    RT(:,j) = a.RT;
    pkVel(:,j) = a.pkVel;
end

% group data into bins of 5 trials
Ntrials = size(RT,1);
for j = 1:Ntrials/bin
    pLengthBin(j,:) = mean(pLength((j-1)*bin+1:(j-1)*bin+bin,:),1);
    movtimeBin(j,:) = mean(movtime((j-1)*bin+1:(j-1)*bin+bin,:),1);
    RTBin(j,:) = mean(RT((j-1)*bin+1:(j-1)*bin+bin,:),1,'omitnan');
    pkVelBin(j,:) = mean(pkVel((j-1)*bin+1:(j-1)*bin+bin,:),1);
end

% average binned data and compute standard error
pLength_mu = mean(pLengthBin,2);
pLength_se = std(pLengthBin,[],2)/sqrt(Nsubj);
movtime_mu = mean(movtimeBin,2);
movtime_se = std(movtimeBin,[],2)/sqrt(Nsubj);
RT_mu = mean(RTBin,2);
RT_se = std(RTBin,[],2)/sqrt(Nsubj);
pkVel_mu = mean(pkVelBin,2);
pkVel_se = std(pkVelBin,[],2)./sqrt(sum(~isnan(pkVelBin),2));

% set variables for plotting
gblock = [3 2 1]; % set order in which to plot groups (10-day, 5-day, then 2-day)
lw = 0.25; % line width for plots
dayStart = [7 47 227 527]; % index of day 1, 2, 5, and 10
dayStartLabels = {'1','2','5','10'};
dayStart2 = 47:60:527; % index for all days

% x-axis for binned trials
% trials{1} = 1:6;
% for i = 1:29
%     trials{i+1} = (i-1)*20 + 7:(i-1)*20 + 26;
% end

%%

plot(pLength_mu)

%% Figure 2D
% f = figure(12); clf; hold on
% set(f,'Position',[200 200 140 140]);
% 
% % vertical lines delineating days
% plot([dayStart2; dayStart2],[0.1 0.5],'Color',[0.8 0.8 0.8])
% 
% % delineate flip block
% rectangle('Position',[87 0.1 20 0.4], 'EdgeColor', 'none', 'FaceColor', [col(1,:) 0.1])
% rectangle('Position',[267 0.1 20 0.4], 'EdgeColor', 'none', 'FaceColor', [col(2,:) 0.1])
% rectangle('Position',[567 0.1 20 0.4], 'EdgeColor', 'none', 'FaceColor', [col(3,:) 0.1])
% 
% % plot baseline data
% for j = 1:3
%     avg = mean(pkVel_mu{j}(trials{1}));
%     plot([7 dayStart(j+1)+60], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
% end
% 
% % plot data from days 1-2
% for i = 2:6
%     for j = 1:3
%         s = shadedErrorBar(trials{i},pkVel_mu{gblock(j)}(trials{i}), pkVel_se{gblock(j)}(trials{i}));
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 3-5
% for i = 6:15
%     for j = 1:2
%         s = shadedErrorBar(trials{i},pkVel_mu{gblock(j)}(trials{i}), pkVel_se{gblock(j)}(trials{i}));
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 6-10
% for i = 15:30
%     s = shadedErrorBar(trials{i},pkVel_mu{3}(trials{i}), pkVel_se{3}(trials{i}));
%     editErrorBar(s,col(3,:),lw);
% end
% xlabel('Day')
% xticks(dayStart)
% xticklabels(dayStartLabels)
% yticks(0.1:0.2:0.5)
% ylabel('Peak velocity (m/s)')
% axis([7 586 0.1 0.5])
% set(gca,'TickDir','out')
% 
% % prints figure for Illustrator
% % print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/peak_vel','-dpdf','-painters')
% 
% %% Figure 2E
% 
% f = figure(13); clf; hold on
% set(f,'Position',[200 200 140 140]);
% 
% % vertical lines delineating days
% plot([dayStart2; dayStart2],[10 40],'Color',[0.8 0.8 0.8])
% 
% % delineate flip block
% rectangle('Position',[87 10 20 30], 'EdgeColor', 'none', 'FaceColor', [col(1,:) 0.1])
% rectangle('Position',[267 10 20 30], 'EdgeColor', 'none', 'FaceColor', [col(2,:) 0.1])
% rectangle('Position',[567 10 20 30], 'EdgeColor', 'none', 'FaceColor', [col(3,:) 0.1])
% 
% % plot baseline data
% for j = 1:3
%     avg = 100*mean(pLength_mu{j}(trials{1}));
%     plot([7 dayStart(j+1)+60], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
% end
% 
% % plot data from days 1-2
% for i = 2:6
%     for j = 1:3
%         s = shadedErrorBar(trials{i},pLength_mu{gblock(j)}(trials{i})*100, pLength_se{gblock(j)}(trials{i})*100);
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 3-5
% for i = 6:15
%     for j = 1:2
%         s = shadedErrorBar(trials{i}, pLength_mu{gblock(j)}(trials{i})*100, pLength_se{gblock(j)}(trials{i})*100);
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 6-10
% for i = 15:30
%     s = shadedErrorBar(trials{i}, pLength_mu{3}(trials{i})*100, pLength_se{3}(trials{i})*100);
%     editErrorBar(s,col(3,:),lw);
% end
% xlabel('Day')
% xticks(dayStart)
% xticklabels(dayStartLabels)
% yticks(10:10:40)
% ylabel('Path length (cm)')
% axis([7 586 10 40])
% set(gca,'TickDir','out')
% 
% % prints figure for Illustrator
% % print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/path_length','-dpdf','-painters')
% 
% %% Figure 2F
% f = figure(14); clf; hold on
% set(f,'Position',[200 200 140 140]);
% 
% % vertical lines delineating days
% plot([dayStart2; dayStart2],[0 6],'Color',[0.8 0.8 0.8])
% 
% % delineate flip block
% rectangle('Position',[87 0 20 6], 'EdgeColor', 'none', 'FaceColor', [col(1,:) 0.1])
% rectangle('Position',[267 0 20 6], 'EdgeColor', 'none', 'FaceColor', [col(2,:) 0.1])
% rectangle('Position',[567 0 20 6], 'EdgeColor', 'none', 'FaceColor', [col(3,:) 0.1])
% 
% % plot baseline data
% for j = 1:3
%     avg = mean(movtime_mu{j}(trials{1}));
%     plot([7 dayStart(j+1)+60], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
% end
% 
% % plot data from days 1-2
% for i = 2:6
%     for j = 1:3
%         s = shadedErrorBar(trials{i},movtime_mu{gblock(j)}(trials{i}), movtime_se{gblock(j)}(trials{i}));
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 3-5
% for i = 6:15
%     for j = 1:2
%         s = shadedErrorBar(trials{i},movtime_mu{gblock(j)}(trials{i}), movtime_se{gblock(j)}(trials{i}));
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 6-10
% for i = 15:30
%     s = shadedErrorBar(trials{i},movtime_mu{3}(trials{i}), movtime_se{3}(trials{i}));
%     editErrorBar(s,col(3,:),lw);
% end
% xlabel('Day')
% xticks(dayStart)
% xticklabels(dayStartLabels)
% yticks(0:2:8)
% ylabel('Movement time (s)')
% axis([7 586 0 6])
% set(gca,'TickDir','out')
% 
% % prints figure for Illustrator
% % print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/move_time','-dpdf','-painters')
% 
% %% Figure 2G
% f = figure(15); clf; hold on 
% set(f,'Position',[200 200 140 140]);
% 
% % vertical lines delineating days
% plot([dayStart2; dayStart2],[0 2],'Color',[0.8 0.8 0.8])
% 
% % delineate flip block
% rectangle('Position',[87 0 20 2], 'EdgeColor', 'none', 'FaceColor', [col(1,:) 0.1])
% rectangle('Position',[267 0 20 2], 'EdgeColor', 'none', 'FaceColor', [col(2,:) 0.1])
% rectangle('Position',[567 0 20 2], 'EdgeColor', 'none', 'FaceColor', [col(3,:) 0.1])
% 
% % plot baseline data
% for j = 1:3
%     avg = mean(RT_mu{j}(trials{1}));
%     plot([7 dayStart(j+1)+60], [avg avg], 'Color', col(j,:), 'LineWidth', 1)
% end
% 
% % plot data from days 1-2
% for i = 2:6
%     for j = 1:3
%         s = shadedErrorBar(trials{i},RT_mu{gblock(j)}(trials{i}), RT_se{gblock(j)}(trials{i}));
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 3-5
% for i = 6:15
%     for j = 1:2
%         s = shadedErrorBar(trials{i},RT_mu{gblock(j)}(trials{i}), RT_se{gblock(j)}(trials{i}));
%         editErrorBar(s,col(gblock(j),:),lw);
%     end
% end
% 
% % plot data from days 6-10
% for i = 15:30
%     s = shadedErrorBar(trials{i},RT_mu{3}(trials{i}), RT_se{3}(trials{i}));
%     editErrorBar(s,col(3,:),lw);
% end
% xlabel('Day')
% xticks(dayStart)
% xticklabels(dayStartLabels)
% yticks(.5:.5:2)
% ylabel('Reaction time (s)')
% axis([7 586 .2 2])
% set(gca,'TickDir','out')
% 
% % prints figure for Illustrator
% % print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/rt','-dpdf','-painters')

end
