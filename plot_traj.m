% plots raw trajectories

function plot_traj(data)

% set variables for plotting
subj = 5; % select which participants to plot from each group
trials = [11:20; 101:110; 681:690]; % trials to plot for 2-day group
titles = {'Baseline','Early','Late'}; 
Nblock = length(titles); % number of blocks
col = [198 156 109
       0 146 69
       100 149 237
       251 176 59]./255; % colors for plotting

%% Figure 2A
figure(1); clf

for j = 1:Nblock % loop over blocks
    subplot(1,3,j); hold on
    a = data{subj};
    
    for i = trials(j,:) % loop over trials
        plot(a.targetAbs(i,1),a.targetAbs(i,2), '.', 'Color', [1 0.4 0.4], 'MarkerSize', 8); % plot targets
        plot(a.L{i}(:,1),a.L{i}(:,2),'Color',col(1,:)) % plot left hand
        plot(a.R{i}(:,1),a.R{i}(:,2),'Color',col(2,:)) % plot right hand
        plot(a.C{i}(:,1),a.C{i}(:,2),'k') % plot cursor
    end
    axis([0.1 1 -0.2 0.65])
    axis square
    title(titles{j})
    if j == 1
        plot([0.4 0.52],[0.5 0.5],'k','LineWidth',4)
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\traj','-dpdf','-painters')

end
