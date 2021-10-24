function plot_heatmap(data)

% set variables for plotting
groups = {'day2', 'day5', 'day10'}; % names of groups
blockNames = {'Baseline','Early','Late','Flip'}; % names of blocks
Ngroup = length(groups); % number of groups
Nsubj = [13 14 5]; % number of subjects in each group
bins = -180:30:180; % bins to generate heat maps
blocks{1} = 1:30; % block indices
blocks{2} = 31:130;
blocks{3} = 131:230;
blocks{4} = 231:330;

% generate color map for heatmap
col1 = [1 1 1];
col2 = [1 0 0];
Nstep = 100;
map = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2)...
    ,Nstep)', linspace(col1(3),col2(3),Nstep)'];

% generate data for heatmap
for j = 1:Ngroup % loop over groups
    
    group = groups{j};
    
    for i = 1:Nsubj(j) % loop over subjects
        
        % get target direction and initial reach direction, offsetting by
        % 90 deg so that positive y-axis is defined as 0
        targ = data.(group){i}.targAng([1:130 end-199:end]) * 180/pi + 90;
        reach = data.(group){i}.initDir_noRot([1:130 end-199:end]) * 180/pi + 90;
        
        % unwrap targ and reach to be between [-180, 180) deg
        targ(targ > 180) = targ(targ > 180) - 360;
        reach(reach > 180) = reach(reach > 180) - 360;
        
        % store data
        targDir{j}(:,i) = targ;
        reachDir{j}(:,i) = reach;
    end
    
    % generate histogram counts of data
    for i = 1:length(bins)-1 % loop over bins
        for k = 1:length(blocks) % loop over blocks

            % find indices of trials where the target direction was within
            % a given bin
            dat = targDir{j}(blocks{k},:);
            idx = (dat >= bins(i)) == (dat < bins(i+1));
            total = sum(idx,'all'); % sum the number of trials in bin
            
            % generate histogram of initial reach direction within a bin
            dat = reachDir{j}(blocks{k},:);
            datBin{j}(:,i,k) = histcounts(dat(idx), bins) ./ total;
        end
    end
end

%% Figure 4B
n = length(bins)-1;
clims = [0 1];
figure(1); clf
for j = 1:Ngroup
    for k = 1:4
        
        % plot data
        subplot(Ngroup, 4, 4*(j-1) + k); hold on
        imagesc(datBin{j}(:,:,k), clims)
        
        % plot diagonal lines 
        plot([0 13],[0 13],'k')
        plot([0 13],[13 0],'--k')
        
        colormap(map)
        xticks(0.5:6:12.5)
        yticks(0.5:6:12.5)
        xticklabels(-180:180:180)
        yticklabels(-180:180:180)
        axis([0.5 n+0.5 0.5 n+0.5])
        axis square
        set(gca,'TickDir','out')
        
        if j == 1
            title(blockNames{k})
        elseif j == 2
            if k == 1
                ylabel('Initial reach direction')
            end
        else
            if k == 2
                xlabel('Target direction')
                c = colorbar;
                c.Ticks = [0 1];
            end
        end
    end
end

% save plot for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/flip_direction','-dpdf','-painters')

%% Figure S2A
f = figure(2); clf
set(f,'Position',[200 200 300 150]);
for i = 1:2
    
    subplot(1,2,i); hold on
    
    % plot data
    imagesc(datBin{1}(:,:,4), clims)
    colormap(map)
    
    % draw von Mises model
    if i == 1
        
        plot([0 13],[0 13],'k')
        plot([0 13],[13 0],'--k')
    
    % draw uniform model
    else
        
        % draw horizontal and vertical lines
        plot([3.5 3.5],[0.5 n+0.5],'k')
        plot([9.5 9.5],[0.5 n+0.5],'k')
        plot([0.5 n+0.5],[3.5 3.5],'k')
        plot([0.5 n+0.5],[9.5 9.5],'k')
        
        % draw diagonal lines
        plot([0.5 6.5],[6.5 12.5],'k')
        plot([0.5 1.5],[8.5 9.5],'k')
        plot([0.5 3.5],[4.5 7.5],'k')
        plot([1.5 3.5],[3.5 5.5],'k')
        plot([3.5 4.5],[11.5 12.5],'k')
        plot([5.5 8.5],[9.5 12.5],'k')
        plot([7.5 9.5],[9.5 11.5],'k')
        plot([3.5 5.5],[1.5 3.5],'k')
        plot([4.5 7.5],[0.5 3.5],'k')
        plot([6.5 12.5],[0.5 6.5],'k')
        plot([8.5 9.5],[0.5 1.5],'k')
        plot([11.5 12.5],[3.5 4.5],'k')
        plot([9.5 12.5],[5.5 8.5],'k')
        plot([9.5 11.5],[7.5 9.5],'k')
    end
    axis([0.5 n+0.5 0.5 n+0.5])
    xticks(0.5:6:12.5)
    yticks(0.5:6:12.5)
    xticklabels(-180:180:180)
    yticklabels(-180:180:180)
    axis square
end

% save plot for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/model_example','-dpdf','-painters')