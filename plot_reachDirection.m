function plot_reachDirection(d)

groups = {'day2', 'day5', 'day10'};
Ngroup = length(groups);
Nsubj = [13 14 5];

% figure(1); clf
% for i = 1:3
%     group = groups{i};
%     for subj = 1:Nsubj(i)
%         targ = d.(group){subj}.targAng * 180/pi + 90;
%         reach = d.(group){subj}.initDir_noRot * 180/pi + 90;
%         
%         targ(targ > 180) = targ(targ > 180) - 360;
%         reach(reach > 180) = reach(reach > 180) - 360;
%         
%         subplot(3,4,4*(i-1)+1); hold on
%         plot([-180 180], [-180 180], 'k')
%         plot([-180 180], [0 0], 'k')
%         plot([0 0], [-180 180], 'k')
%         scatter(targ(1:30), reach(1:30), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.2)
%         if subj == 1
%             xticks(-180:180:180)
%             yticks(-180:180:180)
%             axis([-180 180 -180 180])
%             axis square
%             set(gca,'TickDir','out')
%             if i == 1
%                 title('Baseline')
%             end
%         end
%         
%         subplot(3,4,4*(i-1)+2); hold on
%         plot([-180 180], [-180 180], 'k')
%         plot([-180 180], [0 0], 'k')
%         plot([0 0], [-180 180], 'k')
%         scatter(targ(31:130), reach(31:130), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.2)
%         if subj == 1
%             xticks(-180:180:180)
%             yticks(-180:180:180)
%             axis([-180 180 -180 180])
%             axis square
%             set(gca,'TickDir','out')
%             if i == 1
%                 title('Early')
%             end
%         end
%         
%         subplot(3,4,4*(i-1)+3); hold on
%         plot([-180 180], [-180 180], 'k')
%         plot([-180 180], [0 0], 'k')
%         plot([0 0], [-180 180], 'k')
%         scatter(targ(end-199:end-100), reach(end-199:end-100), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.2)
%         if subj == 1
%             xticks(-180:180:180)
%             yticks(-180:180:180)
%             axis([-180 180 -180 180])
%             axis square
%             set(gca,'TickDir','out')
%             if i == 1
%                 title('Late')
%             end
%         end
%         
%         subplot(3,4,4*(i-1)+4); hold on
%         plot([-180 180], [-180 180], 'k')
%         plot([-180 180], [0 0], 'k')
%         plot([0 0], [-180 180], 'k')
%         scatter(targ(end-99:end), reach(end-99:end), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.2)
%         if subj == 1
%             xticks(-180:180:180)
%             yticks(-180:180:180)
%             axis([-180 180 -180 180])
%             axis square
%             set(gca,'TickDir','out')
%             if i == 1
%                 title('Flip')
%             end
%         end
%     end
% end

bins = -180:30:180;
blocks{1} = 1:30;
blocks{2} = 31:130;
blocks{3} = 131:230;
blocks{4} = 231:330;

for j = 1:Ngroup
    group = groups{j};
    for i = 1:Nsubj(j)
        targ = d.(group){i}.targAng([1:130 end-199:end]) * 180/pi + 90;
        reach = d.(group){i}.initDir_noRot([1:130 end-199:end]) * 180/pi + 90;
        
        targ(targ > 180) = targ(targ > 180) - 360;
        reach(reach > 180) = reach(reach > 180) - 360;
        
        targDir{j}(:,i) = targ;
        reachDir{j}(:,i) = reach;
    end
    
    for i = 1:length(bins)-1
        for k = 1:4
            dat = targDir{j}(blocks{k},:);
            idx = (dat >= bins(i)) == (dat < bins(i+1));
            total = sum(idx,'all');
            
            dat = reachDir{j}(blocks{k},:);
            datBin{j}(:,i,k) = histcounts(dat(idx), bins) ./ total;
        end
    end
end

col1 = [1 1 1];
col2 = [1 0 0];
Nstep = 100;
map = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2)...
    ,Nstep)', linspace(col1(3),col2(3),Nstep)'];
blockNames = {'Baseline','Early','Late','Flip'};

n = length(bins)-1;
clims = [0 1];
figure(1); clf
for j = 1:Ngroup
    for k = 1:4
        subplot(Ngroup, 4, 4*(j-1) + k); hold on
        imagesc(datBin{j}(:,:,k), clims)
        if k == 4
            plot([0 13],[0 13],'k')
            plot([0 13],[13 0],'--k')
        end
        colormap(map)
        xticks(0.5:6:12.5)
        xticklabels(-180:180:180)
        yticklabels(-180:180:180)
        yticks(0.5:6:12.5)
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

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/flip_direction','-dpdf','-painters')

end