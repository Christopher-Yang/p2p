
col = [50 205 50
       255 127 80
       100 149 237]./255;
subj = [4 6 5];
trials{1} = [1:30; 31:60; 201:230];
trials{2} = [1:30; 31:60; 301:330]; 
trials{3} = [1:30; 31:60; 401:430];
groups = {'day2','day5','day10'};
titles = {'Baseline','Early','Late'};
groupNames = {'2-day','5-day','10-day'};

figure(5); clf
for k = 1:3
    s = subj(k);
    for j = 1:3
        subplot(3,3,(k-1)*3+j); hold on
        for i = trials{k}(j,:)
            plot(d.(groups{k}){s}.L{i}(:,1),d.(groups{k}){s}.L{i}(:,2),'Color',col(1,:))
            plot(d.(groups{k}){s}.R{i}(:,1),d.(groups{k}){s}.R{i}(:,2),'Color',col(2,:))
            plot(d.(groups{k}){s}.C{i}(:,1),d.(groups{k}){s}.C{i}(:,2),'Color',col(3,:))
        end
        axis equal
        if k == 1
            title(titles{j})
        end
        if j == 1
            ylabel(groupNames{k})
        end
    end
end
