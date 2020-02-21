%% plot reaction times for x- or y-axis movements

groupName = {'day2','day5','day10'};
Nsubj = [12 8 4];
Ngroup = length(groupName);

for group = 1:Ngroup
    for subj = 1:Nsubj(group)
        a = d.(groupName{group}){subj};
        
        bad.idx = find(a.incorrectReach_x==1);
        bad.idy = find(a.incorrectReach_y==1);
        good.idx = find(a.incorrectReach_x==0);
        good.idy = find(a.incorrectReach_y==0);
        
        bad.x = a.init_x(bad.idx);
        bad.y = a.init_y(bad.idy);
        good.x = a.init_x(good.idx);
        good.y = a.init_y(good.idy);
        
        count = 1;
        for i = 1:300
            if count > length(good.idx)
                goodRT.(groupName{group}).x(i:300,subj) = NaN;
                break
            elseif i == good.idx(count)
                goodRT.(groupName{group}).x(i,subj) = a.time{good.idx(count)}(good.x(count))-a.time{good.idx(count)}(a.go(good.idx(count)));
                count = count+1;
            else
                goodRT.(groupName{group}).x(i,subj) = NaN;
            end
        end
        
        count = 1;
        for i = 1:300
            if count > length(good.idy)
                goodRT.(groupName{group}).y(i:300,subj) = NaN;
                break
            elseif i == good.idy(count)
                goodRT.(groupName{group}).y(i,subj) = a.time{good.idy(count)}(good.y(count))-a.time{good.idy(count)}(a.go(good.idy(count)));
                count = count+1;
            else
                goodRT.(groupName{group}).y(i,subj) = NaN;
            end
        end
        
        count = 1;
        for i = 1:300
            if count > length(bad.idx)
                badRT.(groupName{group}).x(i:300,subj) = NaN;
                break
            elseif i == bad.idx(count)
                badRT.(groupName{group}).x(i,subj) = a.time{bad.idx(count)}(bad.x(count))-a.time{bad.idx(count)}(a.go(bad.idx(count)));
                count = count+1;
            else
                badRT.(groupName{group}).x(i,subj) = NaN;
            end
        end
        
        count = 1;
        for i = 1:300
            if count > length(bad.idy)
                badRT.(groupName{group}).y(i:300,subj) = NaN;
                break
            elseif i == bad.idy(count)
                badRT.(groupName{group}).y(i,subj) = a.time{bad.idy(count)}(bad.y(count))-a.time{bad.idy(count)}(a.go(bad.idy(count)));
                count = count+1;
            else
                badRT.(groupName{group}).y(i,subj) = NaN;
            end
        end
    end
end
rng(1)

% plot for left hand
figure(6); clf;
for i = 1:Ngroup
    subplot(2,3,i)
    plot(1+0.3*(rand(100,Nsubj(i))-0.5),goodRT.(groupName{i}).x(101:200,:)/1000,'k.','MarkerSize',5)
    hold on
    plot(2+0.3*(rand(100,Nsubj(i))-0.5),badRT.(groupName{i}).x(101:200,:)/1000,'k.','MarkerSize',5)
    xticks([])
    ylim([0 4])
    title('Late training')
    if i == 1
        ylabel('Reaction time (s)')
    end
    
    subplot(2,3,i+3)
    plot(1+0.3*(rand(100,Nsubj(i))-0.5),goodRT.(groupName{i}).x(201:300,:)/1000,'k.','MarkerSize',5)
    hold on
    plot(2+0.3*(rand(100,Nsubj(i))-0.5),badRT.(groupName{i}).x(201:300,:)/1000,'k.','MarkerSize',5)
    ylim([0 4])
    title('Post-flip')
    xticks(1:2)
    xticklabels({'Correct reaches','Incorrect reaches'})
    if i == 1
        ylabel('Reaction time (s)')
    end
end

% plot for right hand
figure(7); clf;
for i = 1:Ngroup
    subplot(2,3,i)
    plot(1+0.3*(rand(100,Nsubj(i))-0.5),goodRT.(groupName{i}).y(101:200,:)/1000,'k.','MarkerSize',5)
    hold on
    plot(2+0.3*(rand(100,Nsubj(i))-0.5),badRT.(groupName{i}).y(101:200,:)/1000,'k.','MarkerSize',5)
    xticks([])
    ylim([0 4])
    title('Late training')
    if i == 1
        ylabel('Reaction time (s)')
    end
    
    subplot(2,3,i+3)
    plot(1+0.3*(rand(100,Nsubj(i))-0.5),goodRT.(groupName{i}).y(201:300,:)/1000,'k.','MarkerSize',5)
    hold on
    plot(2+0.3*(rand(100,Nsubj(i))-0.5),badRT.(groupName{i}).y(201:300,:)/1000,'k.','MarkerSize',5)
    ylim([0 4])
    title('Post-flip')
    xticks(1:2)
    xticklabels({'Correct reaches','Incorrect reaches'})
    if i == 1
        ylabel('Reaction time (s)')
    end
end

%% plot reaction times of whole trajectory
groups = fieldnames(d);
graph_names = {'Day 2','Day 5','Day 10'};
names = {'day2','day5','day10'};
Ngroup = length(groups);

col = lines;
col = col(1:7,:);
Ntrials = length(d.(groups{1}){1}.Cr);
Nblock = 2;
Nsubj2 = [length(d.day2) length(d.day5) length(d.day10)];

for i = 1:Ngroup
    Nsubj = length(d.(groups{i}));
    dir = NaN(Ntrials,Nsubj);

    for j = 1:Nsubj
        dir(:,j) = d.(groups{i}){j}.initDir;
    end

    se = (circ_std(dir,[],[],2)/sqrt(Nsubj))*180/pi;
    
    dirBin2 = reshape(dir,[4 size(dir,1)/4 Nsubj]);
    dirBin2 = squeeze(circ_mean(dirBin2,[],1));
    dirBin = circ_mean(dirBin2,[],2)*180/pi;
    seBin = (circ_std(dirBin2,[],[],2)/sqrt(Nsubj))*180/pi;
    
%     dirError(:,i) = circ_mean(dir,[],2)*180/pi;
    dirError{i} = dir;
end

for j = 1:length(names)
    for i = 1:Nsubj2(j)
        early{j}(:,i) = d.(names{j}){i}.RT(1:100);
        before{j}(:,i) = d.(names{j}){i}.RT(101:200);
        after{j}(:,i) = d.(names{j}){i}.RT(201:300);
    end
end

figure(4); clf
for j = 1:3
    subplot(1,3,j); hold on
    plot(abs(dirError{j}(1:100,:)*180/pi), early{j},'.','Color',col(2,:),'MarkerSize',15)
    plot(abs(dirError{j}(201:300,:)*180/pi), after{j},'.','Color',col(1,:),'MarkerSize',15)
    plot(abs(dirError{j}(101:200,:)*180/pi), before{j},'k.','MarkerSize',15)
    if j == 1
        title('2 days of training')
        ylabel('Reaction time (ms)')
    elseif j == 2
        title('5 days of training')
        xlabel('|Reach direction error| (degrees)')
    elseif j == 3
        title('10 days of training')
    end
    axis([0 180 0 3000]) % note there are two outliers for reaction time in 10-day group
    pbaspect([1 1 1])
end