%% plot reaction times for x- or y-axis movements

groupName = {'day2','day5','day10'};
Nsubj = [length(d.day2) length(d.day5) length(d.day10)];
Ngroup = length(groupName);

for group = 1:Ngroup
    for subj = 1:Nsubj(group)
        G = d.(groupName{group}){subj};
        
        bad.idx = find(G.incorrectReach_x==1);
        bad.idy = find(G.incorrectReach_y==1);
        good.idx = find(G.incorrectReach_x==0);
        good.idy = find(G.incorrectReach_y==0);
        
        bad.x = G.init_x(bad.idx);
        bad.y = G.init_y(bad.idy);
        good.x = G.init_x(good.idx);
        good.y = G.init_y(good.idy);
        
        count = 1;
        for i = 1:300
            if count > length(good.idx)
                goodRT.(groupName{group}).x(i:300,subj) = NaN;
                break
            elseif i == good.idx(count)
                goodRT.(groupName{group}).x(i,subj) = G.time{good.idx(count)}(good.x(count))-G.time{good.idx(count)}(G.go(good.idx(count)));
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
                goodRT.(groupName{group}).y(i,subj) = G.time{good.idy(count)}(good.y(count))-G.time{good.idy(count)}(G.go(good.idy(count)));
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
                badRT.(groupName{group}).x(i,subj) = G.time{bad.idx(count)}(bad.x(count))-G.time{bad.idx(count)}(G.go(bad.idx(count)));
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
                badRT.(groupName{group}).y(i,subj) = G.time{bad.idy(count)}(bad.y(count))-G.time{bad.idy(count)}(G.go(bad.idy(count)));
                count = count+1;
            else
                badRT.(groupName{group}).y(i,subj) = NaN;
            end
        end
    end
end

%%
for group = 1:Ngroup
    G = goodRT.(groupName{group}).x(101:200,:);
    G = G(~isnan(G));
    
    B = badRT.(groupName{group}).x(101:200,:);
    B = B(~isnan(B));
    start = 150;
    i = 1;
    while start <= max(G)-100
        x = (G>start-99) + (G<start+100);
        idx = find(x==2);
        Ngood = length(idx);
        
        x = (B>start-99) + (B<start+100);
        idx = find(x==2);
        Nbad = length(idx);
        
        if Ngood+Nbad < 10
            lateError{group}(i) = NaN;
            lateIdx{group}(i) = NaN;
        else
            lateError{group}(i) = Nbad/(Nbad+Ngood);
            lateIdx{group}(i) = 1;
        end
        total{group}(i) = Ngood+Nbad;
        start = start+1;
        i = i+1;
    end
    
    G = goodRT.(groupName{group}).x(201:300,:);
    G = G(~isnan(G));
    
    B = badRT.(groupName{group}).x(201:300,:);
    B = B(~isnan(B));
    start = 150;
    i = 1;
    while start <= max(G)-100
        x = (G>start-99) + (G<start+100);
        idx = find(x==2);
        Ngood = length(idx);
        
        x = (B>start-99) + (B<start+100);
        idx = find(x==2);
        Nbad = length(idx);
        
        if Ngood+Nbad < 10
            flipError{group}(i) = NaN;
            flipIdx{group}(i) = NaN;
        else
            flipError{group}(i) = Nbad/(Nbad+Ngood);
            flipIdx{group}(i) = 1;
        end
        
        start = start+1;
        i = i+1;
    end
end

figure(4); clf
subplot(1,2,1); hold on
for i = 1:3
    plot(150:length(lateError{i})+149,lateError{i})
end
axis([0 2200 0 1])
title('Late learning')
xlabel('Reaction Time (ms)')
ylabel('P(error|RT)')

subplot(1,2,2); hold on
for i = 1:3
    plot(150:length(flipError{i})+149,flipError{i})
end
axis([0 2200 0 1])
title('Post-flip')
xlabel('Reaction Time (ms)')
ylabel('P(error|RT)')
legend({'Day 2','Day 5','Day 10'})

figure(5); clf
for i = 1:Ngroup
    subplot(2,3,i); hold on
    plot(150:length(flipError{i})+149,flipError{i})
    axis([0 2200 0 1])
    if i == 1
        ylabel('P(incorrect|RT)')
    end
    
    subplot(2,3,i+3)
    plot(goodRT.(groupName{i}).x(201:300,:),1+0.3*(rand(100,Nsubj(i))-0.5),'k.','MarkerSize',5)
    hold on
    plot(badRT.(groupName{i}).x(201:300,:),2+0.3*(rand(100,Nsubj(i))-0.5),'k.','MarkerSize',5)
    xlim([0 2200])
    if i == 1
        yticks(1:2)
        yticklabels({'Correct Reaches','Incorrect Reaches'})
    elseif i == 2
        yticks([])
        xlabel('Reaction Time')
    else
        yticks([])
    end
end

%%
rng(1)
% plot for left hand
figure(1); clf;
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
figure(2); clf;
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

figure(3); clf
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