% clear all
% load P2P
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
        dir(:,j) = d.(groups{i}){j}.initDir*180/pi;
        for k = 1:3
            idx = (k-1)*100+1:(k-1)*100+100;
            pd = fitdist(dir(idx,j),'Normal');
            fits{i}{k}(j) = pd.sigma;
        end
    end

    se = (circ_std(dir,[],[],2)/sqrt(Nsubj));
    
%     for j = 1:3
%         pd = fitdist(reshape(dir((j-1)*100+1:(j-1)*100+100,:),[numel(dir{i}((j-1)*100+1:(j-1)*100+100,:)) 1]),'Normal');
%         fits{i}.mu(j) = pd.mu;
%         fits{i}.sigma(j) = pd.sigma;
%     end
    
    dirError{i} = dir;
    bins{i} = d.(groups{i}){1}.targBin;
end



%%

% plot unbinned direction histograms
edges = (-180:15:180)*pi/180;
figure(1)
for i = 1:Ngroup
    subplot(3,Ngroup,i)
    polarhistogram(dirError{i}(1:100,:)*pi/180,edges,'Normalization','pdf')
    title(graph_names{i})
%     rlim([0 1.75])
    
    subplot(3,Ngroup,i+3)
    polarhistogram(dirError{i}(101:200,:)*pi/180,edges,'Normalization','pdf')
%     rlim([0 1.75])
    
    subplot(3,Ngroup,i+6)
    polarhistogram(dirError{i}(201:300,:)*pi/180,edges,'Normalization','pdf')
%     rlim([0 1.75])
end

edges = -180:10:180;
figure(2); clf
for i = 1:Ngroup
    subplot(1,2,1); hold on
%     histogram(dirError{i}(1:100,:),edges,'Normalization','pdf')
    [f,xi] = ksdensity(reshape(dirError{i}(1:100,:),[numel(dirError{i}(1:100,:)) 1]));
    plot(xi,f,'LineWidth',2)
%     pd = fitdist(reshape(dirError{i}(1:100,:),[numel(dirError{i}(1:100,:)) 1]),'Normal');
%     x = -180:0.1:180;
%     y = normpdf(x,pd.mu,pd.sigma);
%     plot(x,y,'k')
    if i == 3
        title('Early Learning')
        xlabel('Reach direction error (degrees)')
        ylabel('Probability density')
        axis([-180 180 0 .04])
        xticks(-180:90:180)
        yticks(0:0.01:0.05)
        box off
        legend({'2-day','5-day','10-day'})
    end
    
    subplot(1,2,2); hold on
%     histogram(dirError{i}(101:200,:),edges,'Normalization','pdf')
    [f,xi] = ksdensity(reshape(dirError{i}(101:200,:),[numel(dirError{i}(101:200,:)) 1]));
    plot(xi,f,'LineWidth',2)
%     pd = fitdist(reshape(dirError{i}(101:200,:),[numel(dirError{i}(101:200,:)) 1]),'Normal');
%     x = -180:0.1:180;
%     y = normpdf(x,pd.mu,pd.sigma);
%     plot(x,y,'k')
    if i == 3
        title('Late Learning')
        xlabel('Reach direction error (degrees)')
        axis([-180 180 0 .04])
        xticks(-180:90:180)
        yticks(0:0.01:0.05)
        box off
    end
    
%     subplot(3,Ngroup,i+6)
% %     histogram(dirError{i}(201:300,:),edges,'Normalization','pdf')
%     ksdensity(reshape(dirError{i}(201:300,:),[numel(dirError{i}(101:200,:)) 1]))
%     axis([-180 180 0 .035])
%     xticks(-180:45:180)
%     if i == 1
%         ylabel('Post-flip')
%     elseif i == 2
%         xlabel('Reach direction error (degrees)')
%     end
%     box off
end

%%
Ntrials2 = 100;
trialsAll = {1:100,101:200,201:300};
for k = 1:3
    trials = trialsAll{k};
    for i = 1:Ngroup
        Nsubj = length(d.(groups{i}));
        for j = 1:Nsubj
            ang = atan2(d.(groups{i}){j}.targetRel(trials,2), d.(groups{i}){j}.targetRel(trials,1));
            angMir = atan2(d.(groups{i}){j}.targetRel(trials,2), -d.(groups{i}){j}.targetRel(trials,1));

            dir = d.(groups{i}){j}.initDir_noRot(trials)';

            error{i}(:,j) = dir-ang;
            errorMir{i}(:,j) = dir-angMir;
        end
        
        for l = 1:numel(error{i})
            while error{i}(l) >= pi
                error{i}(l) = error{i}(l)-2*pi;
            end
            while error{i}(l) < -pi
                error{i}(l) = error{i}(l)+2*pi;
            end
            while errorMir{i}(l) >= pi
                errorMir{i}(l) = errorMir{i}(l)-2*pi;
            end
            while errorMir{i}(l) < -pi
                errorMir{i}(l) = errorMir{i}(l)+2*pi;
            end
        end
        
        towardMir{i}.all(:,k) = sum(abs(error{i})>abs(errorMir{i}),1);
        towardMir{i}.mean(k) = mean(towardMir{i}.all(:,k),1);
    end
end

figure(3); clf; hold on
for i = 1:3
    plot(i,towardMir{1}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i,towardMir{1}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
    
    plot(i+4,towardMir{2}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i+4,towardMir{2}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
    
    plot(i+8,towardMir{3}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i+8,towardMir{3}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
end
xticks([2 6 10])
xticklabels(graph_names)
ylabel('Percent towards mirrored target')
axis([0.5 11.5 0 70])
box off

%%
edges = -180:10:180;
for j = 1:3
    figure(6+j-1); clf
    for i = 1:4
        bin = find(bins{j}==i);
        early = bin(bin<=100);
        post = bin(bin>200);
        a = (bin<=100)+(bin>200);
        a = ~a;
        late = bin(a);
        
        subplot(3,4,i)
        histogram(dirError{j}(early,:),edges,'Normalization','pdf')
        xticks(-180:90:180)
        ylim([0 .04])
        if i == 1
            title('Closest')
            ylabel('Early')
        elseif i == 2
            title('Close')
        elseif i == 3
            title('Far')
        else
            title('Farthest')
        end
        
        subplot(3,4,i+4)
        histogram(dirError{j}(late,:),edges,'Normalization','pdf')
        xticks(-180:90:180)
        ylim([0 .04])
        if i == 1
            ylabel('Late')
        end
        
        subplot(3,4,i+8)
        histogram(dirError{j}(post,:),edges,'Normalization','pdf')
        xticks(-180:90:180)
        ylim([0 .04])
        xlabel('Error (degrees)')
        if i == 1
            ylabel('Post')
        end
    end
end
