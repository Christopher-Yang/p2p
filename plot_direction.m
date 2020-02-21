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

%%

% plot unbinned direction histograms
edges = (-180:15:180)*pi/180;
figure(1)
for i = 1:Ngroup
    subplot(2,Ngroup,i)
    polarhistogram(dirError{i}(1:100,:),edges,'Normalization','pdf')
    title(graph_names{i})
%     rlim([0 1.75])
    
    subplot(2,Ngroup,i+3)
    polarhistogram(dirError{i}(101:200,:),edges,'Normalization','pdf')
%     rlim([0 1.75])
end

edges = -180:10:180;
figure(2); clf
for i = 1:Ngroup
    subplot(2,Ngroup,i)
    histogram(dirError{i}(1:100,:)*180/pi,edges)
%     ksdensity(reshape(dirError{i}(1:100,:)*180/pi,[numel(dirError{i}(1:100,:)) 1]))
%     axis([-180 180 0 250])
    xticks(-180:45:180)
    title(graph_names{i})
    if i == 1
        ylabel('Before flip')
    end
    box off
    
    subplot(2,Ngroup,i+3)
    histogram(dirError{i}(101:200,:)*180/pi,edges)
%     ksdensity(reshape(dirError{i}(101:200,:)*180/pi,[numel(dirError{i}(101:200,:)) 1]))
%     axis([-180 180 0 250])
    xticks(-180:45:180)
    if i == 1
        ylabel('After flip')
    elseif i == 2
        xlabel('Reach direction error (degrees)')
    end
    box off
end

%%
Ntrials2 = 100;
ang = NaN(Ntrials2,Ngroup);
angMir = NaN(Ntrials2,Ngroup);
trialsAll = {1:100,101:200,201:300};
for k = 1:3
    trials = trialsAll{k};
    for i = 1:Ngroup
        Nsubj = length(d.(groups{i}));
        dir = NaN(Ntrials2,Nsubj);
        for j = 1:Nsubj
            targs = [d.(groups{i}){j}.targetRel(trials,1) d.(groups{i}){j}.targetRel(trials,2)];
            targsMir = [-d.(groups{i}){j}.targetRel(trials,1) d.(groups{i}){j}.targetRel(trials,2)];
            for l = 1:Ntrials2
                ang(l) = atan2(targs(j,2),targs(j,1));
                angMir(l) = atan2(targsMir(j,2),targsMir(j,1));
            end
            dir = d.(groups{i}){j}.initDir_noRot(trials);
        
        for j = 1:Ntrials2
            ang(j,i,k) = atan2(targs(j,2),targs(j,1));
            angMir(j,i,k) = atan2(targsMir(j,2),targsMir(j,1));
        end
        for j = 1:Nsubj
            dir(:,j) = d.(groups{i}){j}.initDir_noRot(trials);
        end
        
        error{i} = dir-repmat(ang(:,i),[1 size(dir,2)]);
        errorMir{i} = dir-repmat(angMir(:,i),[1 size(dir,2)]);

        for j = 1:numel(error{i})
            while error{i}(j) >= pi
                error{i}(j) = error{i}(j)-2*pi;
            end
            while error{i}(j) < -pi
                error{i}(j) = error{i}(j)+2*pi;
            end
            while errorMir{i}(j) >= pi
                errorMir{i}(j) = errorMir{i}(j)-2*pi;
            end
            while errorMir{i}(j) < -pi
                errorMir{i}(j) = errorMir{i}(j)+2*pi;
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
    
    plot(i+3,towardMir{2}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i+3,towardMir{2}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
    
    plot(i+6,towardMir{3}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i+6,towardMir{3}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
end
xticks([1.5 4.5 7.5])
xticklabels(graph_names)
ylabel('Percent towards mirrored target')
axis([0.5 9.5 0 70])
box off

%%
[closest,close,far,farthest] = deal(zeros(Ntrials2,Ngroup,2));
for k = 1:2
    for i = 1:Ntrials2
        for j = 1:Ngroup
            a = ang(i,j,k);
            if a>pi/8 && a<3*pi/8
                closest(i,j,k) = 1;
            elseif a<-5*pi/8 && a>-7*pi/8
                closest(i,j,k) = 1;
            elseif a>5*pi/8 && a<7*pi/8
                farthest(i,j,k) = 1;
            elseif a<-pi/8 && a>-3*pi/8
                farthest(i,j,k) = 1;
            elseif a>0 && a<pi/8
                close(i,j,k) = 1;
            elseif a>3*pi/8 && a<pi/2
                close(i,j,k) = 1;
            elseif a<-pi/2 && a>-5*pi/8
                close(i,j,k) = 1;
            elseif a<-7*pi/8 && a>-pi
                close(i,j,k) = 1;
            elseif a>pi/2 && a<5*pi/8
                far(i,j,k) = 1;
            elseif a>7*pi/8 && a<pi
                far(i,j,k) = 1;
            elseif a<0 && a>-pi/8
                far(i,j,k) = 1;
            elseif a<-3*pi/8 && a >-pi/2
                far(i,j,k) = 1;
            end
        end
    end
end

figure(6); clf
for i = 1:Ngroup
    subplot(2,3,i)
    histogram(dirError{i}(find(closest(:,i,1)==1),:)*180/pi)
    subplot(2,3,i+3)
    histogram(dirError{i}(find(closest(:,i,2)==1)+100,:)*180/pi)
end

%%
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
