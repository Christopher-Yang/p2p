% clear all
% load P2P
groups = fieldnames(d);
names = {'Rotation','Mirror Reversal'};

col = lines;
col = col(1:7,:);
Nsubj = length(d.rot);
Ntrials = length(d.(groups{1}){1}.Cr);
Ntrials2 = Ntrials/length(blocks);

for i = 1:length(groups)
    dir = NaN(Ntrials,Nsubj);

    for j = 1:Nsubj
        dir(:,j) = d.(groups{i}){j}.initDir-pi/2;
    end

    dir = unwrap(dir);
    se = 2*(circ_std(dir,[],[],2)/sqrt(Nsubj))*180/pi;
    
    dirBin2 = reshape(dir,[3 size(dir,1)/3 Nsubj]);
    dirBin2 = squeeze(circ_mean(dirBin2,[],1));
    dirBin = circ_mean(dirBin2,[],2)*180/pi;
    seBin = 2*(circ_std(dirBin2,[],[],2)/sqrt(Nsubj))*180/pi;
    
    dir = circ_mean(dir,[],2)*180/pi;
    edges = -50:10:130;
    
    % plot unbinned direction histograms
    figure(1)
    subplot(4,2,(i-1)+1)
    histogram(dir(1:150),edges)
    axis([-60 90 0 75])
    xticks(-45:45:135)
    title(names{i})

    subplot(4,2,(i-1)+3)
    histogram(dir(151:300),edges)
    axis([-60 90 0 75])
    xticks(-45:45:135)
    title('Early Learning')
    ylabel('Counts')

    subplot(4,2,(i-1)+5)
    histogram(dir(301:450),edges)
    axis([-60 90 0 75])
    xticks(-45:45:135)
    title('Mid Learning')

    subplot(4,2,(i-1)+7)
    histogram(dir(451:600),edges)
    axis([-60 90 0 75])
    xticks(-45:45:135)
    title('Late Learning')
    xlabel('Hand Reach Angle')

    % plot reach error binned by 4 trials
    figure
%     for rot_vs_mir
    if i == 1
        x = shadedErrorBar(1:3:150,dirBin(1:50),seBin(1:50));
        editErrorBar(x,[0 0 0],1);
        hold on
        plot([1000 1001],[1000 1001],'b','LineWidth',1)
        x = shadedErrorBar(151:3:300,dirBin(51:100),seBin(51:100));
        editErrorBar(x,[0 0 0],1);
        x = shadedErrorBar(301:3:450,dirBin(101:150),seBin(101:150));
        editErrorBar(x,[0 0 0],1);
        x = shadedErrorBar(451:3:600,dirBin(151:200),seBin(151:200));
        editErrorBar(x,[0 0 0],1);
        plot([-5 1000],[0 0],'--k','LineWidth',1)
        rectangle('Position',[1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
        rectangle('Position',[Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(2,:) 0.1],'EdgeColor','none')
        rectangle('Position',[3*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(3,:) 0.1],'EdgeColor','none')
%         legend('VMR')
    else
        x = shadedErrorBar(1:3:150,dirBin(1:50),seBin(1:50));
        editErrorBar(x,[0 0 0],1);
        hold on
        plot([1000 1001],[1000 1001],'b','LineWidth',1)
        x = shadedErrorBar(151:3:300,dirBin(51:100),seBin(51:100));
        editErrorBar(x,[0 0 0],1);
        x = shadedErrorBar(301:3:450,dirBin(101:150),seBin(101:150));
        editErrorBar(x,[0 0 0],1);
        x = shadedErrorBar(451:3:600,dirBin(151:200),seBin(151:200));
        editErrorBar(x,[0 0 0],1);
        plot([-5 1000],[0 0],'--k','LineWidth',1)
        rectangle('Position',[1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
        rectangle('Position',[Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(2,:) 0.1],'EdgeColor','none')
        rectangle('Position',[3*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(3,:) 0.1],'EdgeColor','none')
%         legend('MR')
    end
    set(gca,'TickDir','out')
    axis([100 600 -60 60])
    xticks([100 151 301 451 600])
    yticks(-45:45:45)
    ylabel(['Reach Direction Error (',char(176),')'])
    xlabel('Trial Number')
end

