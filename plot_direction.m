% clear all
% load P2P
groups = fieldnames(d);
names = {'Rotation','Mirror Reversal'};

col = lines;
col = col(1:7,:);
% groups = {'rot'};
% names = {'Rotation'};
Nsubj = length(d.rot);
    
% a = figure;
% c = figure;
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
%     figure(a)
%     subplot(4,2,(i-1)+1)
%     histogram(dir(1:150),edges)
%     axis([-25 135 0 75])
%     xticks(-45:45:135)
%     title(names{i})
% 
%     subplot(4,2,(i-1)+3)
%     histogram(dir(151:300),edges)
%     axis([-25 135 0 75])
%     xticks(-45:45:135)
%     title('Early Learning')
%     ylabel('Counts')
% 
%     subplot(4,2,(i-1)+5)
%     histogram(dir(301:450),edges)
%     axis([-25 135 0 75])
%     xticks(-45:45:135)
%     title('Mid Learning')
% 
%     subplot(4,2,(i-1)+7)
%     histogram(dir(451:600),edges)
%     axis([-25 135 0 75])
%     xticks(-45:45:135)
%     title('Late Learning')
%     xlabel('Hand Reach Angle')

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

%     no_P2P_15deg
%     if i == 1
%         x = shadedErrorBar(1:30,dir(1:30),se(1:30));
%         editErrorBar(x,[1 0 0],1);
%         hold on
%         plot([1 150],[0 0],'--k','LineWidth',1)
%         x = shadedErrorBar(31:60,dir(31:60),se(31:60));
%         editErrorBar(x,[1 0 0],1);
%         x = shadedErrorBar(61:90,dir(61:90),se(61:90));
%         editErrorBar(x,[1 0 0],1);
%         x = shadedErrorBar(91:120,dir(91:120),se(91:120));
%         editErrorBar(x,[1 0 0],1);
%         x = shadedErrorBar(121:150,dir(121:150),se(121:150));
%         editErrorBar(x,[1 0 0],1);
%         rectangle('Position',[1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
%         rectangle('Position',[Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(2,:) 0.1],'EdgeColor','none')
%         rectangle('Position',[3*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(3,:) 0.1],'EdgeColor','none')
%         rectangle('Position',[4*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(4,:) 0.1],'EdgeColor','none')
%     end
%     axis([1 150 -30 30])
%     xticks(1:30:150)
%     ylabel(['Reach Error (',char(176),')'])
%     xlabel('Trials')


    % plot direction angular error/learning rate
%     figure(c)
% %     subplot(2,1,i)
%     if i == 1
%         x = shadedErrorBar(1:600,dir,se);
%         editErrorBar(x, [1 0 0],1);
%         hold on
% %         x = shadedErrorBar(151:300,dir(151:300),se(151:300));
% %         editErrorBar(x, col(2,:), 1);
% %         x = shadedErrorBar(301:450,dir(301:450),se(301:450));
% %         editErrorBar(x, [0 0 0], 1);
% %         x = shadedErrorBar(451:600,dir(451:600),se(451:600));
% %         editErrorBar(x, col(3,:), 1);
%         plot([-1 Ntrials+1],[0 0 ],'--k','LineWidth',2)
%         rectangle('Position',[1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
%         rectangle('Position',[Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(2,:) 0.1],'EdgeColor','none')
% %         rectangle('Position',[2*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
%         rectangle('Position',[3*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(3,:) 0.1],'EdgeColor','none')
%     else
%         x = shadedErrorBar(1:600,dir,se);
%         editErrorBar(x,[0 0 1],1);
%     end
%     xticks(1:150:600)
%     yticks(-45:45:135)
%     axis([1 Ntrials -90 90])
%     ylabel(['Reach Error (',char(176),')'])
%     xlabel('Trials')
%     title(names{i})
    
%     figure('Name',names{i})
%     for j = 1:600
%         switch j
%             case num2cell(1:150)
%                 subplot(2,2,1)
%                 plot(d.(groups{i}){1}.Cr{j}(:,1),d.(groups{i}){1}.Cr{j}(:,2))
%                 axis([-0.15 0.15 -0.05 0.25])
%                 hold on
%                 pbaspect([1 1 1])
%                 title('Baseline')
%             case num2cell(151:300)
%                 subplot(2,2,2)
%                 plot(d.(groups{i}){1}.Cr{j}(:,1),d.(groups{i}){1}.Cr{j}(:,2))
%                 axis([-0.15 0.15 -0.05 0.25])
%                 hold on
%                 pbaspect([1 1 1])
%                 title('Early')
%             case num2cell(301:450)
%                 subplot(2,2,3)
%                 plot(d.(groups{i}){1}.Cr{j}(:,1),d.(groups{i}){1}.Cr{j}(:,2))
%                 axis([-0.15 0.15 -0.05 0.25])
%                 hold on
%                 pbaspect([1 1 1])
%                 title('Mid')
%             case num2cell(451:600)
%                 subplot(2,2,4)
%                 plot(d.(groups{i}){1}.Cr{j}(:,1),d.(groups{i}){1}.Cr{j}(:,2))
%                 axis([-0.15 0.15 -0.05 0.25])
%                 hold on
%                 pbaspect([1 1 1])
%                 title('Late')
%         end
%     end
end

