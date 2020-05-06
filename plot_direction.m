function plot_direction(d)

groups = {'rot','mir'}; % names of groups
blocks = {'baseline','pert1','pert2','pert3'}; % names of the blocks

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
    se = (circ_std(dir,[],[],2)/sqrt(Nsubj))*180/pi;
    
    dirBin2 = reshape(dir,[3 size(dir,1)/3 Nsubj]);
    dirBin2 = squeeze(circ_mean(dirBin2,[],1));
    dirBin = circ_mean(dirBin2,[],2)*180/pi;
    seBin = (circ_std(dirBin2,[],[],2)/sqrt(Nsubj))*180/pi;
    
    dir = circ_mean(dir,[],2)*180/pi;
    edges = -50:10:130;
    
    % plot reach error binned by 3 trials
    figure(1);
    subplot(1,2,i); hold on
    x = shadedErrorBar(1:150,dir(1:150),se(1:150));
    editErrorBar(x,[0 0 0],1);
    plot([1000 1001],[1000 1001],'b','LineWidth',1)
    x = shadedErrorBar(151:300,dir(151:300),se(151:300));
    editErrorBar(x,[0 0 0],1);
    x = shadedErrorBar(301:450,dir(301:450),se(301:450));
    editErrorBar(x,[0 0 0],1);
    x = shadedErrorBar(451:600,dir(451:600),se(451:600));
    editErrorBar(x,[0 0 0],1);
    plot([-5 1000],[0 0],'--k','LineWidth',1)
    rectangle('Position',[1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
    rectangle('Position',[Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(2,:) 0.1],'EdgeColor','none')
    rectangle('Position',[3*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(3,:) 0.1],'EdgeColor','none')
    if i == 1
        title('Rotation')
        ylabel(['Reach Direction Error (',char(176),')'])
    else
        title('Mirror Reversal')
    end
    set(gca,'TickDir','out')
    axis([100 600 -90 90])
    xticks([100 151 301 451 600])
    yticks(-90:45:90)
    xlabel('Trial Number')
end

end