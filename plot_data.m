% plot aggregate data across subjects
% plot target jump responses (as position or velocity) as well as basic
% trajectory metrics (path length etc)
%
clear all
load P2P
groups = {'rot','mir'};
Nsubj = 10;
g = NaN(10,600);
col = lines;
col = col(1:7,:);

figure
for i = 1:length(groups)
    pathlength = g;
    RT = g;
    pkVel = g;
    for j = 1:Nsubj
        pathlength(j,:) = d.(groups{i}){j}.pathlength;
        RT(j,:) = d.(groups{i}){j}.RT;
        pkVel(j,:) = d.(groups{i}){j}.pkVel;
    end
    
    pathlength_se = 2*std(pathlength)'/sqrt(Nsubj);
    RT_se = 2*std(RT)'/sqrt(Nsubj);
    pkVel_se = 2*std(pkVel)'/sqrt(Nsubj);
    
    pathlength = mean(pathlength)';
    RT = mean(RT)';
    pkVel = mean(pkVel)';
    
    subplot(1,3,1)
    s = shadedErrorBar(1:600,pathlength,pathlength_se);
    editErrorBar(s,col(i,:),1);
    xlabel('Trials')
    ylabel('Path length')
    ylim([0 1.5])
    xticks(0:150:600)
    
    subplot(1,3,2)
    s = shadedErrorBar(1:600,RT,RT_se);
    editErrorBar(s,col(i,:),1);
    xlabel('Trials')
    ylabel('Reaction Time')
    ylim([0 2000])
    xticks(0:150:600)
    
    subplot(1,3,3)
    if i == 1
        plot([0 1],[1 1],'Color',col(1,:))
        hold on
        plot([0 1],[1 1],'Color',col(2,:))
    end
    s = shadedErrorBar(1:600,pkVel,pkVel_se);
    editErrorBar(s,col(i,:),1);
    xlabel('Trials')
    ylabel('Peak Velocity')
    ylim([0 0.005])
    xticks(0:150:600)
end
legend({'Rotation','Mirror'},'Location','southeast')