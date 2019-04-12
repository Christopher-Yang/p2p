
% close all
% figure
subj = 1;
groups = {'rot','mir'};
Ntrials = d.rot{1}.iDir;
for i = 1:length(groups)
    for j = 1:Ntrials
        clf;
        plot(d.(groups{i}){subj}.Cr{j}(:,1),d.(groups{i}){subj}.Cr{j}(:,2))
        hold on
        quiver(d.(groups{i}){subj}.Cr{j}(d.(groups{i}){subj}.iDir(j),1),d.(groups{i}){subj}.Cr{j}(d.(groups{i}){subj}.iDir(j),2),0.1*cos(d.(groups{i}){subj}.initDir(j)),0.1*sin(d.(groups{i}){subj}.initDir(j)))
        axis([-0.1 0.1 -0.05 0.15])
        legend(num2str(j))
        pause;
    end
end