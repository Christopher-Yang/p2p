groups = {'rot','mir'}; % names of groups
Nsubj = length(d.rot); % number of subjects
Ntrials = length(d.(groups{1}){1}.Cr); % number of trials
dir_all = NaN(Ntrials,Nsubj,2);
RT_all = NaN(Ntrials,Nsubj,2);

for i = 1:length(groups)
    % store all reach direction error into dir_all
    for j = 1:Nsubj
        dir_all(:,j,i) = d.(groups{i}){j}.initDir-pi/2;
        RT_all(:,j,i) = d.(groups{i}){j}.RT-100;
    end
end

% unwrap error to be [-pi, pi)
for j = 1:numel(dir_all)
    while abs(dir_all(j)) >= pi
        if dir_all(j) >= pi
            dir_all(j) = dir_all(j) - 2*pi;
        else
            dir_all(j) = dir_all(j) + 2*pi;
        end
    end
end

dir_all = abs(dir_all*180/pi);

%%
clear avg_baseline avg_late

times = 200:10:800;
for i = 1:2
    for j = 1:10
        for k = 1:length(times)
            a = RT_all(1:150,j,i);
            high = times(k)+50;
            low = times(k)-50;
            idx = logical((a < high) + (a > low) - 1);
            
            b = dir_all(1:150,j,i);
            avg_baseline(k,j,i) = mean(b(idx));
            
            a = RT_all(451:600,j,i);
            high = times(k)+50;
            low = times(k)-50;
            idx = logical((a < high) + (a > low) - 1);
            
            b = dir_all(451:600,j,i);
            avg_late(k,j,i) = mean(b(idx));
        end
    end
end

% avg_baseline2 = squeeze(nanmean(avg_baseline,2));
% avg_late2 = squeeze(nanmean(avg_late,2));
% 
% figure(1); clf
% for i = 1:2
%     subplot(1,2,i); hold on
% %     plot(RT_all(1:150,:,i),dir_all(1:150,:,i),'k.')
% %     plot(RT_all(451:600,:,i),dir_all(451:600,:,i),'r.')
%     plot(times,avg_baseline2(:,i),'k','LineWidth',3)
%     plot(times,avg_late2(:,i),'r','LineWidth',3)
% end

avg_baseline2 = squeeze(mean(RT_all(1:150,:,:)));
avg_late2 = squeeze(mean(RT_all(451:600,:,:)));

avg_baseline3 = squeeze(mean(dir_all(1:150,:,:)));
avg_late3 = squeeze(mean(dir_all(451:600,:,:)));

figure(1); clf
for i = 1:2
    subplot(1,2,i); hold on
    plot(avg_baseline2(:,i),avg_baseline3(:,i),'k.','MarkerSize',50)
    plot(avg_late2(:,i),avg_late3(:,i),'r.','MarkerSize',50)
    axis([300 700 0 50])
    
end

%%
for i = 1:10
    RT.rot(:,i) = [mean(d.rot{i}.RT(1:150)) mean(d.rot{i}.RT(451:600))];
    RT.mir(:,i) = [mean(d.mir{i}.RT(1:150)) mean(d.mir{i}.RT(451:600))];
end

figure(1); clf; hold on
plot(1:2, RT.rot,'r','HandleVisibility','off')
plot(1:2, mean(RT.rot,2),'r.','MarkerSize',30)
plot(3:4, RT.mir,'b','HandleVisibility','off')
plot(3:4, mean(RT.mir,2),'b.','MarkerSize',30)
xlim([0.5 4.5])
xticks(1:4)
xticklabels({'Baseline', 'Late', 'Baseline', 'Late'})
ylabel('RT (ms)')
legend({'Rotation', 'Mirror reversal'})

%%

figure(3); clf; hold on
for i = 1:2
    plot(mean(RT_all(:,:,i),2))
end