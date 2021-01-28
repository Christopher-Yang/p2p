
clear tortuosity movtime RT
groups = {'day2','day5','day10'};
Ngroups = length(groups);
Nsubj = [13 14 5];
trials{1} = 1:30;
trials{2} = 31:130;
trials{3} = 131:230;
trials{4} = 231:330;
trials{5} = 331:430;

for i = 1:Ngroups
    for j = 1:Nsubj(i)
        tortuosity{i}(:,j) = d.(groups{i}){j}.pathlength./0.12;
        movtime{i}(:,j) = d.(groups{i}){j}.movtime./1000;
        RT{i}(:,j) = d.(groups{i}){j}.RT;
    end
    tortuosity_mu{i} = mean(tortuosity{i},2);
    tortuosity_sigma{i} = std(tortuosity{i},[],2);
    movtime_mu{i} = mean(movtime{i},2);
    movtime_sigma{i} = std(movtime{i},[],2);
    RT_mu{i} = mean(RT{i},2);
    RT_sigma{i} = std(RT{i},[],2);
end



col = [0 0 0
       1 0 0
       0 0 1];
   
figure(1); clf
subplot(1,3,1); hold on
avg = mean(tortuosity_mu{1}(1:30,:));
plot([151 250],[avg avg],'k','LineWidth',2)
avg = mean(tortuosity_mu{2}(1:30,:));
plot([261 360],[avg avg],'r','LineWidth',2)
avg = mean(tortuosity_mu{3}(1:30,:));
plot([371 470],[avg avg],'b','LineWidth',2)

plot(1:30,tortuosity_mu{3}(1:30,:),'b');
plot(1:30,tortuosity_mu{2}(1:30,:),'r')
plot(1:30,tortuosity_mu{1}(1:30,:),'k')

plot(41:140,tortuosity_mu{3}(31:130,:),'b')
plot(41:140,tortuosity_mu{2}(31:130,:),'r')
plot(41:140,tortuosity_mu{1}(31:130,:),'k')

plot(151:250,tortuosity_mu{3}(131:230,:),'b')
plot(151:250,tortuosity_mu{2}(131:230,:),'r')
plot(151:250,tortuosity_mu{1}(131:230,:),'k')

plot(261:360,tortuosity_mu{3}(231:330,:),'b')
plot(261:360,tortuosity_mu{2}(231:330,:),'r')

plot(371:470,tortuosity_mu{3}(331:430,:),'b')

xticks([41 151 261 371])
xticklabels([1 2 5 10])
ylim([0 8])
xlabel('Day')
ylabel('Tortuosity')
legend({'2-day','5-day','10-day'})


subplot(1,3,2); hold on
avg = mean(movtime_mu{1}(1:30,:));
plot([151 250],[avg avg],'k','LineWidth',2)
avg = mean(movtime_mu{2}(1:30,:));
plot([261 360],[avg avg],'r','LineWidth',2)
avg = mean(movtime_mu{3}(1:30,:));
plot([371 470],[avg avg],'b','LineWidth',2)

plot(1:30,movtime_mu{3}(1:30,:),'b')
plot(1:30,movtime_mu{2}(1:30,:),'r')
plot(1:30,movtime_mu{1}(1:30,:),'k')

plot(41:140,movtime_mu{3}(31:130,:),'b')
plot(41:140,movtime_mu{2}(31:130,:),'r')
plot(41:140,movtime_mu{1}(31:130,:),'k')

plot(151:250,movtime_mu{3}(131:230,:),'b')
plot(151:250,movtime_mu{2}(131:230,:),'r')
plot(151:250,movtime_mu{1}(131:230,:),'k')

plot(261:360,movtime_mu{3}(231:330,:),'b')
plot(261:360,movtime_mu{2}(231:330,:),'r')

plot(371:470,movtime_mu{3}(331:430,:),'b')

xticks([41 151 261 371])
xticklabels([1 2 5 10])
xlabel('Day')
ylim([0 10])
ylabel('Movement time (s)')
legend({'2-day','5-day','10-day'})


subplot(1,3,3); hold on
avg = mean(RT_mu{1}(1:30,:));
plot([151 250],[avg avg],'k','LineWidth',2)
avg = mean(RT_mu{2}(1:30,:));
plot([261 360],[avg avg],'r','LineWidth',2)
avg = mean(RT_mu{3}(1:30,:));
plot([371 470],[avg avg],'b','LineWidth',2)

plot(1:30,RT_mu{3}(1:30,:),'b')
plot(1:30,RT_mu{2}(1:30,:),'r')
plot(1:30,RT_mu{1}(1:30,:),'k')

plot(41:140,RT_mu{3}(31:130,:),'b')
plot(41:140,RT_mu{2}(31:130,:),'r')
plot(41:140,RT_mu{1}(31:130,:),'k')

plot(151:250,RT_mu{3}(131:230,:),'b')
plot(151:250,RT_mu{2}(131:230,:),'r')
plot(151:250,RT_mu{1}(131:230,:),'k')

plot(261:360,RT_mu{3}(231:330,:),'b')
plot(261:360,RT_mu{2}(231:330,:),'r')

plot(371:470,RT_mu{3}(331:430,:),'b')

xticks([41 151 261 371])
xticklabels([1 2 5 10])
xlabel('Day')
ylim([400 2000])
ylabel('Reaction time (ms)')
legend({'2-day','5-day','10-day'})