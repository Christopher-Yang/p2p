% clear all
% load P2P

col = [255 255 102
       255 165 0
       255 69 0]/255;
Nsubj = length(d);
Ntrials = length(d{1}.Cr);
Ntrials2 = Ntrials/length(blocks);

dir = NaN(Ntrials,Nsubj);

for j = 1:Nsubj
    dir(:,j) = d{j}.initDir-pi/2;
end

dir = unwrap(dir);
se = 2*(circ_std(dir,[],[],2)/sqrt(Nsubj))*180/pi;

dirBin2 = reshape(dir,[3 size(dir,1)/3 Nsubj]);
dirBin2 = squeeze(circ_mean(dirBin2,[],1));
dirBin = circ_mean(dirBin2,[],2)*180/pi;
seBin = 2*(circ_std(dirBin2,[],[],2)/sqrt(Nsubj))*180/pi;

dir = circ_mean(dir,[],2)*180/pi;
edges = -30:3:30;

% plot unbinned direction histograms
figure
subplot(4,1,1)
histogram(dir(1:30),edges)
axis([-25 25 0 15])
xticks(-50:5:50)
title('Baseline')

subplot(4,1,2)
histogram(dir(31:60),edges)
axis([-25 25 0 15])
xticks(-50:5:50)
title('Early Learning')
ylabel('Counts')

subplot(4,1,3)
histogram(dir(61:90),edges)
axis([-25 25 0 15])
xticks(-50:5:50)
title('Mid Learning')

subplot(4,1,4)
histogram(dir(91:120),edges)
axis([-25 25 0 15])
xticks(-50:5:50)
title('Late Learning')
xlabel('Hand Reach Angle')

% plot reach direction error
figure
x = shadedErrorBar(1:30,dir(1:30),se(1:30));
editErrorBar(x,[0 0 0],1);
hold on
plot([1 120],[0 0],'--k','LineWidth',1)
x = shadedErrorBar(31:60,dir(31:60),se(31:60));
editErrorBar(x,[0 0 0],1);
x = shadedErrorBar(61:90,dir(61:90),se(61:90));
editErrorBar(x,[0 0 0],1);
x = shadedErrorBar(91:120,dir(91:120),se(91:120));
editErrorBar(x,[0 0 0],1);
rectangle('Position',[Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(1,:) 0.1],'EdgeColor','none')
rectangle('Position',[2*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(2,:) 0.1],'EdgeColor','none')
rectangle('Position',[3*Ntrials2+1 -200 Ntrials2 400],'FaceColor',[col(3,:) 0.1],'EdgeColor','none')
axis([1 120 -30 30])
xticks(1:30:120)
ylabel(['Reach Error (',char(176),')'])
xlabel('Trials')
