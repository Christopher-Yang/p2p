
group = 'day10';
subj = 1;
trial = 2848;

a = d.(group){subj};
idx = a.init_x(trial);

figure(1); clf; hold on
plot(a.targetAbs(trial,1),a.targetAbs(trial,2),'.k','MarkerSize',60)
plot(a.C{trial}(:,1),a.C{trial}(:,2),'k')
plot(a.C{trial}(idx,1),a.C{trial}(idx,2),'.b','MarkerSize',30)
if ~isnan(a.incorrectReach_x(trial))
    idx2 = a.init_x(trial)+20;
    if a.incorrectReach_x(trial) == 1
        plot(a.C{trial}(idx2,1),a.C{trial}(idx2,2),'.r','MarkerSize',30)
    else
        plot(a.C{trial}(idx2,1),a.C{trial}(idx2,2),'.g','MarkerSize',30)
    end
end
axis([0.4 0.8 0.05 0.45])
axis square
