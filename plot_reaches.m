group = 'day10';
subj = 2;
start = 218;

figure(1)
for i = 1:30
    clf; hold on
    t = d.(group){subj}.targetAbs;
    plot(t(i+start,1),t(i+start,2),'k.','MarkerSize',75)
    plot(t(i+1+start,1),t(i+1+start,2),'r.','MarkerSize',75)
    c = d.(group){subj}.C{i+1+start};
    plot(c(:,1),c(:,2),'LineWidth',2)
    axis([0.4 0.8 0.1 0.5])
    xlabel('X position')
    ylabel('Y position')
    pbaspect([1 1 1])
    legend({'Start','Target'})
    pause
end