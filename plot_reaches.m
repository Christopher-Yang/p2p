group = 'day10';
subj = 3;
start = 201;

figure(1)
for i = 1:30
    clf; hold on
    t = d.(group){subj}.targetAbs;
    plot(t(i+start,1),t(i+start,2),'k.','MarkerSize',75) % starting target
    plot(t(i+1+start,1),t(i+1+start,2),'r.','MarkerSize',75) % ending target
    plot(t(i+start,1)-(t(i+1+start,1)-t(i+start,1)),t(i+1+start,2),'ro','MarkerSize',75) %
    plot([t(i+start,1) t(i+start,1)],[0.1 0.5],'--k')
    c = d.(group){subj}.C{i+1+start};
    x = cos(d.(group){subj}.initDir_noRot(i+1+start))*0.1;
    y = sin(d.(group){subj}.initDir_noRot(i+1+start))*0.1;
    time = d.(group){subj}.iDir(i+1+start);
    plot(c(:,1),c(:,2),'LineWidth',2)
%     quiver(c(time,1),c(time,2),x,y,0.5,'LineWidth',2)
    plot([0.5 0.6],[0.4 0.4],'k','LineWidth',3)
    axis([0.4 0.8 0.1 0.5])
    xlabel('X position')
    ylabel('Y position')
    pbaspect([1 1 1])
    legend({'Start','Target'})
    pause
end