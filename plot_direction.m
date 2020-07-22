function plot_direction(d)

groups = {'rot','mir'}; % names of groups
Nsubj = length(d.rot); % number of subjects
Ntrials = length(d.(groups{1}){1}.Cr); % number of trials

% colors for heat map
col1 = [1 1 1]; % RGB for white
col2 = [0 0 1]; % RGB for blue
Nstep = 100;
map = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),...
    Nstep)', linspace(col1(3),col2(3),Nstep)'];

% loop for plotting figure
for i = 1:length(groups)
    dir_all = NaN(Ntrials,Nsubj);
    
    % store all reach direction error into dir_all
    for j = 1:Nsubj
        % reach direction error is centered around pi/2 so subtract off
        % pi/2
        dir_all(:,j) = d.(groups{i}){j}.initDir-pi/2; 
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
    
    edges = -180:10:180; % bins for reach direction error
    bin = 15; % size of bin for trials
    
    % calculate intensity of colors for heat map
    for j = 1:size(dir_all,1)/bin
        dir2 = dir_all(bin*(j-1)+1:bin*(j-1)+bin,:);
        dir2 = reshape(dir2,[numel(dir2) 1]);
        dirBin(:,j) = histcounts(dir2*180/pi,edges);
    end

    % plot reach error
    figure(1);
    subplot(2,1,i); hold on
    imagesc(dirBin,[0 80])
    colormap(map)
    set(gca,'TickDir','out')
    box off
    axis([0.5 150/bin*4+0.5 0.5 36.5])
    xticks(0.5:150/bin:150/bin*3+0.5)
    xticklabels({'Baseline','Early','','Late'})
    colorbar('Ticks',0:20:80)
    if i == 1
        ylabel(['Reach Direction Error (',char(176),')'])
    end
    yticks(0.5:6:36.5)
    yticklabels({'-180','-120','-60','0','60','120','180'})
end

end