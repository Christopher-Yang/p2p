function plot_direction(d)

% set variables for analysis
rng(2);
groups = fieldnames(d);
graph_names = {'2-day','5-day','10-day'};
blocks = {'Baseline','Early','Late'};
Ngroup = length(groups);
Nblock = length(blocks);
allSubj = [length(d.day2) length(d.day5) length(d.day10)];

% indices for dividing up the trials into blocks
trials{1} = 1:30;
for i = 1:29
    trials{i+1} = (i-1)*100 + 31:(i-1)*100 + 130;
end

% blocks for baseline, early, and late for each group
gblocks{1} = [1 2 5];
gblocks{2} = [1 2 14];
gblocks{3} = [1 2 29];
      
% set colors for generating plots
col = [180 180 0
        0 191 255
        255 99 71]./255;
    
% store all reach direction errors into dirError
for i = 1:Ngroup
    Ntrials = length(d.(groups{i}){1}.Cr);
    Nsubj = length(d.(groups{i}));
    dir = NaN(Ntrials,Nsubj);
    for j = 1:Nsubj
        dir(:,j) = d.(groups{i}){j}.initDir*180/pi;
    end
    
    dirError{i} = dir;
    bins{i} = d.(groups{i}){1}.targBin;
end

%% fit mixture model for each participant
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution

% set values for initial parameters for optimization
muInit = 0;
kappaInit = 1;
weightInit = 0.9;

gblocks2{1} = 1:5;
gblocks2{2} = 1:14;
gblocks2{3} = 1:29;

% preallocate variables for parameters
mu_opt = NaN(Nblock,Ngroup,max(allSubj));
kappa_opt = NaN(Nblock,Ngroup,max(allSubj));
weight_opt = NaN(Nblock,Ngroup,max(allSubj));
sd = NaN(length(gblocks2{3}),Ngroup,max(allSubj));

% main loop for EM
for k = 1:Ngroup
    Nsubj = length(d.(groups{k}));
    for m = 1:Nsubj
        Nblock = length(gblocks2{k});
        for j = 1:Nblock
            
            % select trials to analyze and store in samples
            trial = trials{gblocks2{k}(j)};
            samples = dirError{k}(trial,m)*pi/180; % convert error to radians
            samples = samples(~isnan(samples));
            
            % initialize mu and kappa for the VM distribution and the
            % relative weight between the VM and uniform distributions
            mu = muInit;
            kappa = kappaInit;
            weight = weightInit;
            
            log_likelihood = @(params) calc_likelihood(params, samples);
            paramsInit = [mu kappa weight]; % set parameters to current values of mu and kappa
            params_opt = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 200 1]);
            mu = params_opt(1);
            kappa = params_opt(2);
            weight = params_opt(3);
            
            % store fitted parameter values
            mu_opt(j,k,m) = mu;
            kappa_opt(j,k,m) = kappa;
            weight_opt(j,k,m) = weight;
            
            % compute circular standard deviation
            R = besseli(1,kappa)/besseli(0,kappa);
            sd(j,k,m) = sqrt(-2 * log(R)); % circular standard deviation
        end
    end
end

vmFit.mu_opt = mu_opt;
vmFit.kappa_opt = kappa_opt;
vmFit.weight_opt = weight_opt;
vmFit.sd = sd;

% save('variables/vmFit.mat','vmFit')

% points to assess the PDF
delt = pi/64;
x = -pi:delt:pi-delt;

% mean and standard deviation of the circular standard deviation
sd_mu = nanmean(sd,3);
sd_se = nanstd(sd,[],3)./sqrt(repmat(allSubj,[size(sd,1) 1]));

sd2 = permute(sd,[1 3 2]);

y = [sd2(3,1:13,1)'; sd2(5,1:13,1)'; sd2(5,:,2)'; sd2(14,:,2)'; sd2(14,1:5,3)'; sd2(29,1:5,3)']*180/pi;

groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
blockNames([1:13 27:40 55:59],1) = "before";
blockNames([14:26 41:54 60:64],1) = "after";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','sd'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/sd.csv')

%% plot mixture model fit on top of data histograms
subj = 1; % choose which subject to plot fits for
Nblock = 3;

figure(1); clf
for j = 1:Nblock
    for k = 1:Ngroup
        subplot(3, 3, (j-1)*3 + k); hold on
        Nsubj = size(dirError{k},2);
        
        idx = gblocks{k}(j);
        
        % plot histograms
        if subj <= Nsubj
            trial = trials{idx};
            histogram(dirError{k}(trial,subj),x*180/pi,'Normalization','probability');
        end
        
        % plot pdf
        pdf = vmPDF(x, mu_opt(idx,k,subj), kappa_opt(idx,k,subj)); % PDF of von Mises distribution
        mixPDF = weight_opt(idx,k,subj) * pdf + (1-weight_opt(idx,k,subj)) * (1 / (2*pi)); % weight von Mises with uniform distribution
        plot(x*180/pi, mixPDF./sum(mixPDF), 'LineWidth', 2)
        if j == 1
            title(graph_names{k})
        elseif j == 2 && k == 1
            ylabel('Probability')
        elseif j == 3 && k == 2
            xlabel('Reach Direction Error (degrees)')
        end
        axis([-180 180 0 0.3])
    end
end

%%
subj = 1;
block = 28;
group = 3;
trial = trials{block};
delt = pi/16;
x = -pi:delt:pi-delt;

f = figure(2); clf
set(f,'Position',[200 200 170 120]);
polarhistogram(dirError{group}(trial,subj)*pi/180,x,'FaceColor','black','FaceAlpha',.3,'Normalization','pdf'); hold on

delt = pi/64;
x = -pi:delt:pi-delt;
pdf = vmPDF(x, mu_opt(block,group,subj), kappa_opt(block,group,subj)); % PDF of von Mises distribution
mixPDF = weight_opt(block,group,subj) * pdf + (1-weight_opt(block,group,subj)) * (1 / (2*pi)); % weight von Mises with uniform distribution
polarplot(x,mixPDF,'Color',[1 0 0 0.7],'LineWidth',1)
thetaticks(0:90:270)
rticks(1:2)

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/polar','-dpdf','-painters')

%% plot standard deviation of von Mises distribution (x-axis unscaled by practice time)

order = [3 2 1];
Nday = [2 5 10];

f = figure(3); clf; hold on
set(f,'Position',[200 200 140 140]);
for i = 1:Ngroup
    o = order(i);
    for j = 1:Nday(o)
        if j == 1
%             s = shadedErrorBar([-1 1], [sd_mu(1,o) sd_mu(1,o)]*180/pi, [sd_se(1,o) sd_se(1,o)]*180/pi);
%             editErrorBar(s,col(o,:),0.5);
            plot([2 3*Nday(o)-1], [sd_mu(1,o) sd_mu(1,o)]*180/pi, 'Color', col(o,:), 'LineWidth', 1)
            
            idx = 2:3;
        elseif j == Nday(o)
            idx = (j-1)*3 + (1:2);
        else
            idx = (j-1)*3 + (1:3);
        end
        
        s = shadedErrorBar(idx, sd_mu(idx,o)*180/pi, sd_se(idx,o)*180/pi);
        editErrorBar(s,col(o,:),0.5);
        
%         if o == 3 && j > 1
%             plot([idx(1)-0.5 idx(1)-0.5], [0 60], 'k')
%         end
    end
end
set(gca,'TickDir','out')
xticks([2 4 13 28])
xticklabels([1 2 5 10])
xlabel('Day')
ylabel('Circular st dev')
xlim([2 29])
yticks(0:15:60)

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/st_dev','-dpdf','-painters')

%% plot kernel-smoothed PDF
f = figure(4); clf
set(f,'Position',[200 200 500 140]);
for i = 1:Ngroup
    for j = 1:3
        trial = trials{gblocks{4-i}(j)};
        subplot(1,3,j); hold on
        [f,xi] = ksdensity(reshape(dirError{4-i}(trial,:),[numel(dirError{4-i}(trial,:)) 1]));
        plot(xi,f,'LineWidth',1,'Color',col(4-i,:))
        if i == 3
            title(blocks{j})
            axis([-90 90 0 .06])
            xticks(-180:90:180)
            box off
            set(gca,'Tickdir','out')
            yticks(0:0.03:0.06)
            if j == 1
                ylabel('Kernel-smoothed probability density')
            elseif j == 2
                xlabel('Reach direction error (degrees)')
            end
        end
    end
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/ksdensity','-dpdf','-painters')

%% correlate proportion of away trials with standard deviation of von Mises distribution

% set variables for analysis
gblocks2 = [6 15 30]; % index of flip blocks for each group
habit = NaN(max(allSubj), Ngroup);
idx = [5 45]; % for plotting best-fit lines

% store proportion of away trials in "habit"
for i = 1:Ngroup
    Nsubj = length(d.(groups{i}));
    trialIdx = trials{gblocks2(i)};
    
    for j = 1:Nsubj
        a = d.(groups{i}){j}.incorrectReach_x(trialIdx);
        num = nansum(a);
        den = 100-sum(isnan(a));
        habit(j,i) = 100*num/den;
    end
end

figure(5); clf; hold on
for i = 1:Ngroup % plot raw data
    plot(sd2(gblocks{i}(end),:,i)*180/pi, habit(:,i), '.', 'Color', col(i,:), 'MarkerSize', 30)
end

% plot best-fit lines from least-squares regression
p = polyfit(sd2(gblocks{1}(end),1:13,1)'*180/pi, habit(1:13,1), 1);
plot(idx, p(1)*idx + p(2), 'Color', col(1,:))
p = polyfit(sd2(gblocks{2}(end),:,2)'*180/pi, habit(:,2), 1);
plot(idx, p(1)*idx + p(2), 'Color', col(2,:))
p = polyfit(sd2(gblocks{3}(end),1:5,3)'*180/pi, habit(1:5,3), 1);
plot(idx, p(1)*idx + p(2), 'Color', col(3,:))

xlabel('St dev of von Mises')
ylabel('Proportion of away trials (%)')
legend(graph_names, 'Location', 'southeast')
ylim([10 60])

%% plot reach direction histograms binned by target direction
edges = -180:10:180;
for j = 1:Ngroup
    figure(6+j-1); clf
    for i = 1:4
        bin = bins{j}==i;
        binTrials = find(bins{j} == i);
        early = binTrials(logical((binTrials>=trials{2}(1)) + (binTrials<=trials{2}(end))-1));
        late = binTrials(logical((binTrials>=trials{gblocks{j}(end)}(1)) + (binTrials<=trials{gblocks{j}(end)}(end))-1));
        post = binTrials(logical((binTrials>=trials{gblocks{j}(end)+1}(1)) + (binTrials<=trials{gblocks{j}(end)+1}(end))-1));
        
        subplot(3,4,i)
        histogram(dirError{j}(early,:),edges,'Normalization','pdf')
        xticks(-180:90:180)
        ylim([0 .04])
        if i == 1
            title('Closest')
            ylabel('Early')
        elseif i == 2
            title('Close')
        elseif i == 3
            title('Far')
        else
            title('Farthest')
        end
        
        subplot(3,4,i+4)
        histogram(dirError{j}(late,:),edges,'Normalization','pdf')
        xticks(-180:90:180)
        ylim([0 .04])
        if i == 1
            ylabel('Late')
        end
        
        subplot(3,4,i+8)
        histogram(dirError{j}(post,:),edges,'Normalization','pdf')
        xticks(-180:90:180)
        ylim([0 .04])
        xlabel('Error (degrees)')
        if i == 1
            ylabel('Post')
        end
    end
end

end

% function for computing log-likelihod
function neg_log_likelihood = calc_likelihood(params,samples)
    mu = params(1);
    kappa = params(2);
    weight = params(3);
    
    likelihood_unif = ones(size(samples)) .* ((1 - weight) / (2*pi));
    likelihood_vm = weight * exp(kappa * cos(samples-mu)) / (2 * pi * besseli(0,kappa));
    
    likelihood_all = sum([likelihood_unif likelihood_vm],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end