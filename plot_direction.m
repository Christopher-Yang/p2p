% fits mixture model to data and plots initial reach direction error

function plot_direction(data)

% set variables for analysis
rng(2);
groups = fieldnames(data); % names of groups in data
graph_names = {'2-day','5-day','10-day'}; % names of groups for plotting
blocks = {'Baseline','Early','Late'}; % block names
Ngroup = length(groups); % number of groups
Nblock = length(blocks); % number of blocks
allSubj = [length(data.day2) length(data.day5) length(data.day10)]; % number of subjects in each group

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
    Ntrials = length(data.(groups{i}){1}.Cr);
    Nsubj = length(data.(groups{i}));
    dir = NaN(Ntrials,Nsubj);
    for j = 1:Nsubj
        dir(:,j) = data.(groups{i}){j}.initDir*180/pi;
    end
    
    dirError{i} = dir;
end

%% fit mixture model for each participant

% set values for initial parameters for optimization
muInit = 0; % mean of von Mises distribution
kappaInit = 1; % concentration parameter of von Mises
weightInit = 0.9; % relative weight of von Mises and uniform distributions

% blocks to analyze
gblocks2{1} = 1:5;
gblocks2{2} = 1:14;
gblocks2{3} = 1:29;

% preallocate variables for parameters
mu_opt = NaN(Nblock,Ngroup,max(allSubj));
kappa_opt = NaN(Nblock,Ngroup,max(allSubj));
weight_opt = NaN(Nblock,Ngroup,max(allSubj));
sd = NaN(length(gblocks2{3}),Ngroup,max(allSubj));

% maximum likelihood estimation for mixture model
for k = 1:Ngroup % loop over groups
    Nsubj = length(data.(groups{k})); % number of subjects in each group
    for m = 1:Nsubj % loop over subjects
        Nblock = length(gblocks2{k}); % number of blocks
        for j = 1:Nblock % loop over blocks
            
            % select trials to analyze and store in samples
            trial = trials{gblocks2{k}(j)};
            samples = dirError{k}(trial,m)*pi/180; % convert error to radians
            samples = samples(~isnan(samples));
            
            % fit model
            log_likelihood = @(params) calc_likelihood(params, samples);
            paramsInit = [muInit kappaInit weightInit]; % set parameters to current values of mu and kappa
            params_opt = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 200 1]);
            
            % store fitted parameter values
            mu_opt(j,k,m) = params_opt(1);
            kappa_opt(j,k,m) = params_opt(2);
            weight_opt(j,k,m) = params_opt(3);
            
            % compute circular standard deviation
            R = besseli(1,params_opt(2))/besseli(0,params_opt(2));
            sd(j,k,m) = sqrt(-2 * log(R)); % circular standard deviation
        end
    end
end

% mean and standard deviation of the circular standard deviation
sd_mu = nanmean(sd,3);
sd_se = nanstd(sd,[],3)./sqrt(repmat(allSubj,[size(sd,1) 1]));

% make vector of standard deviations for stats in R
sd2 = permute(sd,[1 3 2]); % rearrange sd to make it easier to add to y
y = [sd2(3,1:13,1)'; sd2(5,1:13,1)'; sd2(5,:,2)'; sd2(14,:,2)'; sd2(14,1:5,3)'; sd2(29,1:5,3)']*180/pi;

% write data for analysis in R
groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
blockNames([1:13 27:40 55:59],1) = "before";
blockNames([14:26 41:54 60:64],1) = "after";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','sd'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/sd.csv')

%% Figure 4B

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

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/ksdensity','-dpdf','-painters')

%% Figure 4C

order = [3 2 1];
Nday = [2 5 10];

f = figure(5); clf; hold on
set(f,'Position',[200 200 140 140]);
for i = 1:Ngroup
    o = order(i);
    for j = 1:Nday(o)
        if j == 1
            plot([2 3*Nday(o)-1], [sd_mu(1,o) sd_mu(1,o)]*180/pi, 'Color', col(o,:), 'LineWidth', 1)
            
            idx = 2:3;
        elseif j == Nday(o)
            idx = (j-1)*3 + (1:2);
        else
            idx = (j-1)*3 + (1:3);
        end
        
        if i == 1
            plot([idx(1)-0.5; idx(1)-0.5],[0 60],'Color',[0.8 0.8 0.8])
        end
        
        s = shadedErrorBar(idx, sd_mu(idx,o)*180/pi, sd_se(idx,o)*180/pi);
        editErrorBar(s,col(o,:),0.5);
    end
end
set(gca,'TickDir','out')
xticks([2 4 13 28])
xticklabels([1 2 5 10])
xlabel('Day')
ylabel('Circular st dev')
xlim([2 29])
yticks(0:15:60)

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/st_dev','-dpdf','-painters')

%% check fits of mixture model; can uncomment to generate plots
% subj = 1; % choose which subject to plot fits for
% Nblock = 3; % choose block for plotting
% 
% % parameters for plotting von Mises distribution
% delt = pi/32; % interval for plotting PDF
% x = -pi:delt:pi-delt; % vector of points to plot PDF
% vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
% 
% figure(1); clf
% for j = 1:Nblock
%     for k = 1:Ngroup
%         subplot(3, 3, (j-1)*3 + k); hold on
%         Nsubj = size(dirError{k},2);
%         
%         idx = gblocks{k}(j);
%         
%         % plot histograms
%         if subj <= Nsubj
%             trial = trials{idx};
%             histogram(dirError{k}(trial,subj)*pi/180,x,'Normalization','pdf');
%         end
%         
%         % plot pdf
%         pdf = vmPDF(x, mu_opt(idx,k,subj), kappa_opt(idx,k,subj)); % PDF of von Mises distribution
%         mixPDF = weight_opt(idx,k,subj) * pdf + (1-weight_opt(idx,k,subj)) * (1 / (2*pi)); % weight von Mises with uniform distribution
%         plot(x, mixPDF, 'LineWidth', 2)
%         if j == 1
%             title(graph_names{k})
%             if k == 1
%                 ylabel('Baseline')
%             end
%         elseif j == 2 && k == 1
%             ylabel('Early learning')
%         elseif j == 3 
%             if k == 1
%                 ylabel('Late learning')
%             elseif k == 2
%                 xlabel('Reach direction error (degrees)')
%             end
%         end
%         axis([-pi pi 0 5])
%     end
% end

end

% function for computing log-likelihood
function neg_log_likelihood = calc_likelihood(params,samples)
    mu = params(1);
    kappa = params(2);
    weight = params(3);
    
    likelihood_unif = ones(size(samples)) .* ((1 - weight) / (2*pi));
    likelihood_vm = weight * exp(kappa * cos(samples-mu)) / (2 * pi * besseli(0,kappa));
    
    likelihood_all = sum([likelihood_unif likelihood_vm],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end