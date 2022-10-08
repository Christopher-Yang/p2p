% fits mixture model to data and plots initial reach direction error

function plot_direction(data)

% set variables for analysis
rng(2);
Nblock = 8; % number of blocks
Nsubj = length(data); % number of subjects in each group

% indices for dividing up the trials into blocks
trials{1} = 1:30;
for i = 1:7
    trials{i+1} = (i-1)*100 + 31:(i-1)*100 + 130;
end
      
% set colors for generating plots
col = [180 180 0
        139 69 19
        106 90 205]./255;
    
% store all reach direction errors into dirError
Ntrials = length(data{1}.Cr);
dir = NaN(Ntrials,Nsubj);
for j = 1:Nsubj
    dir(:,j) = data{j}.initDir*180/pi;
end

dirError = dir;

%% fit mixture model for each participant

% set values for initial parameters for optimization
muInit = 0; % mean of von Mises distribution
kappaInit = 1; % concentration parameter of von Mises
weightInit = 0.9; % relative weight of von Mises and uniform distributions

% preallocate variables for parameters
mu_opt = NaN(Nblock, Nsubj);
kappa_opt = NaN(Nblock, Nsubj);
weight_opt = NaN(Nblock, Nsubj);
sd = NaN(Nblock, Nsubj);

% maximum likelihood estimation for mixture model
for m = 1:Nsubj % loop over subjects
    for j = 1:Nblock % loop over blocks
        
        % select trials to analyze and store in samples
        trial = trials{j};
        samples = dirError(trial,m)*pi/180; % convert error to radians
        samples = samples(~isnan(samples));
        
        % fit model
        log_likelihood = @(params) calc_likelihood(params, samples);
        paramsInit = [muInit kappaInit weightInit]; % set parameters to current values of mu and kappa
        params_opt = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 200 1]);
        
        % store fitted parameter values
        mu_opt(j,m) = params_opt(1);
        kappa_opt(j,m) = params_opt(2);
        weight_opt(j,m) = params_opt(3);
        
        % compute circular standard deviation
        R = besseli(1,params_opt(2))/besseli(0,params_opt(2));
        sd(j,m) = sqrt(-2 * log(R)); % circular standard deviation
    end
end

% mean and standard deviation of the circular standard deviation
sd_mu = mean(sd, 2, 'omitnan');
sd_se = std(sd, [], 2, 'omitnan')./sqrt(repmat(Nsubj,[size(sd,1) 1]));

% make vector of standard deviations for stats in R
% sd2 = permute(sd,[1 3 2]); % rearrange sd to make it easier to add to y
% y = [sd2(3,1:13,1)'; sd2(5,1:13,1)'; sd2(5,:,2)'; sd2(14,:,2)']*180/pi;

% write data for analysis in R
% groupNames(1:26,1) = "2-day";
% groupNames(27:54,1) = "5-day";
% blockNames([1:13 27:40],1) = "before";
% blockNames([14:26 41:54],1) = "after";
% subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2])]';
% T = table(groupNames, blockNames, subject, y, 'VariableNames', {'group','block','subject','sd'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/sd.csv')

%%

figure(1); clf; hold on
plot(sd * 180/pi, 'Color', [0 0 0 0.3], 'Linewidth', 0.25)
plot(mean(sd * 180/pi, 2), 'k', 'LineWidth', 3)
xlabel('Block')
ylabel(['Circular st dev (' char(176) ')'])
xticks(1:8)
yticks(0:20:80)
xticklabels({'Baseline', 1, 2, 3, 4, 5, 6, 7})

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\sd','-dpdf','-painters')

%%

gblocks = [1 2 8];

figure(2); clf; hold on
for j = 1:3
    trial = trials{gblocks(j)};
    [f,xi] = ksdensity(reshape(dirError(trial,:),[numel(dirError(trial,:)) 1]));
    plot(xi,f,'LineWidth',3,'Color',col(j,:))
    if j == 3
        axis([-180 180 0 .05])
        xlabel(['Reach direction error (' char(176) ')'])
        ylabel('Probability density')
        xticks(-180:90:180)
        box off
        set(gca,'Tickdir','out')
        yticks(0:0.025:0.05)
        if j == 1
            ylabel('Kernel-smoothed probability density')
        elseif j == 2
            xlabel('Reach direction error (degrees)')
        end
    end
end

print('C:\Users\Chris\Dropbox\Conferences\SFN 2022\ksdensity','-dpdf','-painters')

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