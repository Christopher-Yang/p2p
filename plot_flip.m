% fit mixture model to data and plot habitual behavior from flip block

function plot_flip(data)

% set variables for analysis
rng(2);
groups = fieldnames(data); % names of groups in data
Ngroup = length(groups); % number of groups
allSubj = [length(data.day2) length(data.day5) length(data.day10)]; % number of subjects in each group

% indices for dividing up the trials into blocks
trials{1} = 1:30;
for i = 1:29
    trials{i+1} = (i-1)*100 + 31:(i-1)*100 + 130;
end

% blocks for baseline, early, and late for each group
gblocks{1} = [1 2 5 6];
gblocks{2} = [1 2 14 15];
gblocks{3} = [1 2 29 30];
Nblock = size(gblocks{1},2);

% set colors for generating plots
col = [180 180 0
        0 191 255
        255 99 71]./255;

% store all reach direction errors into dirError
for i = 1:Ngroup

    Ntrials = length(data.(groups{i}){1}.Cr);
    Nsubj = length(data.(groups{i}));

    % preallocate variables
    dir_rot{i} = NaN(100, allSubj(i), 4);
    target_gd{i} = NaN(100, 4);
    target_hab{i} = NaN(100, 4);
    dir = NaN(Ntrials,Nsubj);
    
    % store data from all subjects in dir
    for j = 1:Nsubj
        dir(:,j) = data.(groups{i}){j}.initDir_noRot;
    end
    
    for k = 1:Nblock

        % select data from single blocks
        dir2 = dir(trials{gblocks{i}(k)},:);
        
        % relative direction of each target
        targ = data.(groups{i}){1}.targetRel(trials{gblocks{i}(k)},:);
        
        if k == 1
            idx = 1:30;
        else
            idx = 1:100;
        end
        
        % store data
        target_gd{i}(idx,k) = atan2(targ(:,2), targ(:,1));
        target_hab{i}(idx,k) = atan2(targ(:,2), -targ(:,1));
        dir_rot{i}(idx,:,k) = dir2;
    end
end

%% perform MLE for von Mises model

% probability density function of von Mises distribution
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa)));

% initialize parameters for optimization
weight1_init = 0.33;
weight2_init = 0.33;
kappa_init = 15;

% inequalities for optimization
A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
b = [0 0 1 200 0]';

% optimization loop
for k = 1:Nblock % loop over blocks
    for j = 1:Ngroup % loop over groups
        for i = 1:allSubj(j) % loop over subjects
            
            % indices of data which are not NaNs (rows 31-100 for baseline)
            idx = ~isnan(dir_rot{j}(:,i,k));
            
            % optimization
            log_likelihood = @(params) calc_likelihood(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k), target_hab{j}(idx,k));
            paramsInit = [weight1_init weight2_init kappa_init];
            [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
            
            % placeholder variables for fitted parameters
            weight1 = params_opt(1);
            weight2 = params_opt(2);
            kappa = params_opt(3);
            log_likelihood = -neg_log_likelihood;
            
            % store data
            weight1_opt{j}(i,k) = weight1;
            weight2_opt{j}(i,k) = weight2;
            weight3_opt{j}(i,k) = 1 - weight1 - weight2;
            mix_BIC{j}(i,k) = length(paramsInit) * log(sum(idx)) - 2 * log_likelihood; % calculate BIC
            
            % classify each trial as goal-directed, habitual, or random
            if k == 4 % for flip block
                dat = dir_rot{j}(:,i,k);
                
                % PDF of fitted mixture model
                mix1 = weight1 * vmPDF(dat, target_gd{j}(:,k), kappa);
                mix2 = weight2 * vmPDF(dat, target_hab{j}(:,k), kappa);
                mix3 = (1 - weight1 - weight2) * repelem(1/(2*pi), length(dat))';
                total = mix1 + mix2 + mix3;
                
                % compute responsibilities of each component
                Pr_gd = mix1 ./ total; % responsibility of goal-directed
                Pr_hab = mix2 ./ total; % responsibility of habitual
                Pr_rand = mix3 ./ total; % responsibility of random
                
                % identify which component has most responsibility
                % Pr_idx == 1: goal-directed
                % Pr_idx == 2: habitual
                % Pr_idx == 3: random
                Pr_all = [Pr_gd Pr_hab Pr_rand];
                [~, idx] = max(Pr_all, [], 2);
                idx(isnan(dat)) = NaN;                
                Pr_idx{j}(:,i) = idx;
            end
        end
    end
end

% save weights for analysis with tracking code
% save weight2_opt weight2_opt

%% proportion of habitual reaches based on target direction

for i = 1:Ngroup
    Nsubj = length(data.(groups{i}));
    idx = data.(groups{i}){1}.targBin(end-99:end); % target bins
    
    for j = 1:Nsubj
        dat = Pr_idx{i}(:,j); % reach type
        
        for k = 1:4
            datBin = dat(idx == k);
            Nreaches = length(datBin);
            
            gdBin{i}(k,j) = sum(datBin == 1) / Nreaches;
            habBin{i}(k,j) = sum(datBin == 2) / Nreaches;
        end
    end
end

figure(20); clf
subplot(1,2,1); hold on
for k = 1:4
    for i = 1:Ngroup
        jitter = (rand(length(data.(groups{i})), 1) - 0.5) * 0.5;
        plot((i-1) * 5 + k + jitter, gdBin{i}(k,:), '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility', 'off')
        plot((i-1) * 5 + k, mean(gdBin{i}(k,:)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'LineWidth', 1)
    end
end
yticks(0:0.25:1)
ylabel('Proportion of goal-directed reaches')
xticks([1:4 6:9 11:14])
xticklabels({'Closest', 'Close', 'Far', 'Farthest', 'Closest', 'Close', 'Far', 'Farthest', 'Closest', 'Close', 'Far', 'Farthest'})
xtickangle(45)

subplot(1,2,2); hold on
for k = 1:4
    for i = 1:Ngroup
        jitter = (rand(length(data.(groups{i})), 1) - 0.5) * 0.5;
        plot((i-1) * 5 + k + jitter, habBin{i}(k,:), '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility', 'off')
        plot((i-1) * 5 + k, mean(habBin{i}(k,:)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'LineWidth', 1)
    end
end
yticks(0:0.25:1)
ylabel('Proportion of habitual reaches')
xticks([1:4 6:9 11:14])
xticklabels({'Closest', 'Close', 'Far', 'Farthest', 'Closest', 'Close', 'Far', 'Farthest', 'Closest', 'Close', 'Far', 'Farthest'})
xtickangle(45)
legend({'2-day','5-day','10-day'}, 'Location', 'northwest')

%% perform MLE for uniform model
for k = 1:Nblock % loop over blocks
    for j = 1:Ngroup % loop over groups
        for i = 1:allSubj(j) % loop over subjects
            
            % indices of data which are not NaNs (rows 31-100 for baseline)
            idx = ~isnan(dir_rot{j}(:,i,k));
            
            % optimization
            log_likelihood = @(params) calc_likelihood_unif(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k));
            weightInit = 0.8;
            [weight_opt, neg_log_likelihood] = fmincon(log_likelihood, weightInit, [], [], [], [], 0, 1);
            
            log_likelihood = -neg_log_likelihood;
            
            % store fitted parameters and BIC
            unifWeight{j}(i,k) = weight_opt;
            unif_BIC{j}(i,k) = length(weightInit) * log(sum(idx)) - 2 * log_likelihood;
        end
    end
end

%% check whether there's extinction of habit in 1st vs 2nd half of flip block

for k = 1:2 % loop over first/second half of data
    idx = (k-1)*50 + (1:50); % indices of first or second half of trials
    for j = 1:Ngroup % loop over groups
        for i = 1:allSubj(j) % loop over subjects
            
            % optimization
            log_likelihood = @(params) calc_likelihood(params, dir_rot{j}(idx,i,4), target_gd{j}(idx,4), target_hab{j}(idx,4));            
            paramsInit = [weight1_init weight2_init kappa_init];
            params_opt = fmincon(log_likelihood, paramsInit, A, b);
            
            % store weight
            weight2_opt_half{j}(i,k) = params_opt(2);
        end
    end
end

% store data for statistical analysis in R
y = [weight2_opt_half{1}(:); weight2_opt_half{2}(:)];
groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
half([1:13 27:40],1) = "first";
half([14:26 41:54],1) = "second";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2])]';
T = table(groupNames, half, subject, y, 'VariableNames', {'group','half','subject','habit'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/half.csv')

%% Figure 4C
f = figure(8); clf;
set(f,'Position',[200 200 350 150]);

subplot(1,2,1); hold on
for i = 1:3
    n = allSubj(i);
    plot(repmat((i-1) + [1 5],[n 1]) + 0.5*(rand(n,2) - 0.5), weight1_opt{i}(:,3:4), '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot((i-1) + [1 5], mean(weight1_opt{i}(:,3:4)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'LineWidth', 1)
end
xticks([2 6])
xticklabels({'Late','Flip'})
xlabel('Block')
yticks(0:0.25:1)
ylabel('\alpha_a')
axis([0 8 0 1])
set(gca,'TickDir','out')

subplot(1,2,2); hold on
for i = 1:3
    n = allSubj(i);
    plot(repmat((i-1) + [1 5],[n 1]) + 0.5*(rand(n,2) - 0.5), weight2_opt{i}(:,3:4), '.', 'MarkerSize', 12, 'Color', col(i,:), 'HandleVisibility', 'off')
    plot((i-1) + [1 5], mean(weight2_opt{i}(:,3:4)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', col(i,:), 'LineWidth', 1)
end
xticks([2 6])
xticklabels({'Late','Flip'})
xlabel('Block')
yticks(0:0.25:1)
ylabel('\alpha_m')
axis([0 8 0 1])
set(gca,'TickDir','out')

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/habitMLE','-dpdf','-painters')

% save data for analysis in R
z = [];
for i = 1:2
    z = [z; reshape(weight2_opt{i}(:,3:4), [allSubj(i)*2 1])];
end
groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
blockNames([1:13 27:40],1) = "Late";
blockNames([14:26 41:54],1) = "Flip";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2])]';
T = table(groupNames, blockNames, subject, z, 'VariableNames', {'group','block','subject','habit'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/habit_weight.csv')

%% Figure 4D

% split reaction times by whether they were from goal-directed or habitual
% reaches
for j = 1:Ngroup
    for i = 1:allSubj(j)
        
        RT = data.(groups{j}){i}.RT(end-99:end); % data from flip block
        
        % divide data into goal-directed or habitual
        RT_gd{j}(i) = mean(RT(Pr_idx{j}(:,i) == 1),'omitnan');
        RT_hab{j}(i) = mean(RT(Pr_idx{j}(:,i) == 2),'omitnan');
        
        % store data from all trials
        RT_all{j}(:,i) = RT;
    end
end

f = figure(9); clf; hold on
set(f,'Position',[200 200 150 130]);
for j = 1:Ngroup
    plot(j + 0.5 * (rand(allSubj(j),1) - 0.5), RT_gd{j}, '.', 'Color', col(j,:), 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(j, mean(RT_gd{j}), 'ko', 'MarkerFaceColor', col(j,:), 'MarkerSize', 6, 'LineWidth', 1, 'HandleVisibility', 'off')
    
    plot(j + 4 + 0.5 * (rand(allSubj(j),1) - 0.5), RT_hab{j}, '.', 'Color', col(j,:), 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(j + 4, mean(RT_hab{j}), 'ko', 'MarkerFaceColor', col(j,:), 'MarkerSize', 6, 'LineWidth', 1)
end
xticks([2 6])
xticklabels({'Goal-directed', 'Habitual'})
xlim([0.5 7.5])
yticks(0.3:0.3:1.2)
ylim([0.3 1.2])
ylabel('Reaction time (s)')
set(gca,'Tickdir','out')

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/habit_RT','-dpdf','-painters')

y = [RT_gd{1}'; RT_hab{1}'; RT_gd{2}'; RT_hab{2}'];
groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
reach([1:13 27:40],1) = "gd";
reach([14:26 41:54],1) = "habit";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2])]';
T = table(groupNames, reach, subject, y, 'VariableNames', {'group','reach','subject','RT'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/RT.csv')

%% Figure 4E
f = figure(10); clf; hold on
set(f,'Position',[200 200 140 130]);
for i = 1:3
    plot(1:2, weight2_opt_half{i}','Color',[col(i,:) 0.5],'HandleVisibility','off')
end
for i = 1:3
    plot(1:2, mean(weight2_opt_half{i}),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'MarkerSize',6,'LineWidth',1)
end
set(gca,'TickDir','out')
xlim([0.75 2.25])
xticks([1 2])
xticklabels({'first half','second half'})
ylabel({'Weight'})

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/p2p_half','-dpdf','-painters')

%% Supplementary Figure 1B

f = figure(11); clf; hold on
set(f,'Position',[200 200 200 150]);
for i = 1:Ngroup
    plot(0.5*(rand(allSubj(i),1) - 0.5) + i, mix_BIC{i}(:,4), '.', 'Color', col(i,:), 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(i, mean(mix_BIC{i}(:,4)), 'ok', 'MarkerFaceColor', col(i,:), 'MarkerSize', 6, 'LineWidth', 1)
    
    plot(0.5*(rand(allSubj(i),1) - 0.5) + i + 4, unif_BIC{i}(:,4), '.', 'Color', col(i,:), 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(i + 4, mean(unif_BIC{i}(:,4)), 'ok', 'MarkerFaceColor', col(i,:), 'MarkerSize', 6, 'LineWidth', 1)
end
xticks([2 6])
xticklabels({'VM', 'unif'})
xlabel('Model')
ylim([150 400])
ylabel('BIC')
set(gca,'Tickdir','out')

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/BIC','-dpdf','-painters')

end

%% functions for calculating likelihoods
% calculate likelihood for von Mises mixture model
function neg_log_likelihood = calc_likelihood(params, samples, target_gd, target_hab)
    pdf = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
    
    kappa = params(3);
    weightUnif = 1 - sum(params(1:2));
    
    likelihood_vm1 = params(1) * pdf(samples, target_gd, kappa);
    likelihood_vm2 = params(2) * pdf(samples, target_hab, kappa);
    likelihood_unif = repelem(weightUnif/(2*pi),length(samples))';
    
    likelihood_all = sum([likelihood_vm1 likelihood_vm2 likelihood_unif],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end

% calculate likelihood for uniform mixture model
function neg_log_likelihood = calc_likelihood_unif(weight, samples, target_gd)
    
    unifPDF = 1/pi;
    idx = (target_gd >= 0) == (samples >= 0);
    
    likelihood = [idx ~idx] .* unifPDF;
    likelihood(:,1) = weight .* likelihood(:,1);
    likelihood(:,2) = (1-weight) .* likelihood(:,2);
    likelihood = sum(likelihood, 2);
    neg_log_likelihood = -sum(log(likelihood));
end