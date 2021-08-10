% rotate reach directions into quadrant 1 and plot histograms
function plot_habitMLE(d)

rng(2);
load variables/vmFit
groups = fieldnames(d);
graph_names = {'2-day','5-day','10-day'};
Ngroup = length(groups);
allSubj = [length(d.day2) length(d.day5) length(d.day10)];
tolerance = 0.01; % tolerance for stopping EM loop
xAxis = [0.5 2.5 5 15
         1 3 8 15.5
         1.5 3.5 13 16];

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
col = lines;
col = col(1:7,:);
col2 = [180 180 0
        0 191 255
        255 99 71]./255;

clear target_gd target_hab dir_rot 

% store all reach direction errors into dirError
for i = 1:Ngroup
    dir_rot{i} = NaN(100, allSubj(i), 4);
    target_gd{i} = NaN(100, 4);
    target_hab{i} = NaN(100, 4);

    for k = 1:Nblock
        Ntrials = length(d.(groups{i}){1}.Cr);
        Nsubj = length(d.(groups{i}));
        dir = NaN(Ntrials,Nsubj);
        for j = 1:Nsubj
            dir(:,j) = d.(groups{i}){j}.initDir_noRot;
        end

        dir2 = dir(trials{gblocks{i}(k)},:);
        targ = d.(groups{i}){1}.targetRel(trials{gblocks{i}(k)},:);
        
        dir2 = atan2(sin(dir2), cos(dir2));

        if k == 1
            idx = 1:30;
        else
            idx = 1:100;
        end
        
        target_gd{i}(idx,k) = atan2(targ(:,2), targ(:,1));
        target_hab{i}(idx,k) = atan2(targ(:,2), -targ(:,1));
        dir_rot{i}(idx,:,k) = dir2;
    end
end

% 
% for i = 1:Ngroup
%     for j = 1:allSubj(i)
%         for k = 1:4
%             up = target_gd{i}(:,k) >= 0;
%             down = target_gd{i}(:,k) < 0;
%             
%             up_correct = dir_rot{i}(up,j,k) >= 0;
%             down_correct = dir_rot{i}(down,j,k) < 0;
%             
%             correctCount = sum([up_correct; down_correct]);
%             
%             unif_likelihood{i}(j,k) = correctCount*unif;
%         end
%     end
% end

% % standard MLE
% for k = 1:Nblock
%     for j = 1:Ngroup
%         for i = 1:allSubj(j)
%             idx = ~isnan(dir_rot{j}(:,i,k));
%             
%             log_likelihood = @(params) calc_likelihood_unif(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k));
%             
%             weightInit = 0.8;
%             [weight_opt, neg_log_likelihood] = fmincon(log_likelihood, weightInit, [], [], [], [], 0, 1);
%             
%             log_likelihood = -neg_log_likelihood;
%             
%             unifWeight{j}(i,k) = weight_opt;
%             unif_likelihood{j}(i,k) = log_likelihood;
%             unif_AIC{j}(i,k) = 2 * length(weightInit) - 2 * log_likelihood;
%             unif_BIC{j}(i,k) = length(weightInit) * log(sum(idx)) - 2 * log_likelihood;
%         end
%     end
% end

unifPDF = 1/pi;
weightInit = 0.8;
% EM
for k = 1:Nblock
    for j = 1:Ngroup
        for i = 1:allSubj(j)
            
            weight = weightInit;
            idx = 1; % index for tracking number of EM iterations
            proceed = true; % flag for stopping EM loop
            
            id = ~isnan(dir_rot{j}(:,i,k));
            samples = dir_rot{j}(id,i,k);
            gd = target_gd{j}(id,k);
            same = (gd >= 0) == (samples >= 0);
            
%             log_likelihood = @(params) calc_likelihood_unif_em(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k));
%             
%             weightInit = 0.8;
%             [weight_opt, neg_log_likelihood] = fmincon(log_likelihood, weightInit, [], [], [], [], 0, 1);
            
            % EM loop
            while proceed
                
                denom = weight * (same * unifPDF) + (1 - weight) * (~same * unifPDF);
                Pr_unif1 = weight * (same * unifPDF) ./ denom;
                Pr_unif2 = 1 - Pr_unif1;
                
                resp = [Pr_unif1 Pr_unif2];
                
                % maximization step
                log_likelihood = @(params) calc_likelihood_unif_em(params, samples, gd, resp);
                paramsInit = weight; % set parameters to current values of mu and kappa
                [weight, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, [], [], [], [], 0, 1);
                
                % keep track of log-likelihood to terminate EM loop
                history(idx) = neg_log_likelihood;
                
                % terminate loop if change in log-likelihood is smaller
                % than tolerance,
                if idx > 1 && abs(history(idx) - history(idx-1)) < tolerance
                    proceed = false;
                end
                idx = idx + 1; % increment loop iteration number
            end
            
            log_likelihood = -neg_log_likelihood;
            
            unifWeight{j}(i,k) = weight;
            unif_likelihood{j}(i,k) = log_likelihood;
            unif_AIC{j}(i,k) = 2 * length(weightInit) - 2 * log_likelihood;
            unif_BIC{j}(i,k) = length(weightInit) * log(sum(idx)) - 2 * log_likelihood;
        end
    end
end

%% fit model using fixed kappa
% this isn't right because it's fitting model to random data. it's not
% testing parameter recovery at all; commenting this section out until it's
% needed

% Ntrials = 50;
% 
% up_sim = rand(50, Ntrials) * pi;
% down_sim = rand(50, Ntrials) * -pi;
% test = [up_sim; down_sim];
% 
% A = [-1 0; 0 -1; 1 1];
% b = [0 0 1]';
% 
% up_sim = rand(50, 1) * pi;
% down_sim = rand(50, 1) * -pi;
% targets = [up_sim; down_sim];
% targ_dir = [cos(targets) sin(targets)];
% targ_dir(:,1) = -targ_dir(:,1);
% targets2 = atan2(targ_dir(:,2), targ_dir(:,1));
% 
% weight1 = 0.5;
% weight2 = 0.5;
% kappa = 0:2:20;
% 
% weight1_sim = NaN(Ntrials, length(kappa));
% weight2_sim = NaN(Ntrials, length(kappa));
% weight3_sim = NaN(Ntrials, length(kappa));
% 
% for j = 1:length(kappa)
%     for i = 1:Ntrials
%         log_likelihood = @(params) calc_likelihood2(params, test(:,i), targets, targets2, kappa(j));
%         
%         paramsInit = [weight1 weight2];
%         params_opt = fmincon(log_likelihood, paramsInit, A, b);
%         
%         weight1_sim(i,j) = params_opt(1);
%         weight2_sim(i,j) = params_opt(2);
%         weight3_sim(i,j) = 1 - sum(params_opt);
%     end
% end
% 
% figure(7); clf; hold on
% plot(kappa, mean(weight1_sim))
% plot(kappa, mean(weight2_sim))
% plot(kappa, mean(weight3_sim))
% legend({'P(GD)','P(habit)','P(random)'}, 'Location', 'southeast')
% xlabel('kappa')
% ylabel('probability')

%% MLE on random data

% % bounded uniform distribution
% for i = 1:Ngroup
%     dir_sim{i} = NaN(100, allSubj(i));
%     
%     up = target_gd{i}(:,4) >= 0;
%     down = target_gd{i}(:,4) < 0;
%     
%     up_sim = rand(sum(up), allSubj(i)) * pi;
%     down_sim = rand(sum(down), allSubj(i)) * -pi;
%     
%     dir_sim{i}(up,:) = up_sim;
%     dir_sim{i}(down,:) = down_sim;
% end
% 
% % uniform distribution
% % dir_sim{1} = (rand(100, 13)-0.5) * 2 * pi;
% % dir_sim{2} = (rand(100, 14)-0.5) * 2 * pi;
% % dir_sim{3} = (rand(100, 5)-0.5) * 2 * pi;
% 
% weight1 = 0.33;
% weight2 = 0.33;
% kappaInit = 15;
% 
% for j = 1:Ngroup
%     for i = 1:allSubj(j)
%         log_likelihood = @(params) calc_likelihood(params, dir_sim{j}(:,i), target_gd{j}(:,4), target_hab{j}(:,4));
%         
%         paramsInit = [weight1 weight2 kappaInit];
%         A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
%         b = [0 0 1 200 0]';
%         [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
%         
%         weight1_sim{j}(i) = params_opt(1);
%         weight2_sim{j}(i) = params_opt(2);
%         weight3_sim{j}(i) = 1 - sum(params_opt(1:2));
%         kappa_sim{j}(i) = params_opt(3);
%         mix_likelihood_sim{j}(i) = -neg_log_likelihood;
%     end
% end
% 
% figure(6); clf
% subplot(1,3,1); hold on
% for i = 1:3
%     n = allSubj(i);
%     plot(i + 0.5*(rand(n,1) - 0.5), weight1_sim{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
%     plot(i, mean(weight1_sim{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
% end
% xlabel('Day')
% axis([0 4 0 1])
% ylabel('P(GD)')
% set(gca,'TickDir','out')
% legend(graph_names, 'Location', 'northwest')
% 
% subplot(1,3,2); hold on
% for i = 1:3
%     n = allSubj(i);
%     plot(i + 0.5*(rand(n,1) - 0.5), weight2_sim{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
%     plot(i, mean(weight2_sim{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
% end
% xlabel('Day')
% axis([0 4 0 1])
% ylabel('P(habitual)')
% set(gca,'TickDir','out')
% 
% subplot(1,3,3); hold on
% for i = 1:3
%     n = allSubj(i);
%     plot(i + 0.5*(rand(n,1) - 0.5), weight3_sim{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
%     plot(i, mean(weight3_sim{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
% end
% xlabel('Day')
% axis([0 4 0 1])
% ylabel('P(random)')
% set(gca,'TickDir','out')

%% perform MLE
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa)));
% weight1 = 0.33;
% weight2 = 0.33;
% kappaInit = 15;

weight1_init = 0.33;
weight2_init = 0.33;
kappa_init = 15;

A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
b = [0 0 1 200 0]';

for k = 1:4
    for j = 1:Ngroup
        for i = 1:allSubj(j)
%             if k < 4
%                 kappa = vmFit.kappa_opt(k,j,i);
%             else
%                 kappa = vmFit.kappa_opt(3,j,i);
%             end

%             idx = ~isnan(dir_rot{j}(:,i,k));
%             
%             log_likelihood = @(params) calc_likelihood(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k), target_hab{j}(idx,k));
%             
%             paramsInit = [weight1 weight2 kappaInit];
%             A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
%             b = [0 0 1 200 0]';
%             [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
            
            weight1 = weight1_init;
            weight2 = weight2_init;
            kappa = kappa_init;
            idx = 1; % index for tracking number of EM iterations
            proceed = true; % flag for stopping EM loop
            
            id = ~isnan(dir_rot{j}(:,i,k));
            samples = dir_rot{j}(id,i,k);
            gd = target_gd{j}(id,k);
            hab = target_hab{j}(id,k);
            
            % EM loop
            while proceed
                
                denom = weight1 * vmPDF(samples, gd, kappa) + weight2 * vmPDF(samples, hab, kappa) + (1-weight1-weight2) * (1 / (2*pi));
                Pr_gd = weight1 * vmPDF(samples, gd, kappa) ./ denom;
                Pr_hab = weight2 * vmPDF(samples, hab, kappa) ./ denom;
                Pr_unif = 1 - sum([Pr_gd Pr_hab],2);
                
                resp = [Pr_gd Pr_hab Pr_unif];
                
                % maximization step
                log_likelihood = @(params) calc_likelihood_em(params, samples, gd, hab, resp);
                paramsInit = [weight1 weight2 kappa]; % set parameters to current values of mu and kappa
                [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
                
                % assign optimized values of parameters
                weight1 = params_opt(1);
                weight2 = params_opt(2);
                kappa = params_opt(3);
                
                % keep track of log-likelihood to terminate EM loop
                history(idx) = neg_log_likelihood;
                
                % terminate loop if change in log-likelihood is smaller
                % than tolerance,
                if idx > 1 && abs(history(idx) - history(idx-1)) < tolerance
                    proceed = false;
                end
                idx = idx + 1; % increment loop iteration number
            end
            
            log_likelihood = -neg_log_likelihood;
            
            weight1_opt{j}(i,k) = weight1;
            weight2_opt{j}(i,k) = weight2;
            weight3_opt{j}(i,k) = 1 - weight1 - weight2;
            kappa_opt{j}(i,k) = kappa;
            mix_likelihood{j}(i,k) = log_likelihood;
            mix_AIC{j}(i,k) = 2 * length(paramsInit) - 2 * log_likelihood;
            mix_BIC{j}(i,k) = length(paramsInit) * log(sum(idx)) - 2 * log_likelihood;
            
            if k == 4
                dat = dir_rot{j}(:,i,k);
                mix1 = weight1 * vmPDF(dat, target_gd{j}(:,k), kappa);
                mix2 = weight2 * vmPDF(dat, target_hab{j}(:,k), kappa);
                mix3 = (1 - weight1 - weight2) * repelem(1/(2*pi), length(dat))';
                total = mix1 + mix2 + mix3;
                
                Pr_gd = mix1 ./ total;
                Pr_hab = mix2 ./ total;
                Pr_rand = mix3 ./ total;
                
                Pr_all = [Pr_gd Pr_hab Pr_rand];
                [~, idx] = max(Pr_all, [], 2);
                idx(isnan(dat)) = NaN;
                
                Pr_idx{j}(:,i) = idx;
            end
        end
    end
end

y = [unif_BIC{1}(:,4); unif_BIC{2}(:,4); unif_BIC{3}(:,4); mix_BIC{1}(:,4); mix_BIC{2}(:,4); mix_BIC{3}(:,4)];
model(1:32,1) = {'uniform'};
model(33:64,1) = {'mix'};

group([1:13 33:45],1) = {'2-day'};
group([14:27 46:59],1) = {'5-day'};
group([28:32 60:64],1) = {'10-day'};
% [p, tbl, stats] = anovan(y,{model group},'model',2,'varnames',{'model','group'});

% results = multcompare(stats,'Dimension',[1 2],'CType','hsd')

for i = 1:Ngroup
    weight2_frac{i} = weight2_opt{i} ./ (weight1_opt{i} + weight2_opt{i});
end

% BELOW IS CODE FOR PARAMETER RECOVERY
% t_q1 = pi/4;
% t_q2 = 3*pi/4;
% sigma = 0.9;
% dat = normrnd(t_q1, sigma, [100 1]);
% % dat = (rand(100,1)-0.5) * 2 * pi;
% % dat = [normrnd(t_q1, sigma, [50 1]); (rand(100,1)-0.5) * 2 * pi];
% 
% kappa = 4;
% log_likelihood = @(params) calc_likelihood(params, dat, t_q1, t_q2, kappa);
% 
% weightsInit = [weight1 weight2];
% A = [-1 0; 0 -1; 1 1];
% b = [0 0 1]';
% [params_opt, fval] = fmincon(log_likelihood, weightsInit, A, b, [], [], [0, 0], [1, 1]);
% 
% weight1_opt = params_opt(1);
% weight2_opt = params_opt(2);
% weight3_opt = 1 - sum(params_opt);

% dlmwrite('C:/Users/Chris/Documents/R/habit/data/sd.csv', z)

%% compare  kinematics between reaches towards mirrored and actual targets

for j = 1:Ngroup
    for i = 1:allSubj(j)
        vel = d.(groups{j}){i}.initVel_filt(end-99:end);
        velX = d.(groups{j}){i}.initVel_x(end-99:end);
        RT = d.(groups{j}){i}.RT(end-99:end);
        
        vel_gd{j}(i) = mean(vel(Pr_idx{j}(:,i) == 1));
        vel_hab{j}(i) = mean(vel(Pr_idx{j}(:,i) == 2));
        velX_gd{j}(i) = nanmean(abs(velX(Pr_idx{j}(:,i) == 1)));
        velX_hab{j}(i) = nanmean(abs(velX(Pr_idx{j}(:,i) == 2)));
        RT_gd{j}(i) = mean(RT(Pr_idx{j}(:,i) == 1));
        RT_hab{j}(i) = mean(RT(Pr_idx{j}(:,i) == 2));
    end
end

% z = [vel_gd{1}'; vel_gd{2}'; vel_gd{3}'; vel_hab{1}'; vel_hab{2}'; vel_hab{3}'];
% dlmwrite('C:/Users/Chris/Documents/R/habit/data/vel.csv', z)
% 
% z = [RT_gd{1}'; RT_gd{2}'; RT_gd{3}'; RT_hab{1}'; RT_hab{2}'; RT_hab{3}'];
% dlmwrite('C:/Users/Chris/Documents/R/habit/data/RT.csv', z)

figure(11); clf
subplot(1,3,1); hold on
for j = 1:Ngroup
    plot(j + 0.5 * (rand(allSubj(j),1) - 0.5), vel_gd{j}, '.', 'Color', col2(j,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(j, mean(vel_gd{j}), 'ko', 'MarkerFaceColor', col2(j,:), 'MarkerSize', 10, 'HandleVisibility', 'off')
    
    plot(j + 4 + 0.5 * (rand(allSubj(j),1) - 0.5), vel_hab{j}, '.', 'Color', col2(j,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(j + 4, mean(vel_hab{j}), 'ko', 'MarkerFaceColor', col2(j,:), 'MarkerSize', 10)
end
xticks([2 6])
xticklabels({'Actual', 'Mirrored'})
xlim([0.5 7.5])
yticks(0.05:0.1:0.35)
ylabel('Velocity (m/s)')
set(gca,'Tickdir','out')
legend(graph_names)

subplot(1,3,2); hold on
for j = 1:Ngroup
    plot(j + 0.5 * (rand(allSubj(j),1) - 0.5), velX_gd{j}, '.', 'Color', col2(j,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(j, mean(velX_gd{j}), 'ko', 'MarkerFaceColor', col2(j,:), 'MarkerSize', 10, 'HandleVisibility', 'off')
    
    plot(j + 4 + 0.5 * (rand(allSubj(j),1) - 0.5), velX_hab{j}, '.', 'Color', col2(j,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(j + 4, mean(velX_hab{j}), 'ko', 'MarkerFaceColor', col2(j,:), 'MarkerSize', 10)
end
xticks([2 6])
xticklabels({'Actual', 'Mirrored'})
xlim([0.5 7.5])
% yticks(0.05:0.1:0.35)
ylabel('Velocity (m/s)')
set(gca,'Tickdir','out')
legend(graph_names)

subplot(1,3,3); hold on
for j = 1:Ngroup
    plot(j + 0.5 * (rand(allSubj(j),1) - 0.5), RT_gd{j}, '.', 'Color', col2(j,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(j, mean(RT_gd{j}), 'ko', 'MarkerFaceColor', col2(j,:), 'MarkerSize', 10, 'HandleVisibility', 'off')
    
    plot(j + 4 + 0.5 * (rand(allSubj(j),1) - 0.5), RT_hab{j}, '.', 'Color', col2(j,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    plot(j + 4, mean(RT_hab{j}), 'ko', 'MarkerFaceColor', col2(j,:), 'MarkerSize', 10)
end
xticks([2 6])
xticklabels({'Actual', 'Mirrored'})
xlim([0.5 7.5])
yticks(300:200:1100)
ylabel('RT (ms)')
set(gca,'Tickdir','out')
legend(graph_names)

%% compare fits from flip block with fits from 
figure(1); clf
subplot(1,2,1); hold on
plot([0 100], [0 100], 'k')
for i = 1:Ngroup
    sd_em = permute(vmFit.sd(:,i,1:allSubj(i)), [3 1 2]) * 180/pi;
    R = besseli(1,kappa_opt{i}(:,1:3)) ./ besseli(0,kappa_opt{i}(:,1:3));
    sd_mle = sqrt(-2 * log(R)) * 180/pi; % circular standard deviation
    
    plot(sd_em, sd_mle, '.', 'Color', col2(i,:), 'MarkerSize', 20)
end
axis equal
xlabel('St dev from skill analysis')
ylabel('St dev from habit analysis')

subplot(1,2,2); hold on
for i = 1:Ngroup
    R = besseli(1,kappa_opt{i}) ./ besseli(0,kappa_opt{i});
    sd_mle = sqrt(-2 * log(R)) * 180/pi; % circular standard deviation
        
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), sd_mle, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(sd_mle), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 90])
ylabel('St dev from habit analysis')
set(gca,'TickDir','out')

%% compare log likelihoods
blocks2 = {'Baseline','Early','Late','Flip'};

figure(2); clf
for j = 1:4
    subplot(3,4,j); hold on
    for i = 1:Ngroup
        plot(0.5*(rand(allSubj(i),1) - 0.5) + i, mix_likelihood{i}(:,j), '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
        plot(i, mean(mix_likelihood{i}(:,j)), 'ok', 'MarkerFaceColor', col2(i,:), 'MarkerSize', 10)

        plot(0.5*(rand(allSubj(i),1) - 0.5) + i + 4, unif_likelihood{i}(:,j), '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
        plot(i + 4, mean(unif_likelihood{i}(:,j)), 'ok', 'MarkerFaceColor', col2(i,:), 'MarkerSize', 10, 'HandleVisibility', 'off')
    end
    title(blocks2(j))
    xticks([2 6])
    xticklabels({'VM', 'unif'})
    ylim([-200 50])
    if j == 1
        ylabel('log likelihood')
    end
end
legend(graph_names)

for j = 1:4
    subplot(3,4,j+4); hold on
    for i = 1:Ngroup
        plot(0.5*(rand(allSubj(i),1) - 0.5) + i, mix_AIC{i}(:,j), '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
        plot(i, mean(mix_AIC{i}(:,j)), 'ok', 'MarkerFaceColor', col2(i,:), 'MarkerSize', 10)
        
        plot(0.5*(rand(allSubj(i),1) - 0.5) + i + 4, unif_AIC{i}(:,j), '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
        plot(i + 4, mean(unif_AIC{i}(:,j)), 'ok', 'MarkerFaceColor', col2(i,:), 'MarkerSize', 10)
    end
    xticks([2 6])
    xticklabels({'VM', 'unif'})
    ylim([-100 400])
    if j == 1
        ylabel('AIC')
    end
end

for j = 1:4
    subplot(3,4,j+8); hold on    
    for i = 1:Ngroup
        plot(0.5*(rand(allSubj(i),1) - 0.5) + i, mix_BIC{i}(:,j), '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
        plot(i, mean(mix_BIC{i}(:,j)), 'ok', 'MarkerFaceColor', col2(i,:), 'MarkerSize', 10)

        plot(0.5*(rand(allSubj(i),1) - 0.5) + i + 4, unif_BIC{i}(:,j), '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
        plot(i + 4, mean(unif_BIC{i}(:,j)), 'ok', 'MarkerFaceColor', col2(i,:), 'MarkerSize', 10)
    end
    xticks([2 6])
    xticklabels({'VM', 'unif'})
    ylim([-100 400])
    if j == 1
        ylabel('BIC')
    end
end



%% plot P(habitual)
figure(3); clf
subplot(1,3,1); hold on
for i = 1:3
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), weight1_opt{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(weight1_opt{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 1])
ylabel('P(GD)')
set(gca,'TickDir','out')

subplot(1,3,2); hold on
for i = 1:3
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), weight2_opt{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(weight2_opt{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 1])
ylabel('\alpha_m')
set(gca,'TickDir','out')
legend(graph_names)

subplot(1,3,3); hold on
for i = 1:3
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), weight3_opt{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(weight3_opt{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 1])
ylabel('P(noise)')
set(gca,'TickDir','out')
legend(graph_names)

figure(4); clf; hold on
for i = 1:3
%     weight2_frac = weight2_opt{i} ./ (weight1_opt{i} + weight2_opt{i});
% 
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), weight2_frac{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(weight2_frac{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 0.8])
ylabel('P(habitual) / [P(habitual) + P(GD)]')
set(gca,'TickDir','out')

%%
idx = [5 50];

figure(5); clf; hold on
for i = 1:3
    sd = squeeze(vmFit.sd(3,i,1:allSubj(i)))*180/pi;
    frac = weight2_frac{i}(:,4);
    
    plot(sd, frac, '.', 'MarkerSize', 20, 'Color', col2(i,:))
    p = polyfit(sd, frac, 1);
    plot(idx, p(1)*idx + p(2), 'Color', col2(i,:), 'HandleVisibility', 'off')
end
ylim([-0.1 0.8])
xlabel('Standard deviation of VM (late)')
ylabel('P(habitual) / [P(habitual) + P(GD)]')
legend(graph_names)

end

function neg_log_likelihood = calc_likelihood(params, samples, target_gd, target_hab)
    pdf = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
%     pdf = @(x, mu, sigma) (1/(sigma*sqrt(2*pi)) * exp(-0.5 * ((x - mu)/sigma).^2));
    
    kappa = params(3);
    weightUnif = 1 - sum(params(1:2));
    
    likelihood_vm1 = params(1) * pdf(samples, target_gd, kappa);
    likelihood_vm2 = params(2) * pdf(samples, target_hab, kappa);
    likelihood_unif = repelem(weightUnif/(2*pi),length(samples))';
    
    likelihood_all = sum([likelihood_vm1 likelihood_vm2 likelihood_unif],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end

function neg_log_likelihood = calc_likelihood2(params, samples, target_gd, target_hab, kappa)
    pdf = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
    
    weightUnif = 1 - sum(params);
    
    likelihood_vm1 = params(1) * pdf(samples, target_gd, kappa);
    likelihood_vm2 = params(2) * pdf(samples, target_hab, kappa);
    likelihood_unif = repelem(weightUnif/(2*pi),length(samples))';
    
    likelihood_all = sum([likelihood_vm1 likelihood_vm2 likelihood_unif],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end

function neg_log_likelihood = calc_likelihood_unif(weight, samples, target_gd)
    
    unifPDF = 1/pi;
    idx = (target_gd >= 0) == (samples >= 0);
    
    likelihood = [idx ~idx] .* unifPDF;
    likelihood(:,1) = weight .* likelihood(:,1);
    likelihood(:,2) = (1-weight) .* likelihood(:,2);
    likelihood = sum(likelihood, 2);
    neg_log_likelihood = -sum(log(likelihood));
end

function neg_log_likelihood = calc_likelihood_em(params, samples, target_gd, target_hab, resp)
    pdf = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
    
    kappa = params(3);
    weightUnif = 1 - sum(params(1:2));
    
    likelihood_vm1 = resp(:,1) .* (params(1) * pdf(samples, target_gd, kappa));
    likelihood_vm2 = resp(:,2) .* (params(2) * pdf(samples, target_hab, kappa));
    likelihood_unif = resp(:,3) .* (weightUnif / (2*pi));
    
    likelihood_all = sum([likelihood_vm1 likelihood_vm2 likelihood_unif],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end

function neg_log_likelihood = calc_likelihood_unif_em(weight, samples, target_gd, resp)
    
    unifPDF = 1/pi;
    idx = (target_gd >= 0) == (samples >= 0);
    
    likelihood = [idx ~idx] .* unifPDF;
    likelihood(:,1) = resp(:,1) .* (weight .* likelihood(:,1));
    likelihood(:,2) = resp(:,2) .* ((1-weight) .* likelihood(:,2));
    likelihood = sum(likelihood, 2);
    neg_log_likelihood = -sum(log(likelihood));
end