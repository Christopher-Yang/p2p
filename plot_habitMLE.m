% rotate reach directions into quadrant 1 and plot histograms
function plot_habitMLE(d)

rng(2);
groups = fieldnames(d);
graph_names = {'2-day','5-day','10-day'};
Ngroup = length(groups);
allSubj = [length(d.day2) length(d.day5) length(d.day10)];

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

%% perform MLE for von Mises model
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa)));

weight1_init = 0.33;
weight2_init = 0.33;
kappa_init = 15;

A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
b = [0 0 1 200 0]';

for k = 1:4
    for j = 1:Ngroup
        for i = 1:allSubj(j)

            idx = ~isnan(dir_rot{j}(:,i,k));
            
            log_likelihood = @(params) calc_likelihood(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k), target_hab{j}(idx,k));
            
            paramsInit = [weight1_init weight2_init kappa_init];
            [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
            
            weight1 = params_opt(1);
            weight2 = params_opt(2);
            kappa = params_opt(3);
            log_likelihood = -neg_log_likelihood;
            
            weight1_opt{j}(i,k) = weight1;
            weight2_opt{j}(i,k) = weight2;
            weight3_opt{j}(i,k) = 1 - weight1 - weight2;
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

%% check whether there's extinction of habit in 1st vs 2nd half of flip block
for k = 1:2
    idx = (k-1)*50 + (1:50);
    for j = 1:Ngroup
        for i = 1:allSubj(j)            
            
            log_likelihood = @(params) calc_likelihood(params, dir_rot{j}(idx,i,4), target_gd{j}(idx,4), target_hab{j}(idx,4));
            
            paramsInit = [weight1_init weight2_init kappa_init];
            params_opt = fmincon(log_likelihood, paramsInit, A, b);
            
            weight2 = params_opt(2);
            
            weight2_opt2{j}(i,k) = weight2;
        end
    end
end

y = [weight2_opt2{1}(:); weight2_opt2{2}(:); weight2_opt2{3}(:)];

groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
half([1:13 27:40 55:59],1) = "first";
half([14:26 41:54 60:64],1) = "second";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, half, subject, y, 'VariableNames', {'group','half','subject','habit'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/half.csv')

%% perform MLE for uniform model
for k = 1:Nblock
    for j = 1:Ngroup
        for i = 1:allSubj(j)
            idx = ~isnan(dir_rot{j}(:,i,k));
            
            log_likelihood = @(params) calc_likelihood_unif(params, dir_rot{j}(idx,i,k), target_gd{j}(idx,k));
            
            weightInit = 0.8;
            [weight_opt, neg_log_likelihood] = fmincon(log_likelihood, weightInit, [], [], [], [], 0, 1);
            
            log_likelihood = -neg_log_likelihood;
            
            unifWeight{j}(i,k) = weight_opt;
            unif_BIC{j}(i,k) = length(weightInit) * log(sum(idx)) - 2 * log_likelihood;
        end
    end
end

for i = 1:Ngroup
    weight2_frac{i} = weight2_opt{i} ./ (weight1_opt{i} + weight2_opt{i});
end
%% statistical comparison between von Mises and uniform model

y = [unif_BIC{1}(:,4); mix_BIC{1}(:,4); unif_BIC{2}(:,4); mix_BIC{2}(:,4); unif_BIC{3}(:,4); mix_BIC{3}(:,4)];

groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
modelType([1:13 27:40 55:59],1) = "unif";
modelType([14:26 41:54 60:64],1) = "vm";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, modelType, subject, y, 'VariableNames', {'group','model','subject','bic'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/habitBIC.csv')

%% Figure 4D
f = figure(4); clf;
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

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/habitMLE','-dpdf','-painters')

% save data for analysis in R
z = [];
for i = 1:3
    z = [z; reshape(weight2_opt{i}(:,3:4), [allSubj(i)*2 1])];
end
groupNames(1:26,1) = "2-day";
groupNames(27:54,1) = "5-day";
groupNames(55:64,1) = "10-day";
blockNames([1:13 27:40 55:59],1) = "Late";
blockNames([14:26 41:54 60:64],1) = "Flip";
subject = [repmat(1:13,[1 2]) repmat(14:27,[1 2]) repmat(28:32,[1 2])]';
T = table(groupNames, blockNames, subject, z, 'VariableNames', {'group','block','subject','habit'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/habit_weight1.csv')

z = [];
for i = 1:3
    z = [z; reshape(weight1_opt{i}(:,3:4), [allSubj(i)*2 1])];
end
T = table(groupNames, blockNames, subject, z, 'VariableNames', {'group','block','subject','habit'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/habit_weight2.csv')

%% Figure 4C

times = 0:10:1000;

for j = 1:Ngroup
    for i = 1:allSubj(j)
        RT = d.(groups{j}){i}.RT(end-99:end);
        
        RT_gd{j}(i) = mean(RT(Pr_idx{j}(:,i) == 1));
        RT_hab{j}(i) = mean(RT(Pr_idx{j}(:,i) == 2));
        
        RT_all{j}(:,i) = RT;
    end
    
    Nsubj = length(d.(groups{j}));
    frac{j} = NaN(length(times),Nsubj);
    
    for i = 1:length(times)
        high = RT_all{j} < times(i) + 100;
        low = RT_all{j} > times(i);
        bin = (high + low) == 2;        
        total = sum(bin,1);
        
        if total > 5
            habIdx = Pr_idx{j} == 2;
            habitBin = (habIdx + bin) == 2;

            frac{j}(i,:) = sum(habitBin,1)./total;
        end
    end

    frac_mu{j} = mean(frac{j},2,'omitnan');
    frac_se{j} = std(frac{j},[],2,'omitnan')./sqrt(sum(~isnan(frac{j}),2));

end

figure(8); clf; hold on
for j = 1:Ngroup
    s = shadedErrorBar(times+50,frac_mu{j},frac_se{j});
    editErrorBar(s,col(j,:),1);
end
ylim([0 1])
ylabel('P(habitual)')
xlabel('Reaction time (ms)')
legend(graph_names)

f = figure(7); clf; hold on
set(f,'Position',[200 200 150 130]);
for j = 1:Ngroup
    plot(j + 0.5 * (rand(allSubj(j),1) - 0.5), RT_gd{j}, '.', 'Color', col(j,:), 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(j, mean(RT_gd{j}), 'ko', 'MarkerFaceColor', col(j,:), 'MarkerSize', 6, 'LineWidth', 1, 'HandleVisibility', 'off')
    
    plot(j + 4 + 0.5 * (rand(allSubj(j),1) - 0.5), RT_hab{j}, '.', 'Color', col(j,:), 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(j + 4, mean(RT_hab{j}), 'ko', 'MarkerFaceColor', col(j,:), 'MarkerSize', 6, 'LineWidth', 1)
end
xticks([2 6])
xticklabels({'Actual', 'Mirrored'})
xlim([0.5 7.5])
yticks(300:300:900)
ylim([200 1000])
ylabel('Reaction time (ms)')
set(gca,'Tickdir','out')

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/habit_kinematics','-dpdf','-painters')

%% Figure 4E
f = figure(10); clf; hold on
set(f,'Position',[200 200 140 130]);
for i = 1:3
    plot(1:2, weight2_opt2{i}','Color',[col(i,:) 0.5],'HandleVisibility','off')
end
for i = 1:3
    plot(1:2, mean(weight2_opt2{i}),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'MarkerSize',6,'LineWidth',1)
end
set(gca,'TickDir','out')
xlim([0.75 2.25])
xticks([1 2])
xticklabels({'first half','second half'})
ylabel({'Weight'})

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/p2p_half','-dpdf','-painters')

%% Figure S2B

f = figure(20); clf; hold on
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

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/BIC','-dpdf','-painters')

end

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

function neg_log_likelihood = calc_likelihood_unif(weight, samples, target_gd)
    
    unifPDF = 1/pi;
    idx = (target_gd >= 0) == (samples >= 0);
    
    likelihood = [idx ~idx] .* unifPDF;
    likelihood(:,1) = weight .* likelihood(:,1);
    likelihood(:,2) = (1-weight) .* likelihood(:,2);
    likelihood = sum(likelihood, 2);
    neg_log_likelihood = -sum(log(likelihood));
end
