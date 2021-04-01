function plot_direction(d, vmFit)

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
% gblocks{1} = [1 2 5];
% gblocks{2} = [1 2 6 9 12 14];
% gblocks{3} = [1 2 6 9 12 15 18 21 24 27 29];
      
% set colors for generating plots
col = lines;
col = col(1:7,:);
col2 = [180 180 0
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

if nargin > 1
    mu_opt = vmFit.mu_opt;
    kappa_opt = vmFit.kappa_opt;
    weight_opt = vmFit.weight_opt;
    sd = vmFit.sd;
else
    tolerance = 0.01; % tolerance for stopping EM loop
    
    % set values for initial parameters for EM loop
    muInit = 0;
    kappaInit = 1;
    weightInit = 0.99;
    
    % preallocate variables for parameters
    mu_opt = NaN(Nblock,Ngroup,max(allSubj));
    kappa_opt = NaN(Nblock,Ngroup,max(allSubj));
    weight_opt = NaN(Nblock,Ngroup,max(allSubj));
    sd = NaN(Nblock,Ngroup,max(allSubj));
    
    % main loop for EM
    for k = 1:Ngroup
        Nsubj = length(d.(groups{k}));
        for m = 1:Nsubj
            for j = 1:Nblock
                
                % select trials to analyze and store in samples
                trial = trials{gblocks{k}(j)};
                samples = dirError{k}(trial,m)*pi/180; % convert error to radians
                samples = samples(~isnan(samples));
                
                % initialize mu and kappa for the VM distribution and the
                % relative weight between the VM and uniform distributions
                mu = muInit;
                kappa = kappaInit;
                weight = weightInit;
                idx = 1; % index for tracking number of EM iterations
                proceed = true; % flag for stopping EM loop
                
                % EM loop
                while proceed
                    
                    % expectation step
                    Pr_vm = weight * vmPDF(samples, mu, kappa) ./ (weight * vmPDF(samples, mu, kappa) + (1-weight) * (1 / (2*pi)));
                    Pr_unif = (1-weight) * (1 / (2*pi)) ./ (weight * vmPDF(samples, mu, kappa) + (1-weight) * (1 / (2*pi)));
                    
                    % maximization step
                    log_likelihood = @(params) calc_likelihood(params, samples, Pr_vm);
                    paramsInit = [mu kappa weight]; % set parameters to current values of mu and kappa
                    [params_opt, fval] = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 100 1]);
                    
                    % assign optimized values of parameters
                    mu = params_opt(1);
                    kappa = params_opt(2);
                    weight = params_opt(3);
                    
                    % keep track of log-likelihood to terminate EM loop
                    history{k}{m}{j}(idx) = fval;
                    
                    % terminate loop if change in log-likelihood is smaller
                    % than tolerance,
                    if idx > 1 && abs(history{k}{m}{j}(idx) - history{k}{m}{j}(idx-1)) < tolerance
                        proceed = false;
                    end
                    idx = idx + 1; % increment loop iteration number
                    
%                     % analytical approach to solve MLE
%                     xBar = mean(exp(1j*vmSamples));
%                     R = norm(xBar);
%                     mu = angle(xBar);
%                     kappa = R * (2 - R^2) ./ (1 - R^2);
%                     
%                     if idx > 1 && abs(history{k}{m}{j}(idx) - history{k}{m}{j}(idx-1)) < tolerance
%                         proceed = false;
%                     end
%                     idx = idx + 1;
                end
                
                % store fitted parameter values
                mu_opt(j,k,m) = mu;
                kappa_opt(j,k,m) = kappa;
                weight_opt(j,k,m) = weight;
                
                % compute circular standard deviation
                R = (besseli(1,kappa)/besseli(0,kappa));
                sd(j,k,m) = sqrt(-2 * log(R)); % circular standard deviation
            end
        end
    end
    
    vmFit.mu_opt = mu_opt;
    vmFit.kappa_opt = kappa_opt;
    vmFit.weight_opt = weight_opt;
    vmFit.sd = sd;
    
    save('variables/vmFit.mat','vmFit')
end

% points to assess the PDF
delt = pi/64;
x = -pi:delt:pi-delt;

% mean and standard deviation of the circular standard deviation
sd_mu = nanmean(sd,3);
sd_sd = nanstd(sd,[],3);

sd2 = permute(sd,[1 3 2]);

%% plot mixture model fit on top of data histograms
subj = 2; % choose which subject to plot fits for

figure(1); clf
for j = 1:Nblock
    for k = 1:Ngroup
        subplot(3, 3, (j-1)*3 + k); hold on
        Nsubj = size(dirError{k},2);
        
        % plot histograms
        if subj <= Nsubj
            trial = trials{gblocks{k}(j)};
            histogram(dirError{k}(trial,subj),x*180/pi,'Normalization','probability');
        end
        
        % plot pdf
        pdf = vmPDF(x, mu_opt(j,k,subj), kappa_opt(j,k,subj)); % PDF of von Mises distribution
        mixPDF = weight_opt(j,k,subj) * pdf + (1-weight_opt(j,k,subj)) * (1 / (2*pi)); % weight von Mises with uniform distribution
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

%% plot standard deviation of von Mises distribution (x-axis unscaled by practice time)
% figure(2); clf; hold on
% for i = 1:3
%     plot(repmat([0 5 10]', [1 14]) + (rand(3,14)-0.5) + 1.5*i, sd2(:,:,i)*180/pi, '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
%     plot([0 5 10] + 1.5*i, sd_mu(:,i)*180/pi, 'o', 'MarkerFaceColor', col2(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)
% end
% ylabel('St dev of von Mises distribution')
% xlim([0.5 15.5])
% xticks(3:5:18)
% yticks(0:20:80)
% xticklabels(blocks)
% legend(graph_names)
% set(gca,'Tickdir','out')

xAxis = [0.5 1 1.5
         2.5 3 3.5
         5 11 21];

% same plot as above but x-axis scaled by practice time
figure(2); clf; hold on
for i = 1:Ngroup
    for j = 1:Nblock
        plot(xAxis(j,i) + 0.5*(rand(1,14)-0.5), squeeze(sd2(j,:,i))*180/pi, '.', 'Color', col2(i,:), 'MarkerSize', 20, 'HandleVisibility', 'off')
    end
end
for i = 1:Ngroup
    plot(xAxis(:,i), sd_mu(:,i)*180/pi, 'o', 'MarkerFaceColor', col2(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)
end
xticks([1 3 5 11 21])
xticklabels({'Baseline','1','2','5','10'})
xlabel('Day')
ylabel('St dev of von Mises distribution')
xlim([0 22])
yticks(0:20:100)
legend(graph_names)
set(gca,'Tickdir','out')

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/stdev','-dpdf','-painters')

% plot same data as above but grouped by group, not time during learning
figure(3); clf; hold on
for i = 1:3
    plot((1:3) + (i-1)*4, sd2(:,:,i)*180/pi, 'Color', [col2(i,:) 0.6])
    plot((1:3) + (i-1)*4, sd_mu(:,i)*180/pi, '.', 'Color', col2(i,:), 'MarkerSize', 20)
end
ylabel('St dev of von Mises distribution')
xticks([1:3 5:7 9:11])
xticklabels([blocks blocks blocks])

%% plot kernel-smoothed PDF
figure(4); clf
for i = 1:Ngroup
    for j = 1:3
        trial = trials{gblocks{i}(j)};
        subplot(1,3,j); hold on
        [f,xi] = ksdensity(reshape(dirError{i}(trial,:),[numel(dirError{i}(trial,:)) 1]));
        plot(xi,f,'LineWidth',2,'Color',col2(i,:))
        if i == 3
            title(blocks{j})
            axis([-180 180 0 .06])
            xticks(-180:90:180)
            box off
            set(gca,'Tickdir','out')
            if j == 1
                ylabel('Kernel-smoothed probability density')
                yticks(0:0.02:0.06)
            elseif j == 2
                xlabel('Reach direction error (degrees)')
                yticks([])
            elseif j == 3
                yticks([])
            end
        end
    end
end
legend(graph_names)

% print('C:/Users/Chris/Dropbox/Conferences/CNS 2021/ksdensity','-dpdf','-painters')

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
    plot(sd2(3,:,i)*180/pi, habit(:,i), '.', 'Color', col2(i,:), 'MarkerSize', 30)
end

% plot best-fit lines from least-squares regression
p = polyfit(sd2(3,1:13,1)'*180/pi, habit(1:13,1), 1);
plot(idx, p(1)*idx + p(2), 'Color', col2(1,:))
p = polyfit(sd2(3,:,2)'*180/pi, habit(:,2), 1);
plot(idx, p(1)*idx + p(2), 'Color', col2(2,:))
p = polyfit(sd2(3,1:5,3)'*180/pi, habit(1:5,3), 1);
plot(idx, p(1)*idx + p(2), 'Color', col2(3,:))

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

% %% bootstrap confidence intervals for direction error
% 
% boot = NaN(1000,Nblock,Ngroup);
% for k = 1:Ngroup
%     Nsubj = length(d.(groups{k}));
%     for j = 1:Nblock
%         trial = trials{gblocks(k,j)};
%         nSamples = numel(dirError{k}(trial,:));
%         for i = 1:1000
%             sample = datasample(reshape(dirError{k}(trial,:),[nSamples 1]),13);
%             boot(i,j,k) = std(sample);
%         end
%     end
% end
% 
% col = lines;
% col = col(1:7,:);
% 
% figure(6); clf; hold on
% for j = 1:Ngroup
%     subplot(3,1,j); hold on
%     for i = 1:Nblock
% %         histogram(boot(:,i,j),0:2:120,'Normalization','pdf','FaceColor',col(i,:))
%         [f,xi] = ksdensity(boot(:,i,j));
%         plot(xi,f,'Color',col(i,:),'LineWidth',2)
%     end
%     if j == 2
%         ylabel('Probability')
%     end
%     ylim([0 0.3])
% end
% xlabel('Standard deviation of directional error (degrees)')
% legend(blocks)
% 
% boot = sort(boot,1);
% confInterval = boot([26 975],:,:);
% bootMu = mean(boot,1);
% confDiff = abs(confInterval - repmat(bootMu,[2 1 1]));
% bootMu = squeeze(bootMu);
% 
% figure(7); clf; hold on
% for i = 1:3
%     errorbar(i,bootMu(1,i),confDiff(1,1,i),confDiff(2,1,i),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'LineWidth',2)
%     errorbar(i+4,bootMu(2,i),confDiff(1,2,i),confDiff(2,2,i),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'LineWidth',2,'HandleVisibility','off')
%     errorbar(i+8,bootMu(3,i),confDiff(1,3,i),confDiff(2,3,i),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'LineWidth',2,'HandleVisibility','off')
% end
% xticks(2:4:10)
% xticklabels(blocks)
% xlim([0.5 11.5])
% ylabel('Standard deviation of directional error (degrees)')
% legend(graph_names)
% 
% %%
% Ntrials2 = 100;
% trialsAll = {1:100,101:200,201:300};
% for k = 1:3
%     trials = trialsAll{k};
%     for i = 1:Ngroup
%         Nsubj = length(d.(groups{i}));
%         for j = 1:Nsubj
%             ang = atan2(d.(groups{i}){j}.targetRel(trials,2), d.(groups{i}){j}.targetRel(trials,1));
%             angMir = atan2(d.(groups{i}){j}.targetRel(trials,2), -d.(groups{i}){j}.targetRel(trials,1));
% 
%             dir = d.(groups{i}){j}.initDir_noRot(trials)';
% 
%             error{i}(:,j) = dir-ang;
%             errorMir{i}(:,j) = dir-angMir;
%         end
%         
%         for l = 1:numel(error{i})
%             while error{i}(l) >= pi
%                 error{i}(l) = error{i}(l)-2*pi;
%             end
%             while error{i}(l) < -pi
%                 error{i}(l) = error{i}(l)+2*pi;
%             end
%             while errorMir{i}(l) >= pi
%                 errorMir{i}(l) = errorMir{i}(l)-2*pi;
%             end
%             while errorMir{i}(l) < -pi
%                 errorMir{i}(l) = errorMir{i}(l)+2*pi;
%             end
%         end
%         
%         towardMir{i}.all(:,k) = sum(abs(error{i})>abs(errorMir{i}),1);
%         towardMir{i}.mean(k) = mean(towardMir{i}.all(:,k),1);
%     end
% end
% 
% figure(9); clf; hold on
% for i = 1:3
%     plot(i,towardMir{1}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
%     plot(i,towardMir{1}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
%     
%     plot(i+4,towardMir{2}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
%     plot(i+4,towardMir{2}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
%     
%     plot(i+8,towardMir{3}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
%     plot(i+8,towardMir{3}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
% end
% xticks([2 6 10])
% xticklabels(graph_names)
% ylabel('Percent towards mirrored target')
% axis([0.5 11.5 0 70])
% box off
end

% function for computing log-likelihod
function neg_log_likelihood = calc_likelihood(params,samples,Pr_vm)
    mu = params(1);
    kappa = params(2);
    weight = params(3);
    
    likelihood_unif = (1 - Pr_vm) .* log(1 - weight);
    likelihood_vm = Pr_vm .* (log(weight) + log(exp(kappa * cos(samples-mu)) / (2 * pi * besseli(0,kappa))));
    
    likelihood_all = [likelihood_unif likelihood_vm];
    neg_log_likelihood = -sum(sum(likelihood_all,2),1);
end