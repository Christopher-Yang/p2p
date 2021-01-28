
% set variables for analysis
groups = fieldnames(d);
graph_names = {'2-day','5-day','10-day'};
names = {'day2','day5','day10'};
blocks = {'Baseline','Early','Late'};
Ngroup = length(groups);
Nblock = length(blocks);
gblocks = [1 2 3; 1 2 4; 1 2 5]; % blocks for baseline, early, and late for each group
allSubj = [length(d.day2) length(d.day5) length(d.day10)];

% indices for dividing up the trials into blocks
trials{1} = 1:30;
trials{2} = 31:130;
trials{3} = 131:230;
trials{4} = 231:330;
trials{5} = 331:430;
trials{6} = 431:530;

% set colors for generating plots
col = lines;
col = col(1:7,:);

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

%% fit mixture model for each group
clear mu_opt kappa_opt weight_opt sd history
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
tolerance = 0.01; % tolerance for stopping EM loop

% set values for initial parameters for EM loop
muInit = 0;
kappaInit = 1;
weightInit = 0.99;

% run EM for each group and block
for k = 1:Ngroup % loop over groups
    for j = 1:Nblock % loop over blocks
        
        % select trials to analyze and store in samples
        trial = trials{gblocks(k,j)};
        samples = dirError{k}(trial,:)*pi/180; % convert errors into radians
        samples = reshape(samples, [numel(samples) 1]); % convert matrix to vector
        
        % initialize mu and kappa for the VM distribution and the relative
        % weight between the VM and uniform distributions
        mu = muInit; % initial mu
        kappa = kappaInit; % initial kappa
        weight = weightInit; % initial weight
        idx = 1; % index for tracking iteration number
        proceed = true; % flag for continuing EM
        
        % EM loop
        while proceed
            
            % expectation step
            Pr_vm = weight * vmPDF(samples, mu, kappa) ./ (weight * vmPDF(samples, mu, kappa) + (1-weight) * (1 / (2*pi)));
            
            % maximization step
            log_likelihood = @(params) calc_likelihood(params, samples, Pr_vm); % likelihood function
            paramsInit = [mu kappa weight]; % set parameters to current values of mu and kappa
            [params_opt, fval] = fmincon(log_likelihood, paramsInit, [], [], [], [], [-pi 0 0], [pi 100 1]);
            
            % assign optimized values of parameters
            mu = params_opt(1);
            kappa = params_opt(2);
            weight = params_opt(3);
            
            % keep track of log-likelihood for EM loop termination
            history{k}{j}(idx) = fval;
            
            % check whether change in log-likelihood is below tolerance
            if idx > 1 && abs(history{k}{j}(idx) - history{k}{j}(idx-1)) < tolerance
                proceed = false;
            end
            idx = idx + 1; % increment iteration number

%             % analytical approach to solve MLE
%             xBar = mean(exp(1j*vmSamples));
%             R = norm(xBar);
%             mu = angle(xBar);
%             kappa = R * (2 - R^2) ./ (1 - R^2);
% 
%             if idx > 50
%                 proceed = false;
%             end
%             idx = idx + 1;
        end
        
        % store fitted parameters for each dataset
        mu_opt(j,k) = mu;
        kappa_opt(j,k) = kappa;
        weight_opt(j,k) = weight;
        
        % calculate circular standard deviation
        m1 = (besseli(1,kappa)/besseli(0,kappa))*exp(1j*mu); % first moment
        R = abs(m1);
        sd(j,k) = sqrt(-2 * log(R)); % standard deviation
    end
end

% points to assess the PDF
delt = pi/64;
x = -pi:delt:pi-delt;

% plot the PDFs of mixture model and just the VM distribution
figure(1); clf
for j = 1:Nblock
    
    % mixture model
    subplot(2,3,j); hold on 
    for k = 1:Ngroup
        pdf = vmPDF(x, mu_opt(j,k), kappa_opt(j,k));
        mixPDF = weight_opt(j,k) * pdf + (1-weight_opt(j,k)) * (1 / (2*pi));
        plot(x*180/pi, mixPDF, 'LineWidth', 2)
    end
    title(blocks{j})
    if j == 1
        ylabel({'Probability density (mixture)'})
    elseif j == 2
        xlabel({'Reach direction error (degrees)'})
    end
    axis([-180 180 0 3.5])
    
    % VM distribution
    subplot(2,3,j+3); hold on
    for k = 1:Ngroup
        pdf = vmPDF(x, mu_opt(j,k), kappa_opt(j,k));
        plot(x*180/pi, pdf, 'LineWidth', 2)
    end
    if j == 1
        ylabel({'Probability density (von Mises)'})
    elseif j == 2
        xlabel({'Reach direction error (degrees)'})
    end
    axis([-180 180 0 3.5])
end

% plot the mixture model on top of data histograms
figure(2); clf
for j = 1:Nblock
    for k = 1:Ngroup
        subplot(3, 3, (j-1)*3 + k); hold on
        
        % plot histograms
        trial = trials{gblocks(k,j)};
        histogram(dirError{k}(trial,:),x*180/pi,'Normalization','probability');
        
        % plot pdf
        pdf = vmPDF(x, mu_opt(j,k), kappa_opt(j,k));
        mixPDF = weight_opt(j,k) * pdf + (1-weight_opt(j,k)) * (1 / (2*pi));
        plot(x*180/pi, mixPDF./sum(mixPDF), 'LineWidth', 2)
        
        if j == 1
            title(graph_names{k})
        elseif j == 2 && k == 1
            ylabel('Probability')
        elseif j == 3 && k == 2
            xlabel('Reach Direction Error (degrees)')
        end
        axis([-180 180 0 0.2])
    end
end

% plot mean and standard deviation of the fitted standard deviation from
% the von Mises distribution
figure(3); clf; hold on
b = bar(sd*180/pi);
for i = 1:3
    b(i).FaceColor = col(i,:);
end
xticks(1:3)
xticklabels(blocks)
ylabel('St dev of von Mises distribution (degrees)')
legend(graph_names)

%% fit mixture model for each participant
clear history
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
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
            trial = trials{gblocks(k,j)};
            samples = dirError{k}(trial,m)*pi/180; % convert error to radians
            
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
                
%                 % analytical approach to solve MLE
%                 xBar = mean(exp(1j*vmSamples));
%                 R = norm(xBar);
%                 mu = angle(xBar);
%                 kappa = R * (2 - R^2) ./ (1 - R^2);
% 
%                 if idx > 50
%                     proceed = false;
%                 end
%                 idx = idx + 1;
            end
            
            % store fitted parameter values
            mu_opt(j,k,m) = mu;
            kappa_opt(j,k,m) = kappa;
            weight_opt(j,k,m) = weight;
            
            % compute circular standard deviation
            m1 = (besseli(1,kappa)/besseli(0,kappa))*exp(1j*mu); % first moment
            R = abs(m1);
            sd(j,k,m) = sqrt(-2 * log(R)); % circular standard deviation
        end
    end
end

% mean and standard deviation of the circular standard deviation
sd_mu = nanmean(sd,3);
sd_sd = nanstd(sd,[],3);

% plot the mean and standard deviation of the standard deviations
figure(4); clf; hold on
for i = 1:3
    errorbar([0 4 8]+i, sd_mu(:,i)*180/pi, sd_sd(:,i)*180/pi, '.', 'Color', col(i,:), 'LineWidth', 2, 'MarkerSize', 20)
end
ylabel('St dev of von Mises distribution')
xticks(2:4:10)
xticklabels(blocks)
legend(graph_names)

% plot PDF on top of data histograms
subj = 1;
figure(5); clf
for j = 1:Nblock
    for k = 1:Ngroup
        subplot(3, 3, (j-1)*3 + k); hold on
        Nsubj = size(dirError{k},2);
        
        % plot histograms
        if subj <= Nsubj
            trial = trials{gblocks(k,j)};
            histogram(dirError{k}(trial,subj),x*180/pi,'Normalization','probability');
        end
        
        % plot pdf
        pdf = vmPDF(x, mu_opt(j,k,subj), kappa_opt(j,k,subj));
        mixPDF = weight_opt(j,k,subj) * pdf + (1-weight_opt(j,k,subj)) * (1 / (2*pi));
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
%% bootstrap confidence intervals for direction error
rng(1);

boot = NaN(1000,Nblock,Ngroup);
for k = 1:Ngroup
    Nsubj = length(d.(groups{k}));
    for j = 1:Nblock
        trial = trials{gblocks(k,j)};
        nSamples = numel(dirError{k}(trial,:));
        for i = 1:1000
            sample = datasample(reshape(dirError{k}(trial,:),[nSamples 1]),13);
            boot(i,j,k) = std(sample);
        end
    end
end

col = lines;
col = col(1:7,:);

figure(6); clf; hold on
for j = 1:Ngroup
    subplot(3,1,j); hold on
    for i = 1:Nblock
%         histogram(boot(:,i,j),0:2:120,'Normalization','pdf','FaceColor',col(i,:))
        [f,xi] = ksdensity(boot(:,i,j));
        plot(xi,f,'Color',col(i,:),'LineWidth',2)
    end
    if j == 2
        ylabel('Probability')
    end
    ylim([0 0.3])
end
xlabel('Standard deviation of directional error (degrees)')
legend(blocks)

boot = sort(boot,1);
confInterval = boot([26 975],:,:);
bootMu = mean(boot,1);
confDiff = abs(confInterval - repmat(bootMu,[2 1 1]));
bootMu = squeeze(bootMu);

figure(7); clf; hold on
for i = 1:3
    errorbar(i,bootMu(1,i),confDiff(1,1,i),confDiff(2,1,i),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'LineWidth',2)
    errorbar(i+4,bootMu(2,i),confDiff(1,2,i),confDiff(2,2,i),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'LineWidth',2,'HandleVisibility','off')
    errorbar(i+8,bootMu(3,i),confDiff(1,3,i),confDiff(2,3,i),'-o','Color',col(i,:),'MarkerFaceColor',col(i,:),'LineWidth',2,'HandleVisibility','off')
end
xticks(2:4:10)
xticklabels(blocks)
xlim([0.5 11.5])
ylabel('Standard deviation of directional error (degrees)')
legend(graph_names)

%% plot kernel density estimates
figure(8); clf
for i = 1:Ngroup
    for j = 1:3
        trial = trials{gblocks(i,j)};
        subplot(1,3,j); hold on
        [f,xi] = ksdensity(reshape(dirError{i}(trial,:),[numel(dirError{i}(trial,:)) 1]));
        plot(xi,f,'LineWidth',2)
        if i == 3
            title(blocks{j})
            axis([-180 180 0 .06])
            xticks(-180:90:180)
            yticks(0:0.01:0.06)
            box off
            if j == 1
                ylabel('Probability')
            end
        end
    end
end

%%
Ntrials2 = 100;
trialsAll = {1:100,101:200,201:300};
for k = 1:3
    trials = trialsAll{k};
    for i = 1:Ngroup
        Nsubj = length(d.(groups{i}));
        for j = 1:Nsubj
            ang = atan2(d.(groups{i}){j}.targetRel(trials,2), d.(groups{i}){j}.targetRel(trials,1));
            angMir = atan2(d.(groups{i}){j}.targetRel(trials,2), -d.(groups{i}){j}.targetRel(trials,1));

            dir = d.(groups{i}){j}.initDir_noRot(trials)';

            error{i}(:,j) = dir-ang;
            errorMir{i}(:,j) = dir-angMir;
        end
        
        for l = 1:numel(error{i})
            while error{i}(l) >= pi
                error{i}(l) = error{i}(l)-2*pi;
            end
            while error{i}(l) < -pi
                error{i}(l) = error{i}(l)+2*pi;
            end
            while errorMir{i}(l) >= pi
                errorMir{i}(l) = errorMir{i}(l)-2*pi;
            end
            while errorMir{i}(l) < -pi
                errorMir{i}(l) = errorMir{i}(l)+2*pi;
            end
        end
        
        towardMir{i}.all(:,k) = sum(abs(error{i})>abs(errorMir{i}),1);
        towardMir{i}.mean(k) = mean(towardMir{i}.all(:,k),1);
    end
end

figure(9); clf; hold on
for i = 1:3
    plot(i,towardMir{1}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i,towardMir{1}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
    
    plot(i+4,towardMir{2}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i+4,towardMir{2}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
    
    plot(i+8,towardMir{3}.all(:,i),'.','Color',col(i,:),'MarkerSize',15)
    plot(i+8,towardMir{3}.mean(i),'.','Color',col(i,:),'MarkerSize',40)
end
xticks([2 6 10])
xticklabels(graph_names)
ylabel('Percent towards mirrored target')
axis([0.5 11.5 0 70])
box off

%%
edges = -180:10:180;
for j = 1:3
    figure(6+j-1); clf
    for i = 1:4
        bin = find(bins{j}==i);
        early = bin(bin<=100);
        post = bin(bin>200);
        a = (bin<=100)+(bin>200);
        a = ~a;
        late = bin(a);
        
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