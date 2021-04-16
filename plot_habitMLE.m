% rotate reach directions into quadrant 1 and plot histograms

rng(2);
load variables/vmFit
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
gblocks{1} = [1 2 5 6];
gblocks{2} = [1 2 14 15];
gblocks{3} = [1 2 29 30];

% set colors for generating plots
col = lines;
col = col(1:7,:);
col2 = [180 180 0
        0 191 255
        255 99 71]./255;
    
R(:,:,1) = [-1 0; 0 1];
R(:,:,2) = [-1 0; 0 -1];
R(:,:,3) = [1 0; 0 -1];

clear target_q1 target_q2 dir_rot targets

% store all reach direction errors into dirError
for i = 1:Ngroup
    dir_rot{i} = NaN(100, allSubj(i), 4);
    targets{i} = NaN(100, 2, 4);
    target_q1{i} = NaN(100, 4);
    target_q2{i} = NaN(100, 4);

    for k = 1:4
        Ntrials = length(d.(groups{i}){1}.Cr);
        Nsubj = length(d.(groups{i}));
        dir = NaN(Ntrials,Nsubj);
        for j = 1:Nsubj
            dir(:,j) = d.(groups{i}){j}.initDir_noRot;
        end
        
%         if k == 1
%             dir_rot{i}(:,:,k) = dir(end-199:end-100,:);
%             targets{i}(:,:,k) = d.(groups{i}){1}.targetRel(end-199:end-100,:);
%         else
%             dir_rot{i}(:,:,k) = dir(end-99:end,:);
%             targets{i}(:,:,k) = d.(groups{i}){1}.targetRel(end-99:end,:);
%         end
        
        dir2 = dir(trials{gblocks{i}(k)},:);
        targ = d.(groups{i}){1}.targetRel(trials{gblocks{i}(k)},:);

        quad{1} = (targ(:,1) < 0) + (targ(:,2) > 0) == 2; % quadrant 2
        quad{2} = (targ(:,1) < 0) + (targ(:,2) < 0) == 2; % quadrant 3
        quad{3} = (targ(:,1) > 0) + (targ(:,2) < 0) == 2; % quadrant 4
        
        for m = 1:3
            for j = 1:allSubj(i)
                h = dir2(quad{m},j);
                p = R(:,:,m)*[cos(h) sin(h)]';
                dir2(quad{m},j) = atan2(p(2,:),p(1,:));
            end
        end
        
        t = abs(targ);
        if k == 1
            idx = 1:30;
        else
            idx = 1:100;
        end
        
        target_q1{i}(idx,k) = atan2(t(:,2), t(:,1));
        target_q2{i}(idx,k) = atan2(t(:,2), -t(:,1));
        dir_rot{i}(idx,:,k) = dir2;
        targets{i}(idx,:,k) = targ;
    end
end

%%
figure(1); clf
for i = 1:Ngroup
    for k = 1:4
        dir2 = dir_rot{i}(:,:,k);
        idx = ~isnan(dir2);
        
        subplot(4,3,3*(k-1)+i)
        polarhistogram(dir2(idx),-pi:pi/8:pi,'Normalization','pdf')
        rlim([0 0.85])
        if k == 1
            title(graph_names{i})
        end
    end
end

%%
rng(24);
xAxis = [0.5 2.5 5 15
         1 3 8 15.5
         1.5 3.5 13 16];
     
figure(2); clf; hold on
for i = 1:3
    q2_num = sum(cos(dir_rot{i}) < 0, 1);
    q2_den = sum(~isnan(dir_rot{i}));
    q2_frac = squeeze(q2_num ./ q2_den);
    
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5),q2_frac, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:) ,mean(q2_frac), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
xlim([0 16.5])
ylabel('Fraction of reaches in Q2')
set(gca,'TickDir','out')
legend(graph_names, 'Location', 'northwest')

%% perform MLE

clear weight1_opt weight2_opt weight3_opt
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution

block = 3;
group = 3;
subj = 2;
tolerance = 0.01; % tolerance for stopping EM loop

weight1 = 0.33;
weight2 = 0.33;

for k = 1:4
    for j = 1:Ngroup
        for i = 1:allSubj(j)
            kappa = vmFit.kappa_opt(block,j,i);
            idx = ~isnan(dir_rot{j}(:,i,k));
            
            log_likelihood = @(params) calc_likelihood(params, dir_rot{j}(idx,i,k), target_q1{j}(idx,k), target_q2{j}(idx,k), kappa);
            
            weightsInit = [weight1 weight2];
            A = [-1 0; 0 -1; 1 1];
            b = [0 0 1]';
            [params_opt, fval] = fmincon(log_likelihood, weightsInit, A, b, [], [], [0, 0], [1, 1]);
            
            weight1_opt{j}(i,k) = params_opt(1);
            weight2_opt{j}(i,k) = params_opt(2);
            weight3_opt{j}(i,k) = 1 - sum(params_opt);
        end
    end
end

%% plot P(habitual)
figure(3); clf
subplot(1,2,1); hold on
for i = 1:3
    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), weight2_opt{i}, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(weight2_opt{i}), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 0.8])
ylabel('P(habitual)')
set(gca,'TickDir','out')
legend(graph_names)

subplot(1,2,2); hold on
for i = 1:3
    weight2_frac = weight2_opt{i} ./ (weight1_opt{i} + weight2_opt{i});

    n = allSubj(i);
    plot(repmat(xAxis(i,:),[n 1]) + 0.5*(rand(n,4) - 0.5), weight2_frac, '.', 'MarkerSize', 20, 'Color', col2(i,:), 'HandleVisibility', 'off')
    plot(xAxis(i,:), mean(weight2_frac), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', col2(i,:), 'LineWidth', 1)
end
xticks([1 3 5 8 13 15.5])
xticklabels({'Baseline',1,2,5,10,'Flip'})
xlabel('Day')
axis([0 16.5 0 0.8])
ylabel('P(habitual) / [P(habitual) + P(GD)]')
set(gca,'TickDir','out')

function neg_log_likelihood = calc_likelihood(weights, samples, target_q1, target_q2, kappa)
    vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
    
    weightUnif = 1 - sum(weights);
    likelihood_vm1 = weights(1) * vmPDF(samples, target_q1, kappa);
    likelihood_vm2 = weights(2) * vmPDF(samples, target_q2, kappa);
    likelihood_unif = repelem(weightUnif/(2*pi),length(samples))';
    
    likelihood_all = sum([likelihood_vm1 likelihood_vm2 likelihood_unif],2);
    neg_log_likelihood = -sum(log(likelihood_all));
end