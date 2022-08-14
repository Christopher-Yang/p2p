% performs model recovery/comparison analysis

function modelRecovery(loadAccuracy)

% range of parameters to simulate von Mises model
w1 = 0:0.1:1; % weight of goal-directed von Mises
w2 = 0:0.1:1; % weight of habitual von Mises
kappa = 3; % concentration of von Mises

% load matrix with precomputed accuracy
if loadAccuracy
    load Variables/accuracy

% compute accuracy by model fitting
else
    
    % set variables for analysis
    rng(2);
    paramsInit = [0.33 0.33 10]; % initial parameters for von Mises model
    weightInit = 0.8; % intial weight for uniform model
    Ntrials = 100; % number of trials to simulate in one run
    Nsims = 50; % number of times simulate 
    accuracy = NaN(length(w1),length(w2));
    
    % constraints for optimization
    A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
    b = [0 0 1 200 0]';
        
    for m = 1:length(w1) % loop over all weights of goal-directed distribution
        for k = 1:length(w2)-m+1 % loop over weights of habitual distribution, constrained by choice of w1
            
            % generate targets
            targ_gd = 2 * pi * (rand(100,1)-0.5);
            targ_dir = [cos(targ_gd) sin(targ_gd)];
            targ_dir(:,1) = -targ_dir(:,1);
            targ_hab = atan2(targ_dir(:,2), targ_dir(:,1));
            
            % generate data from uniform model
            negIdx = targ_gd < 0;
            data_unif = pi * rand(100,1);
            data_unif(negIdx) = -data_unif(negIdx);
            
            % generate data from von Mises model
            idx1 = round(w1(m)*100);
            idx2 = round(w2(k)*100);
            
            % simulate goal-directed reaches
            if idx1 == 0
                data1 = [];
            else
                data1 = vmrand(targ_gd(1:idx1), kappa);
            end
            
            % simulate habitual reaches
            if idx2 == 0
                data2 = [];
            else
                data2 = vmrand(targ_hab(idx1+1:idx1+idx2), kappa);
            end
            
            % simulate random reaches and store data in data_mix
            data_mix = [data1; data2];
            data3 = 2*pi*(rand(100-length(data_mix), 1) - 0.5);
            data_mix = [data_mix; data3];
            
            % preallocate confusion matrix
            confusion = zeros(2);
            
            % run simulation Nsims times
            for i = 1:Nsims
                
                % fit mix model to mix data
                mix_likelihood = @(params) calc_likelihood(params, data_mix, targ_gd, targ_hab);
                [~, neg_log_likelihood] = fmincon(mix_likelihood, paramsInit, A, b);
                log_likelihood = -neg_log_likelihood;
                BIC_mix = length(paramsInit) * log(Ntrials) - 2 * log_likelihood;
                
                % fit unif model to mix data
                unif_likelihood = @(params) calc_likelihood_unif(params, data_mix, targ_gd);
                [~, neg_log_likelihood] = fmincon(unif_likelihood, weightInit, [], [], [], [], 0, 1);
                log_likelihood = -neg_log_likelihood;
                BIC_unif = length(weightInit) * log(Ntrials) - 2 * log_likelihood;
                
                % fill in first row of confusion matrix
                if BIC_mix < BIC_unif
                    confusion(1,1) = confusion(1,1) + 1;
                else
                    confusion(1,2) = confusion(1,2) + 1;
                end
                
                % fit mix model to unif data
                mix_likelihood = @(params) calc_likelihood(params, data_unif, targ_gd, targ_hab);
                [~, neg_log_likelihood] = fmincon(mix_likelihood, paramsInit, A, b);
                log_likelihood = -neg_log_likelihood;
                BIC_mix = length(paramsInit) * log(Ntrials) - 2 * log_likelihood;
                
                % fit unif model to unif data
                unif_likelihood = @(params) calc_likelihood_unif(params, data_unif, targ_gd);
                [~, neg_log_likelihood] = fmincon(unif_likelihood, weightInit, [], [], [], [], 0, 1);
                log_likelihood = -neg_log_likelihood;
                BIC_unif = length(weightInit) * log(Ntrials) - 2 * log_likelihood;
                
                % fill in second row of confusion matrix
                if BIC_unif < BIC_mix
                    confusion(2,2) = confusion(2,2) + 1;
                else
                    confusion(2,1) = confusion(2,1) + 1;
                end
            end
            
            accuracy(m,k) = sum(diag(confusion))/sum(confusion,'all');
        end
    end
    
    save variables/accuracy accuracy
end

% store data in more easily plottable forms
x = repmat(w1',[1 length(w2) length(kappa)]);
x = x(:);
y = repmat(w2,[length(w1) 1 length(kappa)]);
y = y(:);

%% Supplementary Figure 1C
f = figure(16); clf
set(f,'Position',[200 200 200 150]);
scatter(x,y,13,accuracy(:),'filled')
colormap(copper)
c = colorbar;
c.Ticks = [0.5 1];
xlabel('Goal-directed component weight')
xticks(0:0.5:1)
ylabel('Habitual component weight')
yticks(0:0.5:1)
axis([0 1 0 1])
axis square
set(gca,'TickDir','out')

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/model_recovery','-dpdf','-painters')

end

% functions for computing likelihood
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