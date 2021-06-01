rng(2);
A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
b = [0 0 1 200 0]';
paramsInit = [0.33 0.33 10];
weightInit = 0.8;

% generate targets
targ_gd = 2 * pi * (rand(100,1)-0.5);
targ_dir = [cos(targ_gd) sin(targ_gd)];
targ_dir(:,1) = -targ_dir(:,1);
targ_hab = atan2(targ_dir(:,2), targ_dir(:,1));

% generate data from uniform model
negIdx = targ_gd < 0;
data_unif = pi * rand(100,1);
data_unif(negIdx) = -data_unif(negIdx);

% generate data from mixture model
mu1 = pi/4;
mu2 = atan2(sin(mu1), -cos(mu1));
w1 = 0.25;
w2 = 0.25;
kappa = 3;

idx1 = round(w1*100);
idx2 = idx1 + round(w2*100);

data1 = vmrand(targ_gd(1:idx1), kappa);
data2 = vmrand(targ_hab(idx1+1:idx2), kappa);
data_mix = [data1; data2];
data3 = 2*pi*(rand(100-length(data_mix), 1) - 0.5);
data_mix = [data_mix; data3];

% fitting
Ntrials = 100;
Nsims = 100;
confusion = zeros(2);

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
    
    if BIC_unif < BIC_mix
        confusion(2,2) = confusion(2,2) + 1;
    else
        confusion(2,1) = confusion(2,1) + 1;
    end
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