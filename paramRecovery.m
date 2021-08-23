function paramRecovery()

rng(2);

weight1 = 0:0.1:1;
weight2 = 0:0.1:1;
kappa = 2:2:10;
R = besseli(1,kappa) ./ besseli(0,kappa);
sd = sqrt(-2 * log(R)) * 180/pi;

% generate targets
Ntrials = 100;

clear corr_weight1 corr_weight2
weight1_opt = NaN(length(weight1), length(weight2), length(kappa));
weight2_opt = NaN(length(weight1), length(weight2), length(kappa));
kappa_opt = NaN(length(weight1), length(weight2), length(kappa));
likelihood = NaN(length(weight1), length(weight2), length(kappa));

paramsInit = [0.33 0.33 10];
A = [-1 0 0; 0 -1 0; 1 1 0; 0 0 1; 0 0 -1];
b = [0 0 1 200 0]';

for k = 1:length(kappa)
    kap = kappa(k);
    
    for i = 1:length(weight1)
        for j = 1:length(weight2)-i+1
            w1 = weight1(i);
            w2 = weight2(j);
            
            targ_gd = 2 * pi * (rand(Ntrials,1)-0.5);
            targ_dir = [cos(targ_gd) sin(targ_gd)];
            targ_dir(:,1) = -targ_dir(:,1);
            targ_hab = atan2(targ_dir(:,2), targ_dir(:,1));
            
            % generate data from mixture model
            idx1 = round(w1*100);
            idx2 = round(w2*100);
            
            if idx1 == 0
                data1 = [];
            else
                data1 = vmrand(targ_gd(1:idx1), kap);
            end
            
            if idx2 == 0
                data2 = [];
            else
                data2 = vmrand(targ_hab(idx1+1:idx1+idx2), kap);
            end
            
            data_sim = [data1; data2];
            data3 = 2*pi*(rand(Ntrials-length(data_sim), 1) - 0.5);
            data_sim = [data_sim; data3];
            
            log_likelihood = @(params) calc_likelihood(params, data_sim, targ_gd, targ_hab);            
            [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
            
            weight1_opt(i,j,k) = params_opt(1);
            weight2_opt(i,j,k) = params_opt(2);
            kappa_opt(i,j,k) = params_opt(3);
            likelihood(i,j,k) = -neg_log_likelihood;
        end
    end

    idx = ~isnan(weight1_opt(:,:,k));
    
    weight1_mat = repmat(weight1',[1 length(weight2)]);
    weight2_mat = repmat(weight2,[length(weight1) 1]);
    weight1_vec = weight1_mat(idx);
    weight2_vec = weight2_mat(idx);
    
    weight1_opt_vec = weight1_opt(:,:,k);
    weight2_opt_vec = weight2_opt(:,:,k);
    weight1_opt_vec = weight1_opt_vec(idx);
    weight2_opt_vec = weight2_opt_vec(idx);
    
    c1 = corrcoef(weight1_vec,weight1_opt_vec);
    c2 = corrcoef(weight2_vec,weight2_opt_vec);
    
    corr_weight1(k) = c1(1,2);
    corr_weight2(k) = c2(1,2);
end

Ndat = sum(sum(~isnan(weight1_opt(:,:,1))));

kappa_vec = repmat(kappa, [Ndat 1]);
kappa_vec = kappa_vec(:);
kappa_opt_vec = kappa_opt(~isnan(kappa_opt));

c3 = corrcoef(kappa_vec,kappa_opt_vec);
corr_kappa = c3(1,2);

figure(2); clf
subplot(1,2,1); hold on
plot(kappa, corr_weight1)
ylim([0 1])

subplot(1,2,2); hold on
plot(kappa, corr_weight2)
ylim([0 1])

%%
k = 2;

R = besseli(1,kappa_opt) ./ besseli(0,kappa_opt);
sd_opt = sqrt(-2 * log(R)) * 180/pi;

weights = weight1_mat + weight2_mat;
idx = (weights <= 1) == (weights >= 0.8);

figure(1); clf
subplot(1,3,1); hold on
plot([0 1], [0 1], 'k')
for i = 1:length(weight1)
    plot(weight1(i), weight1_opt(i,:,k), '.k', 'MarkerSize', 10)
    plot(weight1(i), mean(weight1_opt(i,:,k),'omitnan'), '.r', 'MarkerSize', 30)
end
title(['\kappa = ' num2str(kappa(k))])
axis([0 1 0 1])
axis square
xlabel('weight1 (actual)')
ylabel('weight1 (fit)')

subplot(1,3,2); hold on
plot([0 1], [0 1], 'k')
for i = 1:length(weight2)
    plot(weight2(i), weight2_opt(:,i,k), '.k', 'MarkerSize', 10)
    plot(weight2(i), mean(weight2_opt(:,i,k),'omitnan'), '.r', 'MarkerSize', 30)
end
title(['\kappa = ' num2str(kappa(k))])
axis([0 1 0 1])
axis square
xlabel('weight2 (actual)')
ylabel('weight2 (fit)')

subplot(1,3,3); hold on
plot([0 max(kappa)], [0 max(kappa)], 'k')
for i = 1:length(kappa)
    dat = kappa_opt(:,:,i);
    dat = dat(idx);
    plot(kappa(i), dat(:), '.k', 'MarkerSize', 10)
    plot(kappa(i), nanmean(dat(:)), '.r', 'MarkerSize', 30)
end
% ylim([0 30])
axis square
xlabel('sd (actual)')
ylabel('sd (fit)')

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
