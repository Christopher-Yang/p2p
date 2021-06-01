
rng(2);
vmPDF = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa)));
unifPDF = 1 / 2*pi;

weight1 = 0:0.1:1;
weight2 = 0:0.1:1;
kappa = 1:10;

mu1 = pi/4;
mu2 = atan2(sin(mu1), -cos(mu1));
Nreps = 1;

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
            w3 = 1 - (w1 + w2);
            
            if w1 ~= 0
                data1 = vmrand(mu1,kap,[round(w1*100) 1]);
            else
                data1 = [];
            end
            
            if w2 ~= 0
                data2 = vmrand(mu2,kap,[round(w2*100) 1]);
            else
                data2 = [];
            end
            
            data_sim = [data1; data2];
            data3 = 2*pi*(rand(100-length(data_sim), 1) - 0.5);
            data_sim = [data_sim; data3];
            
            log_likelihood = @(params) calc_likelihood(params, data_sim, mu1, mu2);            
            [params_opt, neg_log_likelihood] = fmincon(log_likelihood, paramsInit, A, b);
            
            weight1_opt(i,j,k) = params_opt(1);
            weight2_opt(i,j,k) = params_opt(2);
            kappa_opt(i,j,k) = params_opt(3);
            likelihood(i,j,k) = -neg_log_likelihood;
        end
    end

    idx = ~isnan(weight1_opt(:,:,k));
    
    weight1_vec = repmat(weight1',[1 length(weight2)]);
    weight2_vec = repmat(weight2,[length(weight1) 1]);
    weight1_vec = weight1_vec(idx);
    weight2_vec = weight2_vec(idx);
    
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

%%
k = 3;

figure(1); clf
subplot(1,3,1); hold on
plot([0 1], [0 1], 'k')
plot(weight1, weight1_opt(:,:,k), '.k', 'MarkerSize', 10)
plot(weight1, nanmean(weight1_opt(:,:,k)'), '.r', 'MarkerSize', 30)
title(['\kappa = ' num2str(kappa(k))])
axis([0 1 0 1])
axis square
xlabel('weight1 (actual)')
ylabel('weight1 (fit)')

subplot(1,3,2); hold on
plot([0 1], [0 1], 'k')
plot(weight2, weight2_opt(:,:,k)', '.k', 'MarkerSize', 10)
plot(weight2, nanmean(weight2_opt(:,:,k)), '.r', 'MarkerSize', 30)
title(['\kappa = ' num2str(kappa(k))])
axis([0 1 0 1])
axis square
xlabel('weight2 (actual)')
ylabel('weight2 (fit)')

subplot(1,3,3); hold on
plot([0 max(kappa)], [0 max(kappa)], 'k')
for i = 1:length(kappa)
    dat = kappa_opt(:,:,i);
    plot(repmat(kappa(i), size(dat)), dat, '.k', 'MarkerSize', 10)
    plot(kappa(i), nanmean(dat(:)), '.r', 'MarkerSize', 30)
end
ylim([0 30])
axis square
xlabel('kappa (actual)')
ylabel('kappa (fit)')

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
