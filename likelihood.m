function likelihood()

rng(2);

% generate targets
Ntrials = 100;
targ_gd = 2 * pi * (rand(Ntrials,1)-0.5);
targ_dir = [cos(targ_gd) sin(targ_gd)];
targ_dir(:,1) = -targ_dir(:,1);
targ_hab = atan2(targ_dir(:,2), targ_dir(:,1));

% generate data from mixture model
mu1 = pi/4;
mu2 = atan2(sin(mu1), -cos(mu1));
w1 = 0.5;
w2 = 0.2;
kappa_real = 3;

weight1 = 0:0.05:1;
weight2 = 0:0.05:1;
kappa = 1:10;

Nsims = 13;
test = NaN(Nsims, 3);
for m = 1:Nsims
    
    % generate data from mixture model
    idx1 = round(w1*Ntrials);
    idx2 = idx1 + round(w2*Ntrials);

    data1 = vmrand(targ_gd(1:idx1), kappa_real);
    data2 = vmrand(targ_hab(idx1+1:idx2), kappa_real);
    data_mix = [data1; data2];
    data3 = 2*pi*(rand(Ntrials-length(data_mix), 1) - 0.5);
    data_mix = [data_mix; data3];

    like = NaN(length(weight1), length(weight2), length(kappa));
    for k = 1:length(kappa)
        kap = kappa(k);
        
        for i = 1:length(weight1)
            for j = 1:length(weight2)-i+1
                like(i,j,k) = calc_likelihood(weight1(i), weight2(j), kap, data_mix, targ_gd, targ_hab);
            end
        end
    end
    
    % k = 25;
    % like_max = max(max(like(:,:,k)));
    %
    % y = max(like(:,:,k));
    % [~,column] = max(y);
    % [~,row] = max(like(:,column,k));
    %
    % [x, y] = meshgrid(weight1,weight2);
    % figure(1); clf
    % surf(y, x, like(:,:,k), 'HandleVisibility', 'off'); hold on
    % xlabel('weight1')
    % ylabel('weight2')
    % zlabel('log likelihood')
    % % zlim([-250 -150])
    % scatter3(w1,w2,like_max,30,[0 0 0],'filled')
    % scatter3(weight1(row),weight2(column),like_max,30,[1 0 0],'filled')
    % legend({'True params','Max likelihood'})
    
    like_vec = like(:);
    for i = 1:numel(like_vec)
        if like_vec(i) < -200
            like_vec(i) = NaN;
        end
    end
    
    [~, idx] = max(like_vec);
    
    x = repmat(weight1',[1 length(weight2) length(kappa)]);
    x = x(:);
    y = repmat(weight2,[length(weight1) 1 length(kappa)]);
    y = y(:);
    z = repmat(permute(kappa, [1 3 2]),[length(weight1) length(weight2) 1]);
    z = z(:);
    
    test(m,:) = [x(idx) y(idx) z(idx)];
end


figure(2); clf
scatter3(x, y, z, 15, like_vec); hold on
scatter3(x(idx), y(idx), z(idx), 100, 'k', 'filled')
colormap(parula);
colorbar
title(['Actual: w_1 = ' num2str(w1) '; w_2 = ' num2str(w2) '; \kappa = ' ... 
    num2str(kappa_real) newline 'Fit: w_1 = ' num2str(x(idx)) '; w_2 = ' ...
    num2str(y(idx)) '; \kappa = ' num2str(z(idx))])
xlabel('weight (GD)')
ylabel('weight (hab)')
zlabel('kappa')
end

function log_likelihood = calc_likelihood(weight1, weight2, kappa, samples, target_gd, target_hab)
    pdf = @(x, mu, kappa) (exp(kappa*cos(x-mu)) / (2 * pi * besseli(0,kappa))); % PDF of von Mises distribution
    
    weightUnif = 1 - (weight1 + weight2);
    
    likelihood_vm1 = weight1 * pdf(samples, target_gd, kappa);
    likelihood_vm2 = weight2 * pdf(samples, target_hab, kappa);
    likelihood_unif = repelem(weightUnif/(2*pi),length(samples))';
    
    likelihood_all = sum([likelihood_vm1 likelihood_vm2 likelihood_unif],2);
    log_likelihood = sum(log(likelihood_all));
end