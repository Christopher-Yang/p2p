clear all
load BimanualSkillData
load Bimanual_compact

pkVel = struct
pkVel.Bi.large = []
pkVel.Bi.small = []
pkVel.L.large = []
pkVel.L.small = []
pkVel.R.large = []
pkVel.R.small = []
nSubj = 8

for subj = 1:nSubj
    jumps = [d_full.Bi{subj}{1}.pkVel, d_full.Bi{subj}{5}.pkVel]
    pkVel.Bi.large = [pkVel.Bi.large, mean(jumps)]
end

for subj = 1:nSubj
    jumps = [d_full.L{subj}{1}.pkVel, d_full.L{subj}{5}.pkVel]
    pkVel.L.large = [pkVel.L.large, mean(jumps)]
end

for subj = 1:nSubj
    jumps = [d_full.R{subj}{1}.pkVel, d_full.R{subj}{5}.pkVel]
    pkVel.R.large = [pkVel.R.large, mean(jumps)]
end

for subj = 1:nSubj
    jumps = [d_full.Bi{subj}{2}.pkVel, d_full.Bi{subj}{4}.pkVel]
    pkVel.Bi.small = [pkVel.Bi.small, mean(jumps)]
end

for subj = 1:nSubj
    jumps = [d_full.L{subj}{2}.pkVel, d_full.L{subj}{4}.pkVel]
    pkVel.L.small = [pkVel.L.small, mean(jumps)]
end

for subj = 1:nSubj
    jumps = [d_full.R{subj}{2}.pkVel, d_full.R{subj}{4}.pkVel]
    pkVel.R.small = [pkVel.R.small, mean(jumps)]
end

figure
subplot(1,2,1)
scatter(pkVel.R.large, d.R.peakResponse_large, 'r', 'filled')
hold on
scatter(pkVel.L.large, d.L.peakResponse_large, 'g', 'filled')
scatter(pkVel.Bi.large, d.Bi.peakResponse_large, 'b', 'filled')
ntitle('Large Jumps')
xlabel('Primary Velocity')
ylabel('Secondary Velocity')
%axis([0.002 0.003 0 0.00002])
legend('Right', 'Left', 'Bimanual')

subplot(1,2,2)
scatter(pkVel.R.small, d.R.peakResponse_small, 'r', 'filled')
hold on
scatter(pkVel.L.small, d.L.peakResponse_small, 'g', 'filled')
scatter(pkVel.Bi.small, d.Bi.peakResponse_small, 'b', 'filled')
ntitle('Small Jumps')
xlabel('Primary Velocity')
ylabel('Secondary Velocity')
%axis([0.002 0.003 0 0.00002])
legend('Right', 'Left', 'Bimanual')