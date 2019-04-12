clear all
load Bimanual_compact

handedness = [1, 0.9, 0.1, 0.6, 1, 1, 0.4, 0.9]

%diffLargeVel = d.R.peakResponse_large - d.L.peakResponse_large
%diffSmallVel = d.R.peakResponse_small - d.L.peakResponse_small
%figure('name','Velocity Difference vs. Handedness')
%scatter(handedness, diffLargeVel, 50, 'r', 'filled'); hold on; scatter(handedness, diffSmallVel, 50, 'r')
%xlabel('Handedness Score')
%ylabel('Right Hand Velocity - Left Hand Velocity')
%legend('Large Velocity', 'Small Velocity')

figure
subplot(1,2,1)
scatter(handedness, d.R.peakResponse_small, 'r', 'filled')
hold on
scatter(handedness, d.L.peakResponse_small, 'g', 'filled')
scatter(handedness, d.Bi.peakResponse_small, 'b', 'filled')
ntitle('Small Jump')
xlabel('Handedness Score')
ylabel('Velocity')
axis([0 1 0 0.00002])
legend('Right', 'Left')

subplot(1,2,2)
scatter(handedness, d.R.peakResponse_large, 'r', 'filled')
hold on
scatter(handedness, d.L.peakResponse_large, 'g', 'filled')
scatter(handedness, d.Bi.peakResponse_large, 'b', 'filled')
ntitle('Large Jump')
xlabel('Handedness Score')
ylabel('Velocity')
axis([0 1 0 0.00002])
legend('Right', 'Left')
%saveas(gcf,'Velocity_vs_Handedness.png')

[resultSmall,pSmall] = ttest(d.R.peakResponse_small,d.L.peakResponse_small)
[resultLarge,pLarge] = ttest(d.R.peakResponse_large,d.L.peakResponse_large)