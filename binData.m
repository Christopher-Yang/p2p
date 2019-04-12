function [binnedData] = binData(data,binSize)
% bins data into bins of given size, then averages
%
% Input:    data    - Nsubj x Ntrials matrix
%           binSize - Number of trials per bin
%
% Output:   binnedData - Nsubj x (Ntrials/binSize)
%           binSTD     - Nsubj x (Ntrials/binSize)
%
% NB - Nsubj must be an integer multiple of binSize
[Nsubj Ntrials] = size(data);

for k=1:Nsubj
    binnedData(k,:) = nanmean(reshape(data(k,:),binSize,Ntrials/binSize));
    %binSTD(k,:) = stdNaN(reshape(data(k,:),binSize,Ntrials/binSize));
    %binSE(k,:) = seNaN(reshape(data(k,:),binSize,Ntrials/binSize));
end
