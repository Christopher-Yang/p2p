function response = getResponses(x,pmax,dt)
% returns mean responses as:
%   response.pos
%   response.vel
%   response.acc
%   response.time

x(isnan(x)) = pmax; % eliminate nans
response.pos = mean(x)';
%dAll.posResponse_large(subj,:,c) = nanmean([-d{subj}.Bi{c}{1}.CrX_post ; d{subj}.Bi{c}{5}.CrX_post]);
response.time = dt*([1:size(x,2)]-1)';
response.vel = diff(savgolayFilt(response.pos',3,9))'/dt;
response.acc = diff(savgolayFilt(response.vel',3,9))'/dt;