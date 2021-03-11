function output = loadSubjData(path)
% Load a single subject's point-to-point data for a given block 
% (blocknames) and starting position of the target (START). The output 
% variable "data" contains cell arrays to be used for analysis in
% processData(). See that function for details on the meaning of each
% array.

disp('Loading...');

% fields from data file to be analyzed
fields = {'time','cursorX','cursorY','begin'};
fields2 = {'extraTime','extraCursorX','extraCursorY','extraBegin'};

allFields = [fields fields2];

fnames = dir(path);
Nsubj = length(fnames)-2;

for i = 1:Nsubj
    disp(['   Subject ' num2str(i)]); % display progress
    
    % read data
    d = readtable([path fnames(2+i).name]);
    
    % only analyze tracking trials
    p2p = strcmp(d.task,'p2p');
    d = d(p2p,:);
    
    % extract useful variables from the data table
    Ntrials = size(d,1);
    frameRate = d.frameRate(1);
    
    % set variables used in the analysis
    frameDrops = NaN(1,Ntrials);
    tolerance = (1/frameRate) + 0.5/frameRate; % tolerance for defining duration of one frame
    count = 1; % index for storing sine parameters in cell array
    
    % organize data into appropriate format
    for k = 1:Ntrials
        trial = [];
        extraData = [];
        
        % format data from strings into matrices
        for j = 1:length(allFields)
            data = d.(allFields{j}){k}(2:end-1); % take of apostrophes and brackets from string
            data = strsplit(data,','); % separate data by commas
            data = cellfun(@str2double,data); % convert data into a matrix
            
            if sum(strcmp(allFields{j},fields))
                trial = [trial data'];
            elseif sum(strcmp(allFields{j},fields2))
                extraData = [extraData data(~isnan(data))'];
            end
        end
        
        Nsamples = size(trial,1);
        
        % insert missing data (extraData) into the the main data
        % (trial)
        for m = 1:size(extraData,1)
            
            % find the time point immediately after and before missing
            % data
            after = find(trial(:,1)>=extraData(m,1),1);
            before = find(trial(:,1)<extraData(m,1),1,'last');
            
            replace = 1; % used to end while loop
            check = 1; % defines search breadth for NaNs
            while replace
                % if NaN is first data point, place missing data there
                if after == 1
                    trial(1,:) = extraData(m,:);
                    replace = 0;
                    
                    % if NaN is last data point, place missing data there
                elseif before == size(trial,1)
                    trial(end,:) = extraData(m,:);
                    replace = 0;
                    
                    % if NaN is between before and after, place the missing
                    % data there
                elseif after-before == 2
                    trial(before+1,:) = extraData(m,:);
                    replace = 0;
                    
                    % if there's 2 NaNs between before and after, fill
                    % missing data in index that's closer in time to before
                    % or after
                elseif after-before == 3
                    if trial(after,1)-extraData(m,1) > extraData(m,1)-trial(before,1)
                        trial(before+1,:) = extraData(m,:);
                    else
                        trial(after-1,:) = extraData(m,:);
                    end
                    replace = 0;
                    
                    % if more than 3 NaNs between before and after, don't
                    % fill data
                elseif after-before > 3
                    disp('   Multiple NaNs between before and after')
                    replace = 0;
                    
                    % if NaN is outside of before and after, place the
                    % missing data as specified below
                else
                    % find the indices of NaN(s)
                    if before-check < 1
                        first = 1;
                    else
                        first = before-check;
                    end
                    if after+check > Nsamples
                        last = Nsamples;
                    else
                        last = after+check;
                    end
                    idxRange = first:last;
                    
                    try
                        missing = idxRange(isnan(trial(idxRange,1)));
                    catch
                        error('Something went wrong')
                    end
                    
                    % do different things based on number of NaNs found
                    switch length(missing)
                        % NaN not found
                        case 0
                            if check > 18 % only search within 300 ms of desired data position
                                disp('NaN could not be found');
                                replace = 0;
                            else % increase search range
                                check = check + 1;
                            end
                            
                            % one NaN found
                        case 1
                            % determine whether Nan is earlier or later
                            % in time
                            if missing < before
                                idx = missing:before;
                            elseif missing > after
                                idx = after:missing;
                            end
                            
                            % create "chunk," which will replace a
                            % subsection of the matrix
                            chunk = trial(idx,:);
                            chunk = sortrows([chunk; extraData(m,:)]);
                            chunk(end,:) = [];
                            trial(idx,:) = chunk;
                            
                            % end the while loop
                            replace = 0;
                            
                            % two NaNs found
                        case 2
                            chunk = sortrows([trial(missing(1)+1:missing(2)-1,:); extraData(m,:)]);
                            
                            % check whether to fill data by pushing
                            % data forward or backward
                            if ~isnan(trial(missing(1)-1,1)) % push data backward
                                if chunk(1,1) - trial(missing(1)-1,1) < tolerance
                                    trial(missing(1):missing(2)-1,:) = chunk;
                                end
                            elseif ~isnan(trial(missing(2)+1,1)) % push data forward
                                if trial(missing(2)+1,1) - chunk(end,1) < tolerance
                                    trial(missing(1)+1:missing(2),:) = chunk;
                                end
                            else % don't insert data if reasonable solution can't be found
                                disp('Could not find reasonable place to insert data');
                            end
                            
                            % end while loop
                            replace = 0;
                    end
                end
            end
        end
        
        % check to see that data was inserted in chronological order
        sortCheck = trial(:,1);
        sortCheck = sortCheck(~isnan(sortCheck));
        if sum(sort(sortCheck) ~= sortCheck)
            error('Data insertion was done incorrectly');
        end
        
        % calculate duration of long windows of nans (>50 ms)
%         nans = find(isnan(trial(:,1)));
%         span = diff(nans);
%         allOnes = span == 1;
%         len = CountOnes(allOnes);
%         longDrops{k} = len(len > 3)*(1/frameRate);
        frameDrops(k) = sum(isnan(trial(:,1)))./frameRate; % calculate frame drops for each trial
        
        % store data in output
        output{i}.time{k} = trial(:,1);
        output{i}.C{k} = trial(:,2:3);
        output{i}.begin{k} = trial(:,4);
        
        count = count + 1;
    end
    
    targetAbs = [d.targetX d.targetY];
    start = [0 0
            targetAbs(1:end-1,1) targetAbs(1:end-1,2)];
    targetRel = targetAbs-start;
    
    for j = 1:Ntrials
        theta = atan2(targetRel(j,2),targetRel(j,1));
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)];

        % rotated cursor trajectory
        C = output{i}.C{j};
        C = C - repmat(start(j,:),[size(C,1) 1]);
        Cr{j} = (R*C')';
    end
    
    output{i}.Cr = Cr;
    output{i}.targetAbs = targetAbs;
    output{i}.targetRel = targetRel;
    output{i}.frameRate = frameRate;
    output{i}.frameDrops = frameDrops;
%     output{i}.longDrops = longDrops;
    output{i}.OS = d.OS{1};
    output{i}.browser = d.browser{1};
    output{i}.cmConvert = d.cmConvert(1);
end
end

% compute the number of consecutive ones in an input vector
function len = CountOnes(v)
n = length(v);
len = zeros(1, ceil(n/2));%, 'uint32');
j = 0;
k = 1;
while k <= n
    if v(k)
        a = k;
        k = k + 1;
        while k <= n && v(k)
            k = k + 1;
        end
        j = j + 1;
        len(j) = k - a;
    end
    k = k + 1;
end
len = len(1:j);
end