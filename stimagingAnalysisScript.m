
% Load tif data

tifDir = dir('E:\Michael\2019_02_04_exp_1\*sid_5_bid_1*.tif');
nTrials = numel(tifDir);

rawData = [];
for iTrial = 1:nTrials
    disp(num2str(iTrial))
    tifPath = fullfile(tifDir(iTrial).folder, tifDir(iTrial).name);    
    [~, rawData{iTrial}, ~] = opentif(tifPath);
       
end

%% Extract fluorescence
meanLum = []; frameCounts = [];
for iTrial = 1:nTrials
    disp(iTrial)
    meanLum{iTrial} = squeeze(mean(mean(rawData{iTrial}, 1), 2));
    frameCounts(iTrial) = numel(meanLum{iTrial});
end
minFrames = min(frameCounts)

%% Trim extra frames off
flData = zeros(minFrames, nTrials);
for iTrial = 1:nTrials 
    flData(:, iTrial) = meanLum{iTrial}(1:minFrames); % --> [frame, trial]
end

%% Save tif data as a .mat file
saveName = regexprep(tifPath, '_00...\.tif', '_tifData.mat');
if ~exist(saveName, 'file')
    save(saveName, 'flData');
end

%%
imgData = flData(1:2:end, :);
stimData = flData(2:2:end, :);
figure(7);clf; plot(imgData);
figure(8);clf; plot(stimData);



% % Drop outlier frames
% imgFrameMean = mean(imgData, 2);
% stimFrameMean = mean(stimData, 2);
% imgFrameStd = std(imgData, 0, 2);
% stimFrameStd = std(stimData, 0, 2);
% imgFrameBounds = [imgFrameMean - (2*imgFrameStd), imgFrameMean + (2*imgFrameStd)];
% stimFrameBounds = [stimFrameMean - (2*stimFrameStd), stimFrameMean + (2*stimFrameStd)];
% 
% imgFrameOutliers = imgData < repmat(imgFrameBounds(:,1), 1, size(imgData, 2)) ...
%         | imgData > repmat(imgFrameBounds(:,2), 1, size(imgData, 2));
% 
% stimFrameOutliers = stimData < repmat(stimFrameBounds(:,1), 1, size(stimData, 2)) ...
%         | stimData > repmat(stimFrameBounds(:,2), 1, size(stimData, 2));
%     
% figure(5);clf; imagesc(imgFrameOutliers');
% figure(6);clf; imagesc(stimFrameOutliers');



%% Identify bad frames
badFrames = zeros(size(flData));
meanFl = mean(flData, 2);

for iTrial = 1:nTrials
    for iFrame = 2:(minFrames - 1)
        currFrame = flData(iFrame, iTrial);
        prevFrame = flData(iFrame - 1, iTrial);
        nextFrame = flData(iFrame + 1, iTrial);
        currMean = meanFl(iFrame);
        prevMean = meanFl(iFrame - 1);
        nextMean = meanFl(iFrame + 1);
        
        prevSmaller = currFrame - prevFrame > 0;
        nextSmaller = currFrame - nextFrame > 0;
        
        prevMeanSmaller = currMean - prevMean > 0;
        nextMeanSmaller = currMean - nextMean > 0;
        
        if prevSmaller ~= prevMeanSmaller || nextSmaller ~= nextMeanSmaller
           badFrames(iFrame, iTrial) = 1; 
        end
        
%         
%         if abs(diff([currFrame, prevMean])) < abs(diff([currFrame, currMean])) ...
%             || abs(diff([currFrame, nextMean])) < abs(diff([currFrame, currMean]))
%         badFrames(iFrame, iTrial) = 1;
%         end
    end
end
figure(10); clf; imagesc(badFrames);

%% Remove bad frames

% imgData(imgFrameOutliers) = nan;
% stimData(stimFrameOutliers) = nan;
% figure(7);clf; plot(imgData);
% figure(8);clf; plot(stimData);

tempData = flData;
tempData(logical(badFrames)) = nan;

imgData = tempData(1:2:end, :);
stimData = tempData(2:2:end, :);
figure(7);clf; plot(imgData);
figure(8);clf; plot(stimData);

%% Plot data

trialList = 1:nTrials;
stimFrames = [17 25];
startFrame = 3;
stimFrames = stimFrames - startFrame + 1;


figure(1);clf; hold on; 
plot(flData(:, trialList));

figure(2);clf; hold on; 
for iTrial = trialList
   plot(optoStimInfo.powerLog{iTrial}+iTrial) 
end

figure(3);clf; hold on
plot(imgData(startFrame:end, trialList));
plot(mean(imgData(startFrame:end, trialList), 2, 'omitnan'), 'linewidth', 4, 'color', 'k')
yL = ylim();
rectPos = [stimFrames(1), yL(1), diff(stimFrames), diff(yL)];
rectangle('Position', rectPos, 'EdgeColor', 'none', 'faceColor', [1 0 0 0.1] )

figure(4);clf; hold on
plot(stimData(startFrame:end, trialList));
plot(mean(stimData(startFrame:end, trialList), 2, 'omitnan'), 'linewidth', 4, 'color', 'k')
yL = ylim();
rectPos = [stimFrames(1), yL(1), diff(stimFrames), diff(yL)];
rectangle('Position', rectPos, 'EdgeColor', 'none', 'faceColor', [1 0 0 0.1] )

