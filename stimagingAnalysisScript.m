
% Load tif data

tifDir = dir('E:\Michael\2019_02_04_exp_1\*sid_5_bid_0*.tif');
nTrials = numel(tifDir);

rawData = [];
for iTrial = 1:nTrials
    disp(num2str(iTrial))
    tifPath = fullfile(tifDir(iTrial).folder, tifDir(iTrial).name);    
    [~, rawData{iTrial}, ~] = opentif(tifPath);
       
end

% Extract fluorescence
meanLum = []; frameCounts = [];
for iTrial = 1:nTrials
    disp(iTrial)
    meanLum{iTrial} = squeeze(mean(mean(rawData{iTrial}, 1), 2));
    frameCounts(iTrial) = numel(meanLum{iTrial});
end
minFrames = min(frameCounts)

% Trim extra frames off
flData = zeros(minFrames, nTrials);
for iTrial = 1:nTrials 
    flData(:, iTrial) = meanLum{iTrial}(1:minFrames); % --> [frame, trial]
end

% Save tif data as a .mat file
saveName = regexprep(tifPath, '_00...\.tif', '_flData.mat');
if ~exist(saveName, 'file')
    save(saveName, 'flData');
end

%%
imgData = flData(1:2:end, :);
stimData = flData(2:2:end, :);
figure(7);clf; plot(imgData);
figure(8);clf; plot(stimData);

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
    end
end
figure(10); clf; imagesc(badFrames);

%% Remove bad frames

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

% figure(2);clf; hold on; 
% for iTrial = trialList
%    plot(optoStimInfo.powerLog{iTrial}+iTrial) 
% end

figure(3);clf; hold on; set(gcf, 'Color', [1 1 1])
plot(imgData(startFrame:end, trialList));
plot(mean(imgData(startFrame:end, trialList), 2, 'omitnan'), 'linewidth', 4, 'color', 'k')
yL = ylim();
rectPos = [stimFrames(1), yL(1), diff(stimFrames), diff(yL)];
rectangle('Position', rectPos, 'EdgeColor', 'none', 'faceColor', [1 0 0 0.1] )

figure(4);clf; hold on; set(gcf, 'Color', [1 1 1])
plot(stimData(startFrame:end, trialList));
plot(mean(stimData(startFrame:end, trialList), 2, 'omitnan'), 'linewidth', 4, 'color', 'k')
yL = ylim();
rectPos = [stimFrames(1), yL(1), diff(stimFrames), diff(yL)];
rectangle('Position', rectPos, 'EdgeColor', 'none', 'faceColor', [1 0 0 0.1] )

