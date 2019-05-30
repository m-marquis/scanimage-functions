

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_05_25_exp_1';
% parentDir = 'E:\Michael\2019_05_25_exp_1';


%% LOAD DATA FROM EXISTING FILE

roiNames = ["stim", "stimCtrl", "imgCtrl", "TypeD", "ANT", "LH", "PPM2"];

[fileName, pathName] = uigetfile(fullfile(parentDir, '*_SI_data.mat'));
load(fullfile(pathName, fileName), 'blockNum', 'cycleCounts', 'nDataFiles', 'roiDataAvg', ...
        'roiMetadata', 'siData');

stimRoiData = roiDataAvg(:, roiMetadata.scanRoiNums(1));
ctrlRoiData = roiDataAvg(:, roiMetadata.scanRoiNums(2));

%% PLOT ROI BOUNDS ON THEIR REFERENCE IMAGES

scanRois = roiMetadata.scanRois;
scanRoiZs = [scanRois.zDepth];
refImgZs = roiMetadata.refImgZs;

% Calculate axes bounds
minCorner = min(roiMetadata.refImgCP);
maxCorner = max(roiMetadata.refImgCP);
xRange = [minCorner(1), maxCorner(1)];
yRange = [minCorner(2), maxCorner(2)];

% Set up figure
nPlots = numel(roiMetadata.refImgZs);
if nPlots == 3
    plotLayout = [2 2];
else
    plotLayout = numSubplots(nPlots);
end
f = figure(1);clf;
f.Color = [1 1 1];

% Plot ref images and rois
for iSlice = 1:nPlots
   subaxis(plotLayout(1), plotLayout(2), iSlice, ...
            'm', 0.03, 'sv', 0.05, 'sh', 0.01, 'holdaxis', 1);
   axis equal, hold on    
   currImg = imgaussfilt(roiMetadata.refImgStack(:,:,iSlice));
   imshow(currImg, [0, max(currImg(:)) * 0.5], 'XData', xRange, 'YData', yRange)
   xlim(xRange); ylim(yRange);
   titleStr = '';
   cm = jet(numel(scanRois));
   for iRoi = 1:numel(scanRois)
       if scanRoiZs(iRoi) == refImgZs(iSlice)
           currRoi = roiMetadata.scanRois(iRoi);
           roiCenter = [currRoi.centerX, currRoi.centerY];
           roiSize = mean([currRoi.sizeX, currRoi.sizeY])/2;
           viscircles(roiCenter, roiSize, 'color', cm(iRoi, :));
           titleStr = [titleStr, ' — ', roiNames{iRoi}, ' (', num2str(iRoi), ') —'];
       end
   end
   title(titleStr)
end





%% Plot fluorescence from photostim/control ROIs 

f = figure(1);clf; hold on
nCyclesTotal = size(roiDataAvg, 1);
xData = header.frameDuration:header.frameDuration:(header.frameDuration * nCyclesTotal);
plot(stimRoiData, 'Color', 'r');
plot(ctrlRoiData, 'Color', 'b');
legend('Photostim', 'Control', 'autoupdate', 'off')

% FIND STIM ON/OFF CYCLES

manualThresh = 930;

stimCycles = stimRoiData > manualThresh;
stimCyclesStr = regexprep(num2str(stimCycles'), ' ', '');
stimOnCycles = regexp(stimCyclesStr, '(?<=0)1');
stimOffCycles = regexp(stimCyclesStr, '(?<=1)0');

% Plot to verify that they're correct
yVal = 911;
figure(f);
stimOnXData = xData(stimOnCycles);
stimOffXData = xData(stimOffCycles);
plot(stimOnCycles, ones(numel(stimOnCycles)) * yVal, 'o', 'color', 'g')
plot(stimOffCycles, ones(numel(stimOffCycles)) * yVal, '*', 'color', 'm')







%% LOAD ANATOMY STACK SLICES AND METADATA
[fileName, filePath] = uigetfile(fullfile(parentDir, '*stack*.mat'));

stackData = load(fullfile(filePath, fileName));
nSlices = stackData.header.SI.hStackManager.numSlices;
zStepSize = stackData.header.SI.hStackManager.stackZStepSize;
zSliceDepth = (0:nSlices - 1) * zStepSize;

roiRefImages = stackData.tifData(:,:,(roiMetadata.zDepth(roiMetadata.scanRoiNums) ./ zStepSize));
for iPlane = 1:size(roiRefImages, 3)
    figure(iPlane); clf; 
    imshow(roiRefImages(:,:,iPlane), [0 2000])
end

%%
currRoi =5;
currRefImage = roiRefImages(:,:,currRoi);

figure(1);clf;hold on
w = size(currRefImage, 2);
h = size(currRefImage, 1);
xPos = [-w/2, w/2];
yPos = [h/2, -h/2];
imagesc(xPos, yPos, currRefImage);colormap('gray')
axis equal
centerXYRel = roiMetadata.centerXY(roiMetadata.scanRoiNums(currRoi), :);
sizeXY = roiMetadata.sizeXY(roiMetadata.scanRoiNums(currRoi), :);

conversionFactor = 15/.59;
centerXpx = centerXYRel(:,1) * conversionFactor;
centerYpx = centerXYRel(:,2) * (-conversionFactor);
plot(centerXpx, centerYpx, 'o');

%% LOAD BEHAVIOR DATA

sid = 0;
annotFileName = 'autoAnnotations.mat';
annotData = load(fullfile(parentDir, ['sid_', num2str(sid)], annotFileName));


%% DIVIDE DATA INTO TRIALS

% Check stim durations
stimCycleDurs = stimOffCycles - stimOnCycles
interStimDurs = [stimOnCycles(1), stimOnCycles(2:end) - stimOffCycles(1:end-1), ...
        nCyclesTotal - stimOffCycles(end)]
    
skipCycles = [2];
analysisWindow = 80;
smWin = 5;

maxStimDur = max(stimCycleDurs);

% Identify analysis cycles for each stim
analysisStartCycles = []; analysisEndCycles = [];
for iStim = 1:numel(stimCycleDurs)
   if ~ismember(iStim, skipCycles)
       analysisStartCycles(end + 1) = stimOnCycles(iStim) - analysisWindow;
       analysisEndCycles(end + 1) = stimOffCycles(iStim) + analysisWindow - 1 + (maxStimDur - stimCycleDurs(iStim));
   end
end

allStimData = [];
for iStim = 1:numel(analysisStartCycles)
    for iRoi = 1:numel(scanRoiNums)
        allStimData(:, iStim, iRoi) = roiDataAvg(analysisStartCycles(iStim):analysisEndCycles(iStim), ...
                scanRoiNums(iRoi)); % --> [cycle, stim, roi]
    end
end

for iRoi = 1:numel(scanRoiNums)
    currData = allStimData(:,:,iRoi); % --> [cycle, stim]
    figure(iRoi*1000); clf; hold on;
    plot(movmean(currData, smWin, 1))
    plot(movmean(mean(currData, 2), smWin), 'linewidth', 2, 'color', 'k')
    title(num2str(iRoi))
end








