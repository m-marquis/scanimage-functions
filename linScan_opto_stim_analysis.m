

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_05_25_exp_1';
% parentDir = 'E:\Michael\2019_05_25_exp_1';


%% LOAD DATA FROM ALL BLOCKS

roiNames(1:3) = repmat({["stim", "stimCtrl", "imgCtrl", "SLP", "ANT", "LH", "PPM2"]}, 1, 3);
roiNames(4:6) = repmat({["stim", "stimCtrl", "c-ANT", "i-ANT", "SLP"]}, 1, 3)
roiNames(7) = repmat({["stim", "stimCtrl", "LH", "SLP"]}, 1, 1);
roiNames(8) = repmat({["stim", "stimCtrl", "MB-CA", "SLP"]}, 1, 1);
roiNames(9:16) = repmat({["stim", "stimCtrl", "SLP"]}, 1,  8);

blockDataFiles = dir(fullfile(parentDir, '*_SI_data.mat'));
nBlocks = numel(blockDataFiles);
clear allBlockData;
for iBlock = 1:nBlocks
   allBlockData(iBlock) = load(fullfile(parentDir, blockDataFiles(iBlock).name)); 
   % .blockNum, .cycleCounts, .nDataFiles, .roiDataAvg, .roiMetadata, .siData
   disp(['Block ', num2str(allBlockData(iBlock).blockNum), ' loaded']);
end
for iBlock = 1:nBlocks
    allBlockData(iBlock).roiNames = roiNames{iBlock};
    allBlockData(iBlock).cycleRate = 1/allBlockData(iBlock).siData.frameDuration;
end

        % DROPPING BLOCK 15 BECAUSE THERE'S NO VIDEO AND 17 BECAUSE THERE'S ONLY 2 TRIALS
        allBlockData([14 16]) = [];

save(fullfile(parentDir, 'allBlockData.mat'), 'allBlockData');

%% PLOT ROI BOUNDS ON THEIR REFERENCE IMAGES

for iBlock = 1:numel(allBlockData)
    
    roiMetadata = allBlockData(iBlock).roiMetadata;
    roiMetadata.refImgCP = roiMetadata.refimgCP;
    
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
    f = figure(iBlock);clf;
    f.Color = [1 1 1];
    f.Position = [-1600 100 1300 800];
    
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
                titleStr = [titleStr, ' — ', allBlockData(iBlock).roiNames{iRoi}, ' (', ...
                        num2str(iRoi), ') —'];
            end
        end
        title(titleStr)
    end
    
    save_figure(f, parentDir, ['Ref_img_plots_block_', num2str(allBlockData(iBlock).blockNum)]);
    
    
end%iBlock

%% LOAD BEHAVIOR DATA

sid = 0;
annotFileName = 'autoAnnotations.mat';
annotData = load(fullfile(parentDir, ['sid_', num2str(sid)], annotFileName));
FRAME_RATE = annotData.frameInfo.FRAME_RATE;

blockBounds = [1:20:size(annotData.flowArr, 2)-1, size(annotData.flowArr, 2)]; 
    blockBounds(end) = [];
    
for iBlock = 1:numel(allBlockData)
    blockInds = blockBounds(iBlock):blockBounds(iBlock + 1) - 1;
    
    allBlockData(iBlock).annotData.annotParams = annotData.annotParams;
    allBlockData(iBlock).annotData.behaviorLabels = annotData.behaviorLabels;
    allBlockData(iBlock).annotData.flowArr = annotData.flowArr(:, blockInds);
    allBlockData(iBlock).annotData.frameInfo = annotData.frameInfo;
    allBlockData(iBlock).annotData.goodTrials = annotData.goodTrials(:, blockInds);
    allBlockData(iBlock).annotData.trialAnnotations = annotData.trialAnnotations(blockInds, :)';
        
    fieldNames = fieldnames(annotData.ftData);    
    for iField = 1:numel(fieldNames)
       sz = size(annotData.ftData.(fieldNames{iField}));
       trialDim = find(sz == size(annotData.flowArr, 2));
       if trialDim == 2
           allBlockData(iBlock).annotData.ftData.(fieldNames{iField}) = ...
               annotData.ftData.(fieldNames{iField})(:, blockInds);
       elseif trialDim == 3
           allBlockData(iBlock).annotData.ftData.(fieldNames{iField}) = ...
               annotData.ftData.(fieldNames{iField})(:, :, blockInds);
       end
    end
end

%% IDENTIFY STIM ON- and OFFSETS

currBlock = 9;
currBlockData = allBlockData([allBlockData.blockNum] == currBlock);
scanRois = currBlockData.roiMetadata.scanRoiNums;
roiDataAvg = currBlockData.roiDataAvg - min(currBlockData.roiDataAvg(:)); % Setting min to zero
stimRoiData = roiDataAvg(:, scanRois(1));
ctrlRoiData = roiDataAvg(:, scanRois(2));
siData = currBlockData.siData;


f = figure(1);clf; hold on
nCyclesTotal = size(roiDataAvg, 1);
xData = siData.frameDuration:siData.frameDuration:(siData.frameDuration * nCyclesTotal);
plot(stimRoiData, 'Color', 'r');
plot(ctrlRoiData, 'Color', 'b');
legend('Photostim', 'Control', 'autoupdate', 'off')

% FIND STIM ON/OFF CYCLES
manualThresh = 30;

stimCycles = stimRoiData > manualThresh;
stimCyclesStr = regexprep(num2str(stimCycles'), ' ', '');
stimOnCycles = regexp(stimCyclesStr, '(?<=0)1');
stimOffCycles = regexp(stimCyclesStr, '(?<=1)0');

% Plot to verify that they're correct
yVal = 11;
figure(f);
stimOnXData = xData(stimOnCycles);
stimOffXData = xData(stimOffCycles);
plot(stimOnCycles, ones(numel(stimOnCycles)) * yVal, 'o', 'color', 'g')
plot(stimOffCycles, ones(numel(stimOffCycles)) * yVal, '*', 'color', 'm')


%% DIVIDE DATA INTO INDIVIDUAL STIM EPOCHS

% Check stim durations
stimCycleDurs = stimOffCycles - stimOnCycles
interStimDurs = [stimOnCycles(1), stimOnCycles(2:end) - stimOffCycles(1:end-1), ...
        nCyclesTotal - stimOffCycles(end)]
    
skipCycles = [];
analysisWindow = [50 50];
targetStimDur = 20;
smWin = 3;

if ~isempty(skipCycles)
   stimOnCycles(skipCycles) = [];
   stimOffCycles(skipCycles) = [];
   stimCycleDurs(skipCycles) = [];
end

% maxStimDur = max(stimCycleDurs);
scanRoiNums = currBlockData.roiMetadata.scanRoiNums;

% Identify analysis cycles for each stim
analysisStartCycles = []; analysisEndCycles = [];
for iStim = 1:numel(stimCycleDurs)
%    if ~ismember(iStim, skipCycles)
       analysisStartCycles(end + 1) = stimOnCycles(iStim) - analysisWindow(1);
       analysisEndCycles(end + 1) = stimOffCycles(iStim) + analysisWindow(2) - 1 + ...
            (targetStimDur - stimCycleDurs(iStim));
%    end
end
analysisStartCycles(analysisStartCycles < 1) = 1;
analysisEndCycles(analysisEndCycles > nCyclesTotal) = nCyclesTotal;
analysisWinSize = max(analysisEndCycles - analysisStartCycles) + 1;

% Identify the behavior frames corresponding to the analysis cycles
currAnnotData = currBlockData.annotData.trialAnnotations;
currAnnotData(1,:) = currAnnotData(2,:); % To get rid of empty frame at the start of each trial
rsAnnotData = currAnnotData(:);
rsFtData = [];
rsFtData(:,1) = currBlockData.annotData.ftData.moveSpeed(:);
rsFtData(:,2) = currBlockData.annotData.ftData.yawSpeed(:);
cyc2frame = sample_lookup(FRAME_RATE, currBlockData.cycleRate);
stimOnFrames = cyc2frame.convert(stimOnCycles);
stimOffFrames = cyc2frame.convert(stimOffCycles);
stimFrameDurs = stimOffFrames - stimOnFrames;
analysisWindowFrames = cyc2frame.convert(analysisWindow);
analysisStartFrames = []; analysisEndFrames = [];
for iStim = 1:numel(stimFrameDurs)
    analysisStartFrames(end + 1) = stimOnFrames(iStim) - analysisWindowFrames(1);
    analysisEndFrames(end + 1) = stimOffFrames(iStim) + analysisWindowFrames(2) - 1 + ...
            (cyc2frame.convert(targetStimDur) - stimFrameDurs(iStim));
end
analysisStartFrames(analysisStartFrames < 1) = 1;
analysisEndFrames(analysisEndFrames > numel(rsAnnotData)) = numel(rsAnnotData);
test = max(analysisEndFrames - analysisStartFrames) + 1;

% Extract and compile data from analysis windows
allStimFlData = []; allStimAnnotData = []; allStimFtData = [];
for iStim = 1:numel(analysisStartCycles)
    
    % GCaMP data
    currStimData = roiDataAvg(analysisStartCycles(iStim):analysisEndCycles(iStim), ...
            scanRoiNums); % --> [cycle, stim, roi]
    preStimCycles = stimOnCycles(iStim) - analysisStartCycles(iStim);
    postStimCycles = analysisEndCycles(iStim) - stimOffCycles(iStim) - ...
            (targetStimDur - stimCycleDurs(iStim)) + 1;        
    if preStimCycles < analysisWindow(1)
       currStimData = cat(1, nan(analysisWindow(1) - preStimCycles, size(currStimData, 2)), ...
                currStimData);
    end
    if postStimCycles < analysisWindow(2)
       currStimData = cat(1, currStimData, nan(analysisWindow(2) - postStimCycles, ...
                size(currStimData, 2)));  
    end
    allStimFlData = cat(3, allStimFlData, currStimData); % --> [cycle, roi, stim]
    
    % Behavior data
    currAnnotData = rsAnnotData(analysisStartFrames(iStim):analysisEndFrames(iStim));
    currFtData = rsFtData(analysisStartFrames(iStim):analysisEndFrames(iStim), :); % -->[frame, var]
    preStimFrames = stimOnFrames(iStim) - analysisStartFrames(iStim);
    postStimFrames = analysisEndFrames(iStim) - stimOffFrames(iStim) - ...
            (cyc2frame.convert(targetStimDur) - stimFrameDurs(iStim)) + 1;
    if preStimFrames < cyc2frame.convert(analysisWindow(1))
        currAnnotData = [nan(cyc2frame.convert(analysisWindow(1)) - preStimFrames, 1); ...
                currAnnotData];
        currFtData = cat(1, nan(cyc2frame.convert(analysisWindow(1)) - preStimFrames, ...
            size(currFtData, 2)), currAnnotData);
    end
    if postStimFrames < cyc2frame.convert(analysisWindow(2))
        currAnnotData = [currAnnotData; nan(cyc2frame.convert(analysisWindow(2)) ...
                - postStimFrames, 1)];
        currFtData = cat(1, currFtData, nan(cyc2frame.convert(analysisWindow(1)) - preStimFrames, ...
                size(currFtData, 2)));     
    end
    allStimAnnotData(:, iStim) = currAnnotData;
    allStimFtData = cat(3, allStimFtData, currFtData); % --> [frame, var, stim]

end% iStim
allStimFlData = permute(allStimFlData, [1 3 2]); % --> [cycle, stim, roi]
allStimFtData = permute(allStimFtData, [1 3 2]); % --> [frame, stim, var]

%% PLOT AVERAGED FLUORESCENCE OVERLAYS

saveFig = 0;

% Set up figure
nPlots = size(allStimFlData, 3);
if nPlots == 3
    plotLayout = [2 2];
else
    plotLayout = numSubplots(nPlots);
end
f = figure(1000);clf;
f.Color = [1 1 1];
f.Position = [-1500 50 1500 800];

% Plot ref images and rois
xData = (1/currBlockData.cycleRate):(1/currBlockData.cycleRate):size(allStimFlData, 1) ...
        / currBlockData.cycleRate;
for iRoi = 1:nPlots
    ax = subaxis(plotLayout(1), plotLayout(2), iRoi, 'sh', 0.05, 'sv', 0.1, 'm', 0.07); hold on
    currData = allStimFlData(:,:,iRoi); % --> [cycle, stim]
    ax = plot_ROI_data(ax, currData, 'stdDevShading', true, 'volOffset', -analysisWindow(1), ...
                    'singleTrialAlpha', 0.5, 'VolumeRate', currBlockData.cycleRate, ...
                    'YAxisLabel', 'Raw F', 'OutlierSD', 1000);
    title(currBlockData.roiNames{iRoi})
    xlabel('Time (sec)')
    ylabel('Raw fluorescence (AU)')
end

if saveFig
    save_figure(f, parentDir, ['Block_', num2str(currBlockData.blockNum), '_fl'])
end

%% CREATE 2D PLOTS OF BEHAVIOR AND GCaMP DATA

saveFig = 0;

infoStruct = []; 
infoStruct.expDate = '';
infoStruct.trialDuration













