

expDate = '2019_06_14_exp_1';
sid = 0;

parentDir = fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate);
% parentDir = 'E:\Michael\2019_05_25_exp_1';

saveDir = fullfile(parentDir, ['sid_', num2str(sid)], 'Analysis');
saveDateStr = regexprep(expDate, 'exp_.*', '');
saveDateStr = regexprep(saveDateStr, '\_', '');
saveDateStr = [saveDateStr, '-', expDate(end)];

%% LOAD DATA FROM ALL BLOCKS

roiNames(1) = repmat({["stim", "stimCtrl", "TypeF", "imgCtrl"]}, 1, 1);
roiNames(2) = repmat({["stim", "stimCtrl", "PPM2soma-1", "PPM2soma-2"]}, 1, 1);
roiNames(3:5) = repmat({["stim", "stimCtrl", "TypeF", "imgCtrl"]}, 1, 3);
roiNames(6:10) = repmat({["stim", "stimCtrl", "TypeF", "imgCtrl"]}, 1, 5);

odorStimRelTimes = {[], [], [], [], [], [], [], [], [], []};

try
    
if exist(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'), 'file')
    load(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'))
    disp('Loaded existing block data file')
else
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
        allBlockData(iBlock).odorStimRelTimes = odorStimRelTimes{iBlock};
    end
    
%     
% % DROPPING BLOCK 15 BECAUSE THERE'S NO VIDEO AND 17 BECAUSE THERE'S ONLY 2 TRIALS
% allBlockData([14 16]) = [];

    
    save(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'), 'allBlockData');
end

catch foldME; rethrow(foldME); end
%% PLOT ROI BOUNDS ON THEIR REFERENCE IMAGES

saveFig = 1;

try
for iBlock = 1:numel(allBlockData)
    
    roiMetadata = allBlockData(iBlock).roiMetadata;
% roiMetadata.refImgCP = roiMetadata.refimgCP; % Only necessary for first couple of experiments
    currBlock = allBlockData(iBlock).blockNum;

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
    f = figure(iBlock + 50);clf;
    f.Color = [1 1 1];
    f.Position = [-1600 100 1300 800];
    
    % Plot ref images and rois
    for iSlice = 1:nPlots
        subaxis(plotLayout(1), plotLayout(2), iSlice, ...
            'm', 0.03, 'sv', 0.05, 'sh', 0.01, 'holdaxis', 1);
        axis equal, hold on
        currImg = imgaussfilt(roiMetadata.refImgStack(:,:,iSlice));
        currImg = currImg - min(currImg(:));
        imshow(currImg, [0, max(currImg(:)) * 0.5], 'XData', xRange, 'YData', yRange)
        xlim(xRange); ylim(yRange);
        titleStr = ['Block ', num2str(currBlock), '  (Z = ', ...
                num2str(roiMetadata.refImgZs(iSlice)), ' um)'];
        cm = jet(numel(scanRois));
        for iRoi = 1:numel(scanRois)
            if scanRoiZs(iRoi) == refImgZs(iSlice)
                currRoi = roiMetadata.scanRois(iRoi);
                roiCenter = [currRoi.centerX, currRoi.centerY];
                roiSize = mean([currRoi.sizeX, currRoi.sizeY])/2;
                viscircles(roiCenter, roiSize, 'color', cm(iRoi, :));
                titleStr = [titleStr, ' — ', allBlockData(iBlock).roiNames{iRoi}, ' (', ...
                        num2str(iRoi), ')'];
            end
        end
        title(titleStr)
    end
    
    if saveFig
        save_figure(f, saveDir, [saveDateStr, '_ref_img_plots_block_', ...
                num2str(allBlockData(iBlock).blockNum)]);
    end
    
end%iBlock

catch foldME; rethrow(foldME); end
%% LOAD BEHAVIOR DATA

annotFileName = 'autoAnnotations.mat';
annotData = load(fullfile(parentDir, ['sid_', num2str(sid)], annotFileName));
FRAME_RATE = annotData.frameInfo.FRAME_RATE;

blockBounds = [1:20:size(annotData.flowArr, 2)-1, size(annotData.flowArr, 2)]; 
% blockBounds = [1, 21, 41, 80]; 
%     blockBounds(end) = [];

try
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

catch foldME; rethrow(foldME); end
%% IDENTIFY STIM ON- and OFFSETS

currBlock = 5;
manualThresh = 10;

disp([allBlockData.blockNum])

try
currBlockData = allBlockData([allBlockData.blockNum] == currBlock);
scanRois = currBlockData.roiMetadata.scanRoiNums;
roiDataAvg = currBlockData.roiDataAvg - min(currBlockData.roiDataAvg(:)); % Setting min to zero

% roiDataAvg = roiDataAvg - 775; roiDataAvg(roiDataAvg < 0) = 0;

stimRoiData = roiDataAvg(:, scanRois(1));
ctrlRoiData = roiDataAvg(:, scanRois(2));
siData = currBlockData.siData;
cycleRate = currBlockData.cycleRate;

f = figure(100);clf; hold on
nCyclesTotal = size(roiDataAvg, 1);
xData = siData.frameDuration:siData.frameDuration:(siData.frameDuration * nCyclesTotal);
plot(stimRoiData, 'Color', 'r');
plot(ctrlRoiData, 'Color', 'b');
legend('Photostim', 'Control', 'autoupdate', 'off')

% FIND STIM ON/OFF CYCLES
stimCycles = stimRoiData > manualThresh;
stimCyclesStr = regexprep(num2str(stimCycles'), ' ', '');
stimOnCycles = regexp(stimCyclesStr, '(?<=0)1');
stimOffCycles = regexp(stimCyclesStr, '(?<=1)0');

% Plot to verify that they're correct
yVal = manualThresh;
figure(f);
stimOnXData = xData(stimOnCycles);
stimOffXData = xData(stimOffCycles);
plot(stimOnCycles, ones(numel(stimOnCycles)) * yVal, 'o', 'color', 'g')
plot(stimOffCycles, ones(numel(stimOffCycles)) * yVal, '*', 'color', 'm')

% Check stim durations
stimCycleDurs = stimOffCycles - stimOnCycles
interStimDurs = [stimOnCycles(1), stimOnCycles(2:end) - stimOffCycles(1:end-1), ...
        nCyclesTotal - stimOffCycles(end)]

catch foldME; rethrow(foldME); end
%% DIVIDE DATA INTO INDIVIDUAL STIM EPOCHS

% Check stim durations
stimCycleDurs = stimOffCycles - stimOnCycles
interStimDurs = [stimOnCycles(1), stimOnCycles(2:end) - stimOffCycles(1:end-1), ...
        nCyclesTotal - stimOffCycles(end)]
    
skipCycles = [];
analysisWindow = [44];
targetStimDur = 10;
smWin = 3;

try
if ~isempty(skipCycles)
   stimOnCycles(skipCycles) = [];
   stimOffCycles(skipCycles) = [];
   stimCycleDurs(skipCycles) = [];
   nCyclesTotal = nCyclesTotal - numel(skipCycles);
end

if numel(analysisWindow) == 1
   analysisWindow = [analysisWindow analysisWindow]; 
end
scanRoiNums = currBlockData.roiMetadata.scanRoiNums;

% Identify analysis cycles for each stim
analysisStartCycles = []; analysisEndCycles = [];
for iStim = 1:numel(stimCycleDurs)
       analysisStartCycles(end + 1) = stimOnCycles(iStim) - analysisWindow(1);
       analysisEndCycles(end + 1) = stimOffCycles(iStim) + analysisWindow(2) - 1 + ...
            (targetStimDur - stimCycleDurs(iStim));
end
analysisStartCycles(analysisStartCycles < 1) = 1;
analysisEndCycles(analysisEndCycles > nCyclesTotal) = nCyclesTotal;
analysisWinSize = max(analysisEndCycles - analysisStartCycles) + 1;

% Identify the behavior frames corresponding to the analysis cycles
if ismember('annotData', fieldnames(allBlockData))
    currAnnotData = currBlockData.annotData.trialAnnotations;
    currAnnotData(1,:) = currAnnotData(2,:); % To get rid of empty frame at the start of each trial
    rsAnnotData = currAnnotData(:);
    rsFtData = [];
    rsFtData(:,1) = currBlockData.annotData.ftData.moveSpeed(:) * 4.5 * FRAME_RATE; % --> mm/sec
    rsFtData(:,2) = currBlockData.annotData.ftData.yawSpeed(:) * FRAME_RATE; % --> rad/sec
    cyc2frame = sample_lookup(FRAME_RATE, cycleRate);
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
end

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
    if ismember('annotData', fieldnames(allBlockData))
        currAnnotData = rsAnnotData(analysisStartFrames(iStim):analysisEndFrames(iStim));
        currFtData = rsFtData(analysisStartFrames(iStim):analysisEndFrames(iStim), :); % -->[frame, var]
        preStimFrames = stimOnFrames(iStim) - analysisStartFrames(iStim);
        postStimFrames = analysisEndFrames(iStim) - stimOffFrames(iStim) - ...
            (cyc2frame.convert(targetStimDur) - stimFrameDurs(iStim)) + 1;
        if preStimFrames < cyc2frame.convert(analysisWindow(1))
            currAnnotData = [nan(cyc2frame.convert(analysisWindow(1)) - preStimFrames, 1); ...
                currAnnotData];
            currFtData = cat(1, nan(cyc2frame.convert(analysisWindow(1)) - preStimFrames, ...
                size(currFtData, 2)), currFtData);
        end
        if postStimFrames < cyc2frame.convert(analysisWindow(2))
            currAnnotData = [currAnnotData; nan(cyc2frame.convert(analysisWindow(2)) ...
                - postStimFrames, 1)];
            currFtData = cat(1, currFtData, nan(cyc2frame.convert(analysisWindow(1)) - postStimFrames, ...
                size(currFtData, 2)));
        end
        allStimAnnotData(:, iStim) = currAnnotData;
        allStimFtData = cat(3, allStimFtData, currFtData); % --> [frame, var, stim]
    end
    
end% iStim
allStimFlData = permute(allStimFlData, [1 3 2]); % --> [cycle, stim, roi]
allStimFtData = permute(allStimFtData, [1 3 2]); % --> [frame, stim, var]

disp('Data separated into stim epochs')

catch foldME; rethrow(foldME); end
%% PLOT AVERAGED FLUORESCENCE AND MOVE SPEED OVERLAYS

saveFig = 0;

try
    
% GCaMP data
nPlots = size(allStimFlData, 3);
if nPlots == 3
    plotLayout = [2 2];
else
    plotLayout = numSubplots(nPlots);
end
f = figure(101);clf;
f.Color = [1 1 1];
f.Position = [-1500 50 1500 800];

% Plot data
xData = (1/cycleRate):(1/cycleRate):size(allStimFlData, 1) / cycleRate;
for iRoi = 1:nPlots
    ax = subaxis(plotLayout(1), plotLayout(2), iRoi, 'sh', 0.05, 'sv', 0.1, 'm', 0.07); hold on
    currData = allStimFlData(:,:,iRoi); % --> [cycle, stim]
    plot_ROI_data(ax, currData, 'stdDevShading', true, 'volOffset', -analysisWindow(1), ...
        'singleTrialAlpha', 0.5, 'VolumeRate', cycleRate, ...
        'YAxisLabel', 'Raw F', 'OutlierSD', 1000);
    if ~isempty(currBlockData.odorStimRelTimes)
        plot_stim_shading(currBlockData.odorStimRelTimes)
        legend(ax.Children(1), 'Odor stim');
    end
    title(['Block ', num2str(currBlock), '  —  ', currBlockData.roiNames{iRoi}])
    xlabel('Time (sec)')
end
if saveFig
    save_figure(f, saveDir, [saveDateStr, '_avg_fl_block_', num2str(currBlock)]);
end


% FicTrac moveSpeed plots
if ismember('annotData', fieldnames(allBlockData))
    smWin = 5;
    f = figure(102);clf
    f.Color = [1 1 1];
    f.Position = [-1500 50 1000 800];
    
    smMoveSpeed = smoothdata(allStimFtData(:,:,1)', 2, 'gaussian', smWin);
    smMoveSpeedNoQui = smMoveSpeed; smMoveSpeedNoQui(allStimAnnotData' ~= 3) = nan;
    
    % First plot with quiescence included
    ax = subaxis(2,1,1, 'sv', 0.15);hold on
    plot(mean(smMoveSpeed, 'omitnan'), ...
        'linewidth', 2, 'color', 'k')
    xlabel('Time (s)')
    ylabel('Movement speed (mm/sec)')
    title(['Block ', num2str(currBlock), '  —  Average movement speed (quiescence excluded)'])
    ax.FontSize = 14;
    ax.Title.FontSize = 12;
    ax.XLabel.FontSize = 14;
    stimStart = analysisWindowFrames(1);
    stimEnd = stimStart + cyc2frame.convert(targetStimDur);
    ax.XTick = (stimStart:FRAME_RATE:(size(smMoveSpeed, 2) + stimStart)) - ...
            (round(stimStart/FRAME_RATE) * FRAME_RATE);
    ax.XTickLabel = (-round(stimStart/FRAME_RATE)):1:(size(smMoveSpeed,2)/FRAME_RATE);
    plot_stim_shading([stimStart, stimEnd]);
    legend(ax.Children(1), 'Opto stim');
    if ~isempty(currBlockData.odorStimRelTimes)
        odorShadeFrames = stimStart + (currBlockData.odorStimRelTimes * FRAME_RATE);
        plot_stim_shading(odorShadeFrames, 'Color', [0 1 0])
        legend(ax.Children(1:2), {'Odor stim', 'Opto stim'});
    end
    
    % Second plot excluding quiescence
    ax = subaxis(2,1,2);hold on
    plot(mean(smMoveSpeedNoQui, 'omitnan'), ...
            'linewidth', 2, 'color', 'k')
    xlabel('Time (s)')
    ylabel('Movement speed (mm/sec)')
    title(['Block ', num2str(currBlock), '  —  Average movement speed (quiescence excluded)'])
    ax.FontSize = 14;
    ax.Title.FontSize = 12;
    ax.XLabel.FontSize = 14;
    stimStart = analysisWindowFrames(1);
    stimEnd = stimStart + cyc2frame.convert(targetStimDur);
    ax.XTick = (stimStart:FRAME_RATE:(size(smMoveSpeed, 2) + stimStart)) - ...
            (round(stimStart/FRAME_RATE) * FRAME_RATE);
    ax.XTickLabel = (-round(stimStart/FRAME_RATE)):1:(size(smMoveSpeed,2)/FRAME_RATE);
    plot_stim_shading([stimStart, stimEnd]);
    legend(ax.Children(1), 'Opto stim');
    if ~isempty(currBlockData.odorStimRelTimes)
        odorShadeFrames = stimStart + (currBlockData.odorStimRelTimes * FRAME_RATE);
        plot_stim_shading(odorShadeFrames, 'Color', [0 1 0])
        legend(ax.Children(1:2), {'Odor stim', 'Opto stim'});
    end
    
    if saveFig
        save_figure(f, saveDir, [saveDateStr, '_avg_moveSpeed_block_', ...
                num2str(currBlockData.blockNum)])
    end
end

catch foldME; rethrow(foldME); end
%% CREATE 2D PLOTS OF BEHAVIOR AND GCaMP DATA

saveFig = 0;

try
if ismember('annotData', fieldnames(allBlockData))
    
    % Behavior data
    f = figure(103);clf;
    ax = axes();
    f.Color = [1 1 1];
    f.Position = [-1050 45 800 700];
    cMap = [rgb('Indigo'); ...
            rgb('Orange'); ...
            rgb('Green');
            rgb('Cyan'); ...
            ];
    allStimAnnotData(end, end) = 4; % to keep the color mapping consistent
    plot_2D_summary(allStimAnnotData', FRAME_RATE, 'plotAxes', ax, 'colormap', cMap);
    hold on;
    ax.Title.String = {['Block ', num2str(currBlock) , ' behavior annotation'], ...
            'Indigo: quiescence,  cyan: locomotion,  orange: isolated movement'};
    ax.FontSize = 14;
    ax.Title.FontSize = 14;
    ax.XLabel.FontSize = 14;
    stimStart = analysisWindowFrames(1);
    stimEnd = stimStart + cyc2frame.convert(targetStimDur);
    ax.XTick = (stimStart:FRAME_RATE:(size(allStimAnnotData, 1) + stimStart)) - ...
            (round(stimStart/FRAME_RATE) * FRAME_RATE);
    ax.XTickLabel = (-round(stimStart/FRAME_RATE)):1:(size(allStimAnnotData,1)/FRAME_RATE);
    plot(ax, [stimStart, stimStart], ylim(ax), 'Color', 'r', 'linewidth', 2);
    plot(ax, [stimEnd, stimEnd], ylim(ax), 'Color', 'r', 'linewidth', 2);
    if ~isempty(currBlockData.odorStimRelTimes)
        odorShadeFrames = stimStart + (currBlockData.odorStimRelTimes * FRAME_RATE);
        plot(ax, [odorShadeFrames(1), odorShadeFrames(1)], ylim(ax), 'Color', 'g', 'linewidth', 2);
        plot(ax, [odorShadeFrames(2), odorShadeFrames(2)], ylim(ax), 'Color', 'g', 'linewidth', 2);
        legend(ax.Children([1 3]), {'Odor stim', 'Opto stim'});
    end
    if saveFig
        save_figure(f, saveDir, [saveDateStr, '_2D_annotData_block_', ...
                num2str(currBlockData.blockNum)])
    end
    
    % FicTrac data
    smWin = 5;
    speedData = smoothdata(allStimFtData(:,:,1)', 2, 'gaussian', smWin); % --> [frame, stim]
    f = figure(104);clf;
    ax = axes();
    f.Color = [1 1 1];
    f.Position = [-1050 45 800 700];
    % speedData(speedData > 15) = 15;
    titleStr = ['Block ', num2str(currBlock), '  —  Move Speed (mm/sec)'];
    plot_2D_summary(speedData, FRAME_RATE, 'plotAxes', ax, 'titleStr', titleStr);
    hold on
    colorbar
    ax.FontSize = 14;
    ax.Title.FontSize = 14;
    ax.XLabel.FontSize = 14;
    stimStart = analysisWindowFrames(1);
    stimEnd = stimStart + cyc2frame.convert(targetStimDur);
    ax.XTick = (stimStart:FRAME_RATE:(size(speedData, 2) + stimStart)) - ...
        (round(stimStart/FRAME_RATE) * FRAME_RATE);
    ax.XTickLabel = (-round(stimStart/FRAME_RATE)):1:(size(speedData,2)/FRAME_RATE);
    plot(ax, [stimStart, stimStart], ylim(ax), 'Color', 'r', 'linewidth', 2);
    plot(ax, [stimEnd, stimEnd], ylim(ax), 'Color', 'r', 'linewidth', 2);
    if ~isempty(currBlockData.odorStimRelTimes)
        odorShadeFrames = stimStart + (currBlockData.odorStimRelTimes * FRAME_RATE);
        plot(ax, [odorShadeFrames(1), odorShadeFrames(1)], ylim(ax), 'Color', 'g', 'linewidth', 2);
        plot(ax, [odorShadeFrames(2), odorShadeFrames(2)], ylim(ax), 'Color', 'g', 'linewidth', 2);
        legend(ax.Children([1 3]), {'Odor stim', 'Opto stim'});
    end
    if saveFig
        save_figure(f, saveDir, [saveDateStr, '_2D_moveSpeed_block_', ...
                num2str(currBlockData.blockNum)])
    end
end

% GCaMP data
smWin = 5;
for iRoi = 1:size(allStimFlData, 3)
    currRoiName = currBlockData.roiNames{iRoi};
    currFlData = allStimFlData(:,:,iRoi); % --> [cycle, stim]
    currFlData = smoothdata(currFlData, 1, 'gaussian', smWin, 'includenan')';
    f = figure(iRoi + 110);clf;
    ax = axes();
    f.Color = [1 1 1];
    f.Position = [-1050 45 800 700];
    titleStr = ['Block ', num2str(currBlock), '  —  ', currRoiName];
    plot_2D_summary(currFlData, cycleRate, 'titleStr', titleStr, 'plotAxes', ax);
    hold on; colorbar;
    ax.FontSize = 14;
    ax.Title.FontSize = 14;
    ax.XLabel.FontSize = 14;
    stimStart = analysisWindow(1);
    stimEnd = stimStart + targetStimDur;
    ax.XTick = (stimStart:cycleRate:(size(currFlData, 2) + stimStart)) - ...
            (round(stimStart/cycleRate) * cycleRate);
    ax.XTickLabel = (-round(stimStart/cycleRate)):1:(size(currFlData,2)/cycleRate);
    plot(ax, [stimStart, stimStart], ylim(ax), 'Color', 'r', 'linewidth', 2);
    plot(ax, [stimEnd, stimEnd], ylim(ax), 'Color', 'r', 'linewidth', 2);
    if ~isempty(currBlockData.odorStimRelTimes)
    odorShadeCycles = stimStart + (currBlockData.odorStimRelTimes * cycleRate);
    plot(ax, [odorShadeCycles(1), odorShadeCycles(1)], ylim(ax), 'Color', 'g', 'linewidth', 2);
    plot(ax, [odorShadeCycles(2), odorShadeCycles(2)], ylim(ax), 'Color', 'g', 'linewidth', 2);
    legend(ax.Children([1 3]), {'Odor stim', 'Opto stim'});
end
    if saveFig
        save_figure(f, saveDir, [saveDateStr, '_2D_fl_', currRoiName, '_block_', ...
                num2str(currBlockData.blockNum)])
    end
end

catch foldME; rethrow(foldME); end


%%

currRoi = 4;
smWin = 5;

cycleRate = allBlockData(currBlock).cycleRate;
ftDataLin = as_vector(allStimFtData(:,:,1));
flDataLin = as_vector(allStimFlData(:,:, currRoi));
frameX = (1/FRAME_RATE):(1/FRAME_RATE):(numel(ftDataLin)/FRAME_RATE);
cycleX = (1/cycleRate):(1/cycleRate):(numel(flDataLin)/cycleRate);

ftDataSm = smoothdata(ftDataLin,  1, 'gaussian', smWin);
flDataSm = smoothdata((flDataLin - min(flDataLin)),  1, 'gaussian', smWin);

ftDataNorm = ftDataSm / max(ftDataSm);
flDataNorm = flDataSm / max(flDataSm);

figure(1); clf; hold on; ax = gca();
plot_stim_shading
plot(frameX, ftDataNorm, 'linewidth', 2);
plot(cycleX, flDataNorm, 'linewidth', 2);





