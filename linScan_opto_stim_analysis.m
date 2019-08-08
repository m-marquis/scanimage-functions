

expDate = '2019_08_05_exp_2';
sid = 0;

parentDir = fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate);
% parentDir = 'E:\Michael\2019_05_25_exp_1';

saveDir = fullfile('D:\Dropbox (HMS)\2P Data\Analysis', expDate, ['sid_', num2str(sid)]);
saveDateStr = regexprep(expDate, 'exp_.*', '');
saveDateStr = regexprep(saveDateStr, '\_', '');
saveDateStr = [saveDateStr, '-', expDate(end)];

%% LOAD DATA FROM ALL BLOCKS

try
    
if exist(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'), 'file')
    load(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'))
    disp('Loaded existing block data file')
    
    FRAME_RATE = allBlockData(1).annotData.frameInfo.FRAME_RATE;
    
else
    disp('Error: no block data file found!')
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
    f.Name = ['Block ', num2str(currBlock), ' ROIs'];
    
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
                roiStartPos = roiCenter - ([abs(currRoi.sizeX), abs(currRoi.sizeY)]/2); 
                rectangle('position', [roiStartPos, abs(currRoi.sizeX), abs(currRoi.sizeY)], ...
                        'edgecolor', cm(iRoi, :), 'linewidth', 2);
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


%% SELECT CURRENT BLOCK

currBlock = 8;

currBlockData = allBlockData([allBlockData.blockNum] == currBlock);
allStimFlData = currBlockData.stimSepData.flData;
allStimFtData = currBlockData.stimSepData.ftData;
allStimAnnotData = currBlockData.stimSepData.annotData;
cycleRate = currBlockData.cycleRate;
cyc2frame = currBlockData.stimSepData.cyc2frame;
analysisWindow = currBlockData.stimSepData.analysisWindow;
analysisWindowFrames = cyc2frame.convert(analysisWindow);
stimStart = currBlockData.stimSepData.stimStart;
stimEnd = currBlockData.stimSepData.stimEnd;
stimOnCycles = currBlockData.stimSepData.stimOnCycles;
stimOffCycles = currBlockData.stimSepData.stimOffCycles;
targetStimDur = currBlockData.stimSepData.targetStimDur;
skipCycles = currBlockData.stimSepData.skipCycles;

%% PLOT AVERAGED FLUORESCENCE AND MOVE SPEED OVERLAYS

saveFig = 0;
smWin = 3;

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
f.Position = [-1918 167 1050 600];

% Plot data
xData = (1/cycleRate):(1/cycleRate):size(allStimFlData, 1) / cycleRate;
for iRoi = 1:nPlots
    ax = subaxis(plotLayout(1), plotLayout(2), iRoi, 'sh', 0.05, 'sv', 0.13, 'm', 0.07); hold on
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
% suptitle(['Stim power: ', num2str(stimPower), '%,  Imaging power: ', num2str(imgPower), '%']);
if saveFig
    save_figure(f, saveDir, [saveDateStr, '_avg_fl_block_', num2str(currBlock)]);
end


% FicTrac moveSpeed plots
if ismember('annotData', fieldnames(allBlockData))
    smWin = 5;
    f = figure(102);clf
    f.Color = [1 1 1];
    f.Position = [-997 25 1000 950];
    
    smMoveSpeed = smoothdata(allStimFtData(:,:,1)', 2, 'gaussian', smWin);
    smMoveSpeedNoQui = smMoveSpeed; smMoveSpeedNoQui(allStimAnnotData' ~= 3) = nan;
    smYawSpeed = smoothdata(allStimFtData(:,:,2)', 2, 'gaussian', smWin);
    smYawSpeedNoQui = smYawSpeed; smYawSpeedNoQui(allStimAnnotData' ~=3) = nan;
    
    % First plot with quiescence included
    ax = subaxis(3,1,1, 'sv', 0.15);hold on
    plot(mean(smMoveSpeed, 'omitnan'), ...
        'linewidth', 2, 'color', 'k')
    xlabel('Time (s)')
    ylabel('Movement speed (mm/sec)')
    title(['Block ', num2str(currBlock), '  —  Average movement speed (quiescence included)'])
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
    ax = subaxis(3,1,2);hold on
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
    
    % Third plot with yaw speed (excluding quiescence)
    ax = subaxis(3,1,3);hold on
    plot(mean(smYawSpeed, 'omitnan'), ...
            'linewidth', 2, 'color', 'k')
    xlabel('Time (s)')
    ylabel('Yaw speed (rad/sec)')
    title(['Block ', num2str(currBlock), '  —  Average yaw velocity (quiescence excluded)'])
    ax.FontSize = 14;
    ax.Title.FontSize = 12;
    ax.XLabel.FontSize = 14;
    stimStart = analysisWindowFrames(1);
    stimEnd = stimStart + cyc2frame.convert(targetStimDur);
    ax.XTick = (stimStart:FRAME_RATE:(size(smYawSpeed, 2) + stimStart)) - ...
        (round(stimStart/FRAME_RATE) * FRAME_RATE);
    ax.XTickLabel = (-round(stimStart/FRAME_RATE)):1:(size(smYawSpeed,2)/FRAME_RATE);
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

saveFig = 1;
smWin = 3;

try
if ismember('annotData', fieldnames(allBlockData))
    
    % Behavior data
    f = figure(103);clf;
    ax = axes();
    f.Color = [1 1 1];
    f.Position = [-1860 45 800 700];
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
    f.Position = [-1860 45 800 700];
    % speedData(speedData > 15) = 15;
    titleStr = ['Block ', num2str(currBlock), '  —  Move Speed (mm/sec)'];
    speedData(speedData > max(speedData(:) * 0.75)) = max(speedData(:) * 0.75);
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

groupName = 'TypeF stim';
plotBlocks = [0 1 3 5 7 8]
% plotBlocks = [0:2:allBlockData(end - 1).blockNum, allBlockData(end - 1).blockNum];
% plotBlocks = 0:2:allBlockData(end - 1).blockNum;

% groupName = 'TypeF stim';
% plotBlocks = 1:2:7

stimLines = [8 12];
saveFig = 0;
smWin = 5;

plotAnnotData = []; plotFtData = []; plotBlockBounds = 1;
for iBlock = 1:numel(plotBlocks)
    currBlockData = allBlockData([allBlockData.blockNum] == plotBlocks(iBlock)).stimSepData;
    plotAnnotData = cat(1, plotAnnotData, currBlockData.annotData');
    currFtData = smoothdata(currBlockData.ftData(:,:,1)' .* 4.5 .* FRAME_RATE, 2, 'gaussian', smWin);
    plotFtData = cat(1, plotFtData, currFtData);
%     plotAnnotData = cat(1, plotAnnotData, allBlockData([allBlockData.blockNum] == plotBlocks(iBlock)).annotData.trialAnnotations');
%     currFtData = smoothdata(allBlockData([allBlockData.blockNum] == plotBlocks(iBlock)).annotData.ftData.moveSpeed' .* ...
%             4.5 .* FRAME_RATE, 2, 'gaussian', smWin);
%     plotFtData = cat(1, plotFtData, currFtData);
   plotBlockBounds(end + 1) = plotBlockBounds(end) + size(currFtData, 1); 
end
plotBlockBounds(end) = plotBlockBounds(end - 1);


stimLines = [currBlockData.stimStart, currBlockData.stimEnd] ./ FRAME_RATE;


% Annotations
f = figure(201); clf;
f.Color = [1 1 1];
f.Position = [-1050 45 800 700];
cMap = [rgb('Indigo'); ...
rgb('Orange'); ...
rgb('Green');
rgb('Cyan'); ...
];
plotAnnotData(end, end) = 4; % to keep the color mapping consistent
plot_2D_summary(plotAnnotData, FRAME_RATE, 'plotAxes', gca, 'colormap', cMap);
ax = gca();
ax.Title.String = {[groupName, ' behavior annotation'], ...
            'Indigo: quiescence,  cyan: locomotion,  orange: isolated movement'};
hold on;
xL = xlim();
for iBound = 1:numel(plotBlockBounds)
    currY = plotBlockBounds(iBound) - 0.5;
    plot(xL, [currY, currY], 'color', 'k', 'linewidth', 3);
end
xlim(xL);
draw_stim_lines(stimLines, {'red'}, 'FrameRate', FRAME_RATE)

if saveFig
    save_figure(f, saveDir, [saveDateStr, '_', regexprep(groupName, ' ', '_'), '_2D_annotation_summary']);
end


% FicTrac
smWin = 5;
speedData = plotFtData;
f = figure(202); clf;
ax = axes();
f.Color = [1 1 1];
f.Position = [-1860 45 800 700];
speedData(speedData > max(speedData(:) * 0.75)) = max(speedData(:) * 0.75);
plot_2D_summary(speedData, FRAME_RATE, 'plotAxes', ax);
colorbar
hold on;
xL = xlim();
for iBound = 1:numel(plotBlockBounds)
    currY = plotBlockBounds(iBound) - 0.5;
    plot(xL, [currY, currY], 'color', 'k', 'linewidth', 3);
end
xlim(xL);
ax.Title.String = [groupName, ' moveSpeed (mm/sec)'];
draw_stim_lines(stimLines, {'red'}, 'FrameRate', FRAME_RATE)

if saveFig
    save_figure(f, saveDir, [saveDateStr, '_', regexprep(groupName, ' ', '_'), '_2D_moveSpeed_summary']);
end

f = figure(203);clf
f.Color = [1 1 1];
f.Position = [-997 25 1000 950];

smMoveSpeed = plotFtData;
smMoveSpeedNoQui = smMoveSpeed; smMoveSpeedNoQui(plotAnnotData ~= 3) = nan;

% First plot with quiescence included
ax = subaxis(3,1,1, 'sv', 0.15);hold on
xData = (1/FRAME_RATE):(1/FRAME_RATE):(size(smMoveSpeed, 2)/FRAME_RATE);
plot(xData(3:end), mean(smMoveSpeed(:, 3:end), 'omitnan'), ...
    'linewidth', 2, 'color', 'k')
xlabel('Time (s)')
ylabel('Movement speed (mm/sec)')
title([groupName, '  —  Average movement speed (quiescence included)'])
ax.FontSize = 14;
ax.Title.FontSize = 12;
ax.XLabel.FontSize = 14;
plot_stim_shading(stimLines);
legend(ax.Children(1), 'Opto stim');
% xlim([xData(3) xData(end)]);

% Second plot excluding quiescence
ax = subaxis(3,1,2);hold on
plot(xData, mean(smMoveSpeedNoQui, 'omitnan'), ...
        'linewidth', 2, 'color', 'k')
xlabel('Time (s)')
ylabel('Movement speed (mm/sec)')
title([groupName, '  —  Average movement speed (quiescence excluded)'])
ax.FontSize = 14;
ax.Title.FontSize = 12;
ax.XLabel.FontSize = 14;
plot_stim_shading(stimLines);
legend(ax.Children(1), 'Opto stim');

% Third plot with locomotion probability
ax = subaxis(3,1,3);hold on
plot(xData(3:end), mean(plotAnnotData(:, 3:end) == 3, 'omitnan'), 'linewidth', 2, 'color', 'k')
xlabel('Time (s)')
ylabel('Locomotion probability')
title([groupName, '  — Locomotion probability'])
ax.FontSize = 14;
ax.Title.FontSize = 12;
ax.XLabel.FontSize = 14;
plot_stim_shading(stimLines);
legend(ax.Children(1), 'Opto stim');



if saveFig
    save_figure(f, saveDir, [saveDateStr, '_', regexprep(groupName, ' ', '_'), '_avg_moveSpeed']);
end