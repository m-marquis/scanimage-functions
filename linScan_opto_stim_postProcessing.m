expDate = '2019_08_02_exp_1';
sid = 0;

parentDir = fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate);
% parentDir = 'E:\Michael\2019_05_25_exp_1';

saveDir = fullfile('D:\Dropbox (HMS)\2P Data\Analysis', expDate, ['sid_', num2str(sid)]);
saveDateStr = regexprep(expDate, 'exp_.*', '');
saveDateStr = regexprep(saveDateStr, '\_', '');
saveDateStr = [saveDateStr, '-', expDate(end)];

%% LOAD DATA FROM ALL BLOCKS

clear roiNames;

roiNames(1) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(2) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(3) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(4) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(5) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(6) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(7) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(8) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(9) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(10) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
roiNames(11) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
% roiNames(12) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
% roiNames(13) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
% roiNames(14) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
% roiNames(15) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
% roiNames(16) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};
% roiNames(17) = {["stim-1", "stimCtrl-1", "stim-2", "stimCtrl-2", "Img"]};

odorStimRelTimes = repmat({[]}, 1, numel(roiNames));

try
    
if exist(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'), 'file')
    load(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'))
    disp('Loaded existing block data file')
else
    blockDataFiles = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*_SI_data.mat']));
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
      
    save(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'), 'allBlockData');
end



catch foldME; rethrow(foldME); end


%% LOAD BEHAVIOR DATA

saveFig = 1;

annotFileName = 'autoAnnotations.mat';
annotData = load(fullfile(parentDir, ['sid_', num2str(sid)], annotFileName));
FRAME_RATE = annotData.frameInfo.FRAME_RATE;

% blockBounds = [1:20:size(annotData.flowArr, 2)-1, size(annotData.flowArr, 2)]; 
blockBounds = [1:10:size(annotData.flowArr, 2)-1, size(annotData.flowArr, 2)]; 
% blockBounds = [1, 21, 41, 61, 71, 81, 91, 101, 111, 121, 130]; 
% blockBounds = [1, 6, 11, 16, 21, 26, 46, 66, 86, 105]; 
% blockBounds = [1:10:41 46 51:10:81 86 90] + 15;
% blockBounds = [1:10:81 101:10:141 150]

try
    
if ~isfield(annotData, 'blockBounds')
    annotData.blockBounds = blockBounds;
    save(fullfile(parentDir, ['sid_', num2str(sid)], annotFileName), '-struct', 'annotData');
    disp('Saving new block boundaries...')
else
    disp('Using previously saved block boundaries') 
    blockBounds = annotData.blockBounds;
end
    
    
    
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


% Plot 2D summaries of annotation and FicTrac data throughout the entire experiment

% Annotations
f = figure(1); clf;
f.Color = [1 1 1];
f.Position = [-1050 45 800 700];
cMap = [rgb('Indigo'); ...
rgb('Orange'); ...
rgb('Green');
rgb('Cyan'); ...
];
allAnnotData = annotData.trialAnnotations;
allAnnotData(end, end) = 4; % to keep the color mapping consistent
plot_2D_summary(allAnnotData, FRAME_RATE, 'plotAxes', gca, 'colormap', cMap);
ax = gca();
ax.Title.String = {'Full experiment behavior annotation', ...
            'Indigo: quiescence,  cyan: locomotion,  orange: isolated movement'};
hold on;
xL = xlim();
for iBound = 1:numel(blockBounds)
    currY = blockBounds(iBound) - 0.5;
    plot(xL, [currY, currY], 'color', 'k', 'linewidth', 3);
end
xlim(xL);

if saveFig
    save_figure(f, saveDir, [saveDateStr, '_AllTrials_raw_2D_annotation_summary']);
end


% FicTrac
smWin = 5;
speedData = smoothdata(annotData.ftData.moveSpeed' .* 4.5 .* FRAME_RATE, ...
        2, 'gaussian', smWin); % --> [frame, stim]
f = figure(2); clf;
ax = axes();
f.Color = [1 1 1];
f.Position = [-1860 45 800 700];
speedData(speedData > max(speedData(:) * 0.75)) = max(speedData(:) * 0.75);
plot_2D_summary(speedData, FRAME_RATE, 'plotAxes', ax);
colorbar
hold on;
xL = xlim();
for iBound = 1:numel(blockBounds)
    currY = blockBounds(iBound) - 0.5;
    plot(xL, [currY, currY], 'color', 'k', 'linewidth', 3);
end
xlim(xL);
ax.Title.String = 'Full experiment moveSpeed (mm/sec)';

if saveFig
    save_figure(f, saveDir, [saveDateStr, '_AllTrials_raw_2D_moveSpeed_summary']);
end






catch foldME; rethrow(foldME); end

%% IDENTIFY STIM ON- and OFFSETS

currBlock = 10;
manualThresh = 15;


try
currBlockData = allBlockData([allBlockData.blockNum] == currBlock);
scanRois = currBlockData.roiMetadata.scanRoiNums;
roiDataAvg = currBlockData.roiDataAvg - min(currBlockData.roiDataAvg(:)); % Setting min to zero

% roiDataAvg = roiDataAvg - 775; roiDataAvg(roiDataAvg < 0) = 0;

stimRoiData = roiDataAvg(:, scanRois(1));
ctrlRoiData = roiDataAvg(:, scanRois(2));
siData = currBlockData.siData;
cycleRate = currBlockData.cycleRate;

stimPower = currBlockData.siData.SI.hUserFunctions.userFunctionsCfg(1).Arguments{2};
imgPower = currBlockData.siData.SI.hBeams.powers;

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
stimCycleDurs = stimOffCycles - stimOnCycles;
interStimDurs = [stimOnCycles(1), stimOnCycles(2:end) - stimOffCycles(1:end-1), ...
        nCyclesTotal - stimOffCycles(end)];

disp(' ')
disp(['stimCycleDurs =  ', num2str(stimCycleDurs)])
disp(['interStimDurs =  ', num2str(interStimDurs)])

catch foldME; rethrow(foldME); end
%% DIVIDE DATA INTO INDIVIDUAL STIM EPOCHS

% Check stim durations
disp(' ')
disp(['stimCycleDurs = ', num2str(stimCycleDurs)])
disp(['interStimDurs = ', num2str(interStimDurs)])
    
skipCycles = [];
analysisWindow = [47];
targetStimDur = 23;
smWin = 3;

try
if ~isempty(skipCycles)
   stimOnCycles(skipCycles) = [];
   stimOffCycles(skipCycles) = [];
   stimCycleDurs(skipCycles) = [];
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

stimStart = analysisWindowFrames(1);
stimEnd = stimStart + cyc2frame.convert(targetStimDur);

allStimFlData = permute(allStimFlData, [1 3 2]); % --> [cycle, stim, roi]
allStimFtData = permute(allStimFtData, [1 3 2]); % --> [frame, stim, var]

% Update main data structure with this block's info
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.flData = allStimFlData;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.FtData = allStimFtData;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.annotData = allStimAnnotData;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.stimStart = stimStart;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.stimEnd = stimEnd;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.stimOnCycles = stimOnCycles;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.stimOffCycles = stimOffCycles;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.analysisWindow = analysisWindow;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.targetStimDur = targetStimDur;
allBlockData([allBlockData.blockNum] == currBlock).stimSepData.skipCycles = skipCycles;

disp('Data separated into stim epochs')

catch foldME; rethrow(foldME); end

%% SAVE PROCESSED BLOCK DATA

save(fullfile(parentDir, ['sid_', num2str(sid)], 'allBlockData.mat'), 'allBlockData')

