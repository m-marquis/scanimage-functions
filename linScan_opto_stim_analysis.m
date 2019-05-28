

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_05_24_exp_1';


%% PROCESS RAW DATA

try
fileName = uigetfile(fullfile(parentDir, '*_00001.meta.txt'));
baseFileName = regexp(fileName, '.*(?=_00001.meta.txt)', 'match');
baseFileName = baseFileName{:};

% Figure out how many trials there are
allFiles = dir(fullfile(parentDir, [baseFileName, '*.pmt.dat']));
nFiles = numel(allFiles);

% Load data
metadataFileName = fullfile(parentDir,[baseFileName, '_00001']);

% Get metadata from first trial
[header, ~, ~, roiGroup] = scanimage.util.readLineScanDataFiles_MM(metadataFileName);

% Parse ROI data
sampRate = header.sampleRate;
allRois = roiGroup.rois;
roiMetadata = [];
scanRoiNums = [];
nRois = numel(allRois);
for iRoi = 1:nRois
    currRoi = allRois(iRoi);
    currRoiName = currRoi.scanfields.shortDescription;
    currRoiDur = currRoi.scanfields.duration;
    
    roiMetadata.name{iRoi} = currRoiName;
    roiMetadata.zDepth(iRoi) = currRoi.zs;
    roiMetadata.duration(iRoi) = currRoiDur;
    roiMetadata.durationInSamples(iRoi) = floor(currRoiDur * sampRate);
    roiMetadata.centerXY(iRoi, :) = currRoi.scanfields.centerXY;
    roiMetadata.sizeXY(iRoi, :) = currRoi.scanfields.sizeXY;
    roiMetadata.stimParams{iRoi} = currRoi.scanfields.stimparams;
    roiMetadata.transformParams(iRoi) = currRoi.scanfields.transformParams;
    
    if ~strcmp(currRoiName(7:end), 'pause') && ...
                ~strcmp(currRoiName(7:end), 'park')
        scanRoiNums(end + 1) = iRoi;
    end
end
roiMetadata.scanRoiNums = scanRoiNums;
scanRois = allRois(scanRoiNums);

% Extract fluorescence data averaged across each ROI
roiDataAvg = []; cycleCounts = [];
for iFile = 1:nFiles
    
    % Load data for current file
    disp(['Loading file ', num2str(iFile)])
    currBaseFileName = fullfile(parentDir, [baseFileName, '_', pad(num2str(iFile), 5, 'left', '0')]);
    [~, pmtData, ~, ~] = scanimage.util.readLineScanDataFiles_MM(currBaseFileName, metadataFileName);
    
    cycleCounts(end + 1) = size(pmtData, 3);
    
    % Separate PMT data from current file
    nSamples = header.samplesPerFrame;
    currRoiDataAvg = [];
    for iRoi = 1:nRois
        if iRoi == 1
            startSample = 1;
        else
            startSample = sum(roiMetadata.durationInSamples(1:iRoi-1));
        end
        if iRoi == numel(roiMetadata.duration)
            endSample = nSamples;
        else
            endSample = sum(roiMetadata.durationInSamples(1:iRoi));
        end
        currRoiDataAvg(:, iRoi) = squeeze(mean(pmtData(startSample:endSample, :, :), 1)); % --> [cycle, ROI]
    end
    roiDataAvg = cat(1, roiDataAvg, currRoiDataAvg); % --> [cycle, ROI]    
end

% Separate actual ROIs from pauses and parks
stimRoiData = roiDataAvg(:,scanRoiNums(1));       % --> [cycle]
ctrlRoiData = roiDataAvg(:,scanRoiNums(2));       % --> [cycle]
imgCtrlRoiData = roiDataAvg(:,scanRoiNums(3));    % --> [cycle]
imgRoiData = roiDataAvg(:,scanRoiNums(4:end));    % --> [cycle, ROI]

nCyclesTotal = numel(stimRoiData);
disp(['Total cycles = ' num2str(nCyclesTotal)])
disp(['Cycle counts = ', num2str(cycleCounts)])

% Save data for easy access
save(fullfile(parentDir, [baseFileName, '_SI_data']), 'header', 'nFiles', 'cycleCounts', ...
        'roiDataAvg', 'roiMetadata');
         
catch ME
    rethrow(ME); 
end

%% Load any high-res stacks used for ROI creation and save averaged versions
tifs = dir(fullfile(parentDir, '*stack*.tif'));
for iFile = 1:numel(tifs)
    disp(['Reading tif stack #', num2str(iFile)])
   [header, aout] = opentif(fullfile(parentDir, tifs(iFile).name));
   tifData = squeeze(mean(aout, 6)); % --> [y, x, slice]
   saveFile = fullfile(parentDir, ['avg_', tifs(iFile).name, 'f']);
   saveastiff(uint32(tifData), saveFile);
end

%% LOAD DATA FROM EXISTING FILE

[fileName, pathName] = uigetfile(fullfile(parentDir, '*_SI_data.mat'));
load(fullfile(pathName, fileName), 'cycleCounts', 'nFiles', 'header', 'roiDataAvg', ...
        'roiMetadata');

stimRoiData = roiDataAvg(:, roiMetadata.scanRoiNums(1));
ctrlRoiData = roiDataAvg(:, roiMetadata.scanRoiNums(2));

%% Plot fluorescence from photostim/control ROIs 

f = figure(1);clf; hold on
nCyclesTotal = size(roiDataAvg, 1);
xData = header.frameDuration:header.frameDuration:(header.frameDuration * nCyclesTotal);
plot(stimRoiData, 'Color', 'r');
plot(ctrlRoiData, 'Color', 'b');
legend('Photostim', 'Control', 'autoupdate', 'off')

%% FIND STIM ON/OFF CYCLES

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






























