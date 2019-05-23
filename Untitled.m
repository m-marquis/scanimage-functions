parentDir = 'E:\Michael\Test';
baseFileName = 'cdata_20190523_194100_sid_0_bid_1_dur_50_nTrials_10_00001';

skipLastFile = 1;

% Figure out how many trials there are
allFiles = dir(fullfile(parentDir, [baseFileName, '*.pmt.dat']));
if skipLastFile
    nTrials = numel(allFiles) - 1;
else
    nTrials = numel(allFiles);
end

% Load data
metadataFileName = fullfile(parentDir,[baseFileName, '_00001']);
allPmtData = [];
allScannerPosData = [];
for iFile = 1:nTrials
   currBaseFileName = fullfile(parentDir, [baseFileName, '_', pad(num2str(iFile), 5, 'left', '0')]);
   [header, pmtData, scannerPosData, roiGroup] = scanimage.util.readLineScanDataFiles_MM(currBaseFileName, metadataFileName);
   allPmtData = cat(4, allPmtData, pmtData);
   allScannerPosData = cat(4, allScannerPosData, scannerPosData.G);
end
allPmtData = squeeze(allPmtData); % --> [samp, cycle, trial]

%% Separate data by ROIs

sampRate = header.sampleRate;

% Get handles for stimulus and control ROIs
allRois = roiGroup.rois;
scanRois = []; scanRoiNums = []; roiNames = [];
for iRoi = 1:numel(allRois)
    currRoiName = allRois(iRoi).scanfields.shortDescription;
    roiNames{iRoi} = currRoiName;
    roiDurs(iRoi) = allRois(iRoi).scanfields.duration;
    roiDursSamples(iRoi) = floor(roiDurs(iRoi) * sampRate); 
    if ~strcmp(currRoiName(7:end), 'pause') && ...
                ~strcmp(currRoiName(7:end), 'park')
        scanRoiNums(end + 1) = iRoi;
    end
end
scanRois = allRois(scanRoiNums);
hStimRoi = scanRois(1).scanfields;
hControlRoi = scanRois(2).scanfields;
hImageRois = scanRois(3:end);

% Separate PMT data
nSamples = size(allPmtData, 1);
roiDataSep = [];  roiDataAvg = [];
for iRoi = 1:numel(roiDurs)
   
        if iRoi == 1
            startSample = 1;
        else
            startSample = sum(roiDursSamples(1:iRoi-1));
        end
        if iRoi == numel(roiDurs)
            endSample = nSamples;
        else
            endSample = sum(roiDursSamples(1:iRoi)); 
        end
        
        roiDataSep{iRoi} = allPmtData(startSample:endSample, :, :);     % --> [samp, cycle, trial]
        roiDataAvg(:, :, iRoi) = squeeze(mean(roiDataSep{iRoi}, 1));    % --> [cycle, trial, ROI]
end

% Separate actual ROIs from pauses and parks
stimRoiData = roiDataAvg(:,:,scanRoiNums(1));       % --> [cycle, trial]
ctrlRoiData = roiDataAvg(:,:,scanRoiNums(2));       % --> [cycle, trial]
imgRoiData = roiDataAvg(:,:,scanRoiNums(3:end));    % --> [cycle, trial, ROI]


%%

currTrial = 2;
figure; hold on
nTrials = size(stimRoiData, 2);
cm = jet(nTrials);
xData = header.frameDuration:header.frameDuration:(header.frameDuration * header.numFrames);
for iTrial = 1:nTrials
plot(xData', stimRoiData(:, iTrial), 'Color', cm(iTrial, :));
plot(xData', ctrlRoiData(:, iTrial), 'Color', cm(iTrial, :));
end
% figure;
% plot([xData', xData'], squeeze(imgRoiData(:,currTrial, :)));

%% For each trial, identify the frames when the laser stim turned on and off

















