

% parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_05_25_exp_1';
parentDir = 'E:\Michael\2019_05_25_exp_1';

%% PROCESS RAW DATA FOR ALL BLOCKS
    
metadataFiles = dir(fullfile(parentDir, '*_00001.meta.txt'));
nBlocks = numel(metadataFiles);

for iBlock = 1:nBlocks
try   
    % Get current file names
    currMetaFileName = metadataFiles(iBlock).name;
    baseFileName = regexp(currMetaFileName, '.*(?=_00001.meta.txt)', 'match');
    baseFileName = baseFileName{:};
    refImgFileName = [baseFileName, '_00001.ref.dat'];
    
    blockNum = str2double(regexp(baseFileName, '(?<=bid_).*(?=_dur)', 'match'));
    disp(['Processing block ', num2str(blockNum)])  
    
    % Figure out how many files the data was split into for this block
    currDataFiles = dir(fullfile(parentDir, [baseFileName, '*.pmt.dat']));
    nDataFiles = numel(currDataFiles);
    
    % Load metadata and roiGroup data for this block
    blockFileStem = fullfile(parentDir, regexprep(currMetaFileName, '.meta.txt', ''));
    [siData, ~, ~, roiGroup] = scanimage.util.readLineScanDataFiles_MM(blockFileStem);
    
    % Fix reference image file extension and load that data as well
    newRefImgFileName = regexprep(refImgFileName, '.ref.dat', '.mat');
    copyfile(fullfile(parentDir, refImgFileName), fullfile(parentDir, newRefImgFileName));
    refImgData = load(fullfile(parentDir, newRefImgFileName));        
    
    % Process ref image data into a more convenient form
    roiPlanes = ismember(refImgData.contextImageZs{end}, roiGroup.zs);
    refImgCells = [refImgData.contextImageImgs{end}{roiPlanes}];
    refImgStack = [];
    for iPlane = 1:numel(refImgCells)
        refImgStack = cat(3, refImgStack, refImgCells{iPlane}{1}'); % --> [y, x, slice]
    end
    refImgCP = refImgData.contextImageRoiCPs{end}{1}{1};
    refImgZs = refImgData.contextImageZs{end}(roiPlanes);
    
    % Parse ROI data
    sampRate = siData.sampleRate;
    allRois = roiGroup.rois;
    roiMetadata = [];
    scanRoiNums = [];
    nRois = numel(allRois);
    for iRoi = 1:nRois
        currRoi = allRois(iRoi);
        currRoiName = currRoi.scanfields.shortDescription;
        currRoiDur = currRoi.scanfields.duration;
        
        roiMetadata.allRois(iRoi).name = currRoiName;
        roiMetadata.allRois(iRoi).zDepth = currRoi.zs;
        roiMetadata.allRois(iRoi).duration = round(currRoiDur, 4); % To prevent floating point issue
        roiMetadata.allRois(iRoi).durationInSamples = floor(round(currRoiDur, 4) * sampRate);
        roiMetadata.allRois(iRoi).centerX = currRoi.scanfields.centerXY(1);
        roiMetadata.allRois(iRoi).centerY = currRoi.scanfields.centerXY(2);
        roiMetadata.allRois(iRoi).sizeX = currRoi.scanfields.sizeXY(1);
        roiMetadata.allRois(iRoi).sizeY = currRoi.scanfields.sizeXY(2);
        roiMetadata.allRois(iRoi).stimParams = currRoi.scanfields.stimparams;
        roiMetadata.allRois(iRoi).transformParams = currRoi.scanfields.transformParams;
        
        if ~strcmp(currRoiName(7:end), 'pause') && ...
                ~strcmp(currRoiName(7:end), 'park')
            scanRoiNums(end + 1) = iRoi;
        end
    end% iRoi
    roiMetadata.scanRoiNums = scanRoiNums;
    roiMetadata.scanRois = roiMetadata.allRois(scanRoiNums);
    roiMetadata.refImgStack = refImgStack;
    roiMetadata.refimgCP = refImgCP;
    roiMetadata.refImgZs = refImgZs;
    
    % Extract fluorescence data averaged across each ROI
    roiDataAvg = []; cycleCounts = [];
    for iFile = 1:nDataFiles
        
        % Load data for current file
        disp(['Loading file ', num2str(iFile)])
        currBaseFileName = fullfile(parentDir, [baseFileName, '_', ...
                pad(num2str(iFile), 5, 'left', '0')]);
        [~, pmtData, ~, ~] = scanimage.util.readLineScanDataFiles_MM(currBaseFileName, ... 
                fullfile(parentDir, regexprep(currMetaFileName, '.meta.txt', '')));
        
        cycleCounts(end + 1) = size(pmtData, 3);
        
        % Separate PMT data from current file
        nSamples = siData.samplesPerFrame;
        currRoiDataAvg = [];
        for iRoi = 1:nRois
            if iRoi == 1
                startSample = 1;
            else
                startSample = sum([roiMetadata.allRois(1:iRoi-1).durationInSamples]);
            end
            if iRoi == numel([roiMetadata.allRois.duration])
                endSample = nSamples;
            else
                endSample = sum([roiMetadata.allRois(1:iRoi).durationInSamples]);
            end
            currRoiDataAvg(:, iRoi) = squeeze(mean(pmtData(startSample:endSample, :, :), 1)); % --> [cycle, ROI]
        end
        roiDataAvg = cat(1, roiDataAvg, currRoiDataAvg); % --> [cycle, ROI]
        
    end%iFile
    
    % Separate actual ROIs from pauses and parks
    stimRoiData = roiDataAvg(:,scanRoiNums(1));       % --> [cycle]
    ctrlRoiData = roiDataAvg(:,scanRoiNums(2));       % --> [cycle]
    imgCtrlRoiData = roiDataAvg(:,scanRoiNums(3));    % --> [cycle]
    imgRoiData = roiDataAvg(:,scanRoiNums(4:end));    % --> [cycle, ROI]
    
    nCyclesTotal = numel(stimRoiData);
    disp(['Total cycles = ' num2str(nCyclesTotal)])
    disp(['Cycle counts = ', num2str(cycleCounts)])
    
    % Save data for easy access
    save(fullfile(parentDir, [baseFileName, '_SI_data']), 'siData', 'nDataFiles', 'cycleCounts', ...
        'roiDataAvg', 'roiMetadata', 'blockNum');
    
catch ME
    disp(['Warning: failed to process data for block ', num2str(blockNum)])
    disp(['with error message: ' ME.message]);
    
end    

end%iBlock

%% Load any high-res stacks used for ROI creation and save averaged versions
tifs = dir(fullfile(parentDir, '*stack*.tif'));
for iFile = 1:numel(tifs)
    disp(['Reading tif stack #', num2str(iFile)])
   [header, aout] = opentif(fullfile(parentDir, tifs(iFile).name));
   tifData = squeeze(mean(aout, 6)); % --> [y, x, slice]
   saveFile = fullfile(parentDir, ['avg_', tifs(iFile).name, 'f']);
   try
   saveastiff(uint32(tifData), saveFile);
   catch
   end
   
   % Save a .mat version too, with metadata
   save(fullfile(regexprep(saveFile, '.tiff', '.mat')), 'tifData', 'header', '-v7.3')
end
