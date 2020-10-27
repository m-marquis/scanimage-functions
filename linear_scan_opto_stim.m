function linear_scan_opto_stim(src, event, varargin)
% Custom user function to allow online control of laser power in multiple linear scan ROIs at specific times in
% each trial by counting the acquired frames and updating ROI power accordingly.
%
% First argument is a 2-element vector specifying [startTime, endTime] in seconds.
%
% Second argument specifies what laser power % to use for the photostimulation (laser power at the 
%       beginning of the acquisition will be used as the imaging laser power).
%  
% Third argument specifies the number of photostim ROIs (defaults to one if you omit this argument). 
% For each stim ROI, the ROI immediately after that is assumed to be the cognate control ROI. For
% example, if you pass [2] as this argument, ROIs 1 and 3 will be photostims and ROIs 2 and 4 will 
% be the corresponding control ROIs.

hSI = src.hSI;

switch event.EventName
    
    case 'acqModeStart'
        try
            disp('Starting grab')
            
            % Get handles for stimulus and control ROIs
            allRois = hSI.hRoiManager.roiGroupLineScan;
            scanRois = []; scanRoiNums = [];
            for iRoi = 1:numel(allRois.rois)
                currRoiName = allRois.rois(iRoi).scanfields.shortDescription;
                if ~strcmp(currRoiName(7:end), 'pause') && ...
                            ~strcmp(currRoiName(7:end), 'park')
                    scanRoiNums(end + 1) = iRoi;
                end
            end
            scanRois = allRois.rois(scanRoiNums);
            
             if numel(varargin) > 2
                nStimRois = varargin{3};
            else
                nStimRois = 1;
             end
            hSI.extCustomProps.nStimRois = nStimRois;
            for iStimRoi = 1:nStimRois
                currRoiNum = (2 * iStimRoi) - 1;
                hStimRois(iStimRoi) = scanRois(currRoiNum).scanfields;
                hControlRois(iStimRoi) = scanRois(currRoiNum + 1).scanfields;
            end
            hImageRois = scanRois((2 * nStimRois) + 1:end);
            
            % Setup laser power as acquisition starts
            hSI.extCustomProps.stimROIPower = varargin{2}; % Get stim power
            hSI.extCustomProps.nFramesAcq = 0;
            hSI.hBeams.beamsOff();
            hSI.extCustomProps.stimOn = 0;
            for iStimRoi = 1:nStimRois
                hStimRois(iStimRoi).powers = 0.1;
                hControlRois(iStimRoi).powers = hSI.extCustomProps.stimROIPower;
            end
            for iRoi = 1:numel(hImageRois)
                hImageRois(iRoi).scanfields.powers = hSI.extCustomProps.imagingPower;
            end
            hSI.hBeams.updateBeamBufferAsync(true); % Important to avoid ScanImage bug
            
            % Calculate trial start timing
            stimTimes = varargin{1}; % [startTime, endTime]
            cps = hSI.hRoiManager.scanFrameRate;
            cpt = hSI.extCustomProps.cyclesPerTrial;
            trialStartCycles = [];
            for iTrial = 1:hSI.extCustomProps.nTrials
                if iTrial == 1
                    trialStartCycles(iTrial) = 1;
                else
                    trialStartCycles(iTrial) = 1 + (cpt * (iTrial - 1));
                end
            end
            hSI.extCustomProps.trialStartCycles = trialStartCycles;
            hSI.extCustomProps.pendingTrialStartCycles = trialStartCycles;
            hSI.extCustomProps.nextTrialNum = 1;
            disp(['trialStartCycles: ', num2str(trialStartCycles)]);
            
            
            % Calculate stim timing
            if ~isempty(stimTimes)
                stimStartCycles = [];
                stimEndCycles = [];
                relStimStartCycle = ceil(stimTimes(1) * cps) - 1;
                relStimEndCycle = ceil(stimTimes(2) * cps) - 1;
                if relStimStartCycle < 5
                    relStimStartCycle = 5; % Leaving some padding for safety
                    disp('WARNING: stim start cycle too close to beginning of trial!')
                end
                if relStimEndCycle >= cpt - 5
                    relStimEndCycle = cpt - 6; % Leaving some padding for safety
                    disp('WARNING: stim end cycle too close to end of trial!')
                end
                for iTrial = 1:hSI.extCustomProps.nTrials
                    stimStartCycles(iTrial) = trialStartCycles(iTrial) + relStimStartCycle;
                    stimEndCycles(iTrial) = trialStartCycles(iTrial) + relStimEndCycle;
                end
                disp(['stimStartCycles: ', num2str(stimStartCycles)])
                disp(['stimEndCycles: ', num2str(stimEndCycles)])
                hSI.extCustomProps.stimStartCycles = stimStartCycles;
                hSI.extCustomProps.stimEndCycles = stimEndCycles;
                hSI.extCustomProps.pendingStimStartCycles = stimStartCycles;
                hSI.extCustomProps.pendingStimEndCycles = stimEndCycles;
                hSI.extCustomProps.pendingStimStartCycles(end + 1) = ...
                        cpt * (hSI.extCustomProps.nTrials + 1); % So it doesn't break on the last trial
                hSI.extCustomProps.pendingStimEndCycles(end + 1) = ...
                        cpt * (hSI.extCustomProps.nTrials + 1); % So it doesn't break on the last trial
                disp(['Photostim at ', num2str(hSI.extCustomProps.stimROIPower), '% power'])
                disp(['Stim on time: ', num2str(stimTimes(1)), ' sec'])
                disp(['Stim off time: ', num2str(stimTimes(2)), ' sec'])
                disp(['Stim on from cycles ', num2str(relStimStartCycle), ' to ', ...
                        num2str(relStimEndCycle), ' of each trial'])
            else
                hSI.extCustomProps.stimStartCycles = [];
                hSI.extCustomProps.stimEndCycles = [];
                hSI.extCustomProps.pendingStimStartCycles = [];
                hSI.extCustomProps.pendingStimEndCycles = [];
                disp('No photostimulation in this block');
            end
            disp('--------------------------------------------')
            %         disp('Starting trial 1...')
        catch ME
            disp(ME.message);
            rethrow(ME)
        end
    case 'acqModeDone'
        
        disp('Block finished')
        
        % Save log file for the block that just finished
        optoStimInfo = hSI.extCustomProps;        
        baseName = hSI.hScan2D.logFileStem;
        saveDir = hSI.hScan2D.logFilePath;
        save(fullfile(saveDir, [baseName, '.mat']), 'optoStimInfo');
        
        % Clear extCustomProps
        hSI.hBeams.powers = hSI.extCustomProps.imagingPower;
        hSI.extCustomProps = [];
        
    case 'frameAcquired'
        try
           
            % Get handles for stimulus and control ROIs
            allRois = hSI.hRoiManager.roiGroupLineScan.rois;
            scanRois = []; scanRoiNums = [];
            for iRoi = 1:numel(allRois)
                currRoiName = allRois(iRoi).scanfields.shortDescription;
                if ~strcmp(currRoiName(7:end), 'pause') && ...
                            ~strcmp(currRoiName(7:end), 'park')
                    scanRoiNums(end + 1) = iRoi;
                end
            end
            scanRois = allRois(scanRoiNums);
            nStimRois = hSI.extCustomProps.nStimRois;
            for iStimRoi = 1:nStimRois
                currRoiNum = (2 * iStimRoi) - 1;
                hStimRois(iStimRoi) = scanRois(currRoiNum).scanfields;
                hControlRois(iStimRoi) = scanRois(currRoiNum + 1).scanfields;
            end
            
            % Get current internal frame count
            currCycleCount = hSI.hScan2D.hAcq.frameCounter;
            
            % Modify laser power if necessary
            nextStimStartCycle = hSI.extCustomProps.pendingStimStartCycles(1);
            nextStimEndCycle = hSI.extCustomProps.pendingStimEndCycles(1);
            if ~isempty(nextStimStartCycle)
                if currCycleCount >= nextStimEndCycle
                    
                    % Switch laser to control ROI
                    hSI.extCustomProps.pendingStimEndCycles = ...
                            hSI.extCustomProps.pendingStimEndCycles(2:end);
                    for iStimRoi = 1:nStimRois
                        hStimRois(iStimRoi).powers = 0.3;
                        hControlRois(iStimRoi).powers = hSI.extCustomProps.stimROIPower;
                    end
                    hSI.hBeams.updateBeamBufferAsync(true);
                    disp(['Setting laser to ', num2str(hSI.extCustomProps.stimROIPower), ...
                        '% power in control ROIs'])
                    
                elseif currCycleCount >= nextStimStartCycle
                    
                    % Switch laser to stim ROI
                    hSI.extCustomProps.pendingStimStartCycles = ...
                        hSI.extCustomProps.pendingStimStartCycles(2:end);
                    for iStimRoi = 1:nStimRois
                        hStimRois(iStimRoi).powers = hSI.extCustomProps.stimROIPower;
                        hControlRois(iStimRoi).powers = 0.3;
                    end
                    hSI.hBeams.updateBeamBufferAsync(true);
                    disp(['Setting laser to ', num2str(hSI.extCustomProps.stimROIPower), ...
                            '% power in stim ROIs'])                    
                end
                
            end
            
            % Notify user if starting new trial
            if hSI.extCustomProps.nextTrialNum <= hSI.extCustomProps.nTrials
                nextTrialStartCycle = hSI.extCustomProps.trialStartCycles(hSI.extCustomProps.nextTrialNum);
                if currCycleCount >= nextTrialStartCycle
                    disp(['Starting trial ', num2str(hSI.extCustomProps.nextTrialNum), '...']);
                    hSI.extCustomProps.nextTrialNum = hSI.extCustomProps.nextTrialNum + 1;
                end
            end
            
        catch ME
            disp(ME.message);
            rethrow(ME)
        end

end%case


end