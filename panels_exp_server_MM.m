function panels_exp_server_MM()
% Gets information from Wombat about an experimental session and runs a triggered imaging 
% acqusition using the scanimage API. Modeled on the "fly_tracker_server" functions.


% Get hSI from the base workspace
hSI = evalin('base','hSI');             

% Start a server to listen on the socket to stim control computer (Wombat)
clear t;
PORT = 30000;
t = tcpip('0.0.0.0', PORT, 'NetworkRole', 'server');
set(t, 'InputBufferSize', 30000); 
set(t, 'TransferDelay', 'off');
disp(['About to listen on port: ' num2str(PORT)]);
fopen(t);
pause(1.0);
disp('Client connected');

while 1
    
    % Wait for something to be written to the server
    while (t.BytesAvailable == 0)
        pause(0.1);
    end
    
    % Read the message
    data = fscanf(t, '%s');
    data = strtrim(data);
    disp(data);

    % Break out of loop if the termination message is received
    if strcmp(data,'END_OF_SESSION')
        break;
    end
    
    % Parse message string (format should be: 'expID_expName_trial_0_dur_20')
    expID = regexp(data, '.*(?=_.*_trial)', 'match', 'once');
    expName = regexp(data, '(?<=_).*(?=_trial)', 'match', 'once');
    trialNum = str2double(regexp(data, '(?<=_.*_trial_).*(?=_dur)', 'match', 'once'));
    trialDuration = str2double(regexp(data, '(?<=_dur_).*', 'match', 'once'));
    
    % Make sure scanimage is saving to the correct directory
    parentDir = hSI.hScan2D.logFilePath;
    expDirName = [expID, '_', expName];
    if contains(parentDir, expDirName)
        % Scanimage path is already set to the correct directory
    else
        % Update ScanImage path
        expDir = fullfile(parentDir, expDirName);
        if ~isdir(expDir)
            % Experiment directory needs to be created
            mkdir(expDir)   
        end
        hSI.hScan2D.logFilePath = expDir;
    end
    
    % Set base file name
    hSI.hScan2D.logFileStem = [expID, '_', datestr(now, 'HHMMSS'), '_trial_', ...
            pad(num2str(trialNum), 3, 'left', '0')];
   
    % Set number of volumes to acquire
    volumeRate = hSI.hRoiManager.scanVolumeRate;
    disp(['Volume rate: ' num2str(volumeRate), ' VPS']);
    hSI.hFastZ.numVolumes = int32(ceil(trialDuration * volumeRate));
    disp(['Number of volumes per trial: ' num2str(hSI.hFastZ.numVolumes)]); 
    
    % Adjust other scanimage settings
    hSI.hScan2D.logFileCounter = 1;         % Set file counter to 1
    hSI.hChannels.loggingEnable = true;     % Enable logging
    hSI.hChannels.channelSave = 1;          % Make sure we're only logging the green PMT 
    hSI.extTrigEnable = true;               % Enable external trigger

    pause(1); % Can't remember if I had a good reason for putting this here, but it can't hurt
    
    % Start the grab
    hSI.startGrab();                        
    
    % Signal back to the client that it can start the trial
    fprintf(t, 'SI51_Acq_1');
end

% Clean up
hSI.extTrigEnable = false;         
hSI.hChannels.loggingEnable = false;

% Close the socket
fclose(t);
delete(t);

end