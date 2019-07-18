function [] = fly_tracker_server_MM_chrimson()
% fly_tracker_server 
% example for using the ScanImage API to set up a grab
hSI = evalin('base','hSI');             % get hSI from the base workspace

%%%%
% Start a server to listen on the socket to fly tracker
%%%%
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
    
    % Read the message     
    while (t.BytesAvailable == 0)
        pause(0.1);
    end
        
    data = fscanf(t, '%s');
    data = strtrim(data);
    disp(data);

    if( strcmp(data,'END_OF_SESSION') == 1 )
        break;
    end
    
    % Extract trial info
    durCell = regexp(data, '(?<=dur_).*(?=_nTrials)', 'match');
    blockDur = str2num(durCell{:});
    nTrialCell = regexp(data, '(?<=nTrials_).*', 'match');
    nTrials = str2num(nTrialCell{:});
    hSI.hScan2D.logFileStem = data; % set the base file name for the Tiff file     
    
     % Configure acquitions
    hSI.hChannels.loggingEnable = true;     % Make sure logging is enabled  
    hSI.hScan2D.logFileCounter = 1;         % Set file counter to 1    
    hSI.extTrigEnable = true;               % Enable external trigger 
    hSI.extCustomProps.nTrials = nTrials;   % Set total number of trials
    hSI.hFastZ.enable = true;
    hSI.acqsPerLoop = 1;
    
    % Set up scan cycles
    cycles_per_second = hSI.hRoiManager.scanFrameRate;
    disp(['Cycles/sec: ' num2str( cycles_per_second )]);
    hSI.hStackManager.framesPerSlice = int32(ceil((blockDur) * cycles_per_second) + 1);
    disp(['Total number of cycles in block: ' num2str( hSI.hStackManager.framesPerSlice)]);
    hSI.extCustomProps.cyclesPerTrial = floor(hSI.hStackManager.framesPerSlice ./ nTrials);
    disp(['Cycles per trial: ', num2str(hSI.extCustomProps.cyclesPerTrial)]);
    hSI.hScan2D.logFramesPerFile = hSI.extCustomProps.cyclesPerTrial;
    
    % Make sure user functions are enabled for laser power modulation
    hSI.extCustomProps.imagingPower = hSI.hBeams.powers;
    for iFun = 1:numel(hSI.hUserFunctions.userFunctionsCfg)
        hSI.hUserFunctions.userFunctionsCfg(iFun).Enable = true;
    end
    
    % Start acquisitions
    pause(1.0);
    hSI.startGrab();
    disp('Grab started');
    % Signal back to the fly tracker client that it can start daq and image
    % acquisition.
    fprintf(t, 'SI51_Acq_1');
    disp('Printed ready signal to server');
end

% Clean up
hSI.extTrigEnable = false;

%close the socket
fclose(t);
delete(t);

end