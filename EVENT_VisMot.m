%% EVENTS FILE FOR TASK VISMOT

function EVENT_VisMot
    %Define path for input and output dir
    inputDir = 'Y:\projekty\ELICIT\derivatives\spm_processed_dataset\EVENTS\VisMot';       
    outputDir = 'Y:\projekty\ELICIT\derivatives\spm_processed_dataset\EVENTS\VisMot'; 

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %Find all logs which need to reformat
    logFiles = dir(fullfile(inputDir, 'ELICIT_Visual_motor_checkerboard-*.txt'));

    %For loop across all log files
    for f = 1:length(logFiles)
        logPath = fullfile(logFiles(f).folder, logFiles(f).name);

        %Find ID of subject
        subjID = regexp(logFiles(f).name, 'ELICIT_Visual_motor_checkerboard-(\d+)-\d+', 'tokens');
        subjID = subjID{1}{1};

        %Load file
        fid = fopen(logPath, 'r');
        rawLines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        lines = rawLines{1};

        %Find "PulseTime" to find out TR 
        pulseTimes = [];
        for i = 1:length(lines)
            if contains(lines{i}, 'PulseTime_')
                val = str2double(extractAfter(lines{i}, ': '));
                if ~isnan(val)
                    pulseTimes(end+1) = val;
                end
            end
        end

        if numel(pulseTimes) < 2
            warning('Two PulseTime are required in file %s.', logFiles(f).name);
            continue;
        end

        TR = pulseTimes(2) - pulseTimes(1);
        if abs(TR - 1800) <= 5
            trType = 'tr1800';
        elseif abs(TR - 800) <= 5
            trType = 'tr800';
        else
            trType = sprintf('tr%d', round(TR));
            warning('Unknown TR');
        end

        pulseZero = pulseTimes(1);

        %Define the structure
        events = struct('onset', {}, 'duration', {}, 'trial_type', {}, 'stim_type', {}, 'button_press', {});

        %Work with logs
        i = 1;
        while i <= length(lines)
            if contains(lines{i}, '*** LogFrame Start ***')
                blockName = '';
                blockOnset = NaN;
                keyTimes = [];
                keyNames = [];

                j = i + 1;
                while j <= length(lines) && ~contains(lines{j}, '*** LogFrame End ***')
                    line = strtrim(lines{j});
                    if startsWith(line, 'Block: ')
                        blockName = strtrim(extractAfter(line, 'Block: '));
                    elseif startsWith(line, 'BlockOnset:')
                        blockOnset = str2double(extractAfter(line, 'BlockOnset:'));
                    elseif startsWith(line, 'Time_')
                        tval = str2double(extractAfter(line, ': '));
                    if ~isnan(tval)
                        keyTimes(end+1) = tval;
                    end
                    elseif startsWith(line, 'Key_')
                        kval = strtrim(extractAfter(line, ': '));
                        keyNames{end+1} = kval;
                    end
                    j = j + 1;
                end

                %Rest block
                if strcmpi(blockName, 'baseline') && ~isnan(blockOnset)
                    events(end+1).onset = (blockOnset - pulseZero)/1000;
                    events(end).duration = NaN; % 
                    events(end).trial_type = 'rest block';
                    events(end).stim_type = 'n/a';
                    events(end).button_press = 'n/a';

                    %Button press in rest block
                    for k = 1:length(keyTimes)
                        pressedKey = keyNames{k};
                        %Transform h/j/m/n -> g/f/v/b
                        switch pressedKey
                            case 'h', pressedKey = 'g';
                            case 'j', pressedKey = 'f';
                            case 'm', pressedKey = 'v';
                            case 'n', pressedKey = 'b';
                        end
                        events(end+1).onset = (keyTimes(k) - pulseZero)/1000;
                        events(end).duration = 0;
                        events(end).trial_type = 'button pressed';
                        events(end).stim_type = 'n/a';
                        events(end).button_press = pressedKey;
                    end
                end

                %Active block
                if strcmpi(blockName, 'active') && ~isnan(blockOnset)
                    stimDuration = 500;  
                    stimDelay = 16; %vertical blank      
                    nStim = 40; 
                    stimCodes = repmat(1:4, 1, 10); 

                    %Start of active block
                    events(end+1).onset = (blockOnset - pulseZero)/1000;
                    events(end).duration = NaN; 
                    events(end).trial_type = 'active block';
                    events(end).stim_type = 'n/a';
                    events(end).button_press = 'n/a';

                    stimOnset = blockOnset;
                    for s = 1:nStim
                        if s > 1
                            stimOnset = stimOnset + stimDuration + stimDelay;
                        end
                        stimEnd = stimOnset + stimDuration;

                        %Stimulus
                        events(end+1).onset = (stimOnset - pulseZero)/1000;
                        events(end).duration = stimDuration/1000;
                        events(end).trial_type = 'stimulus display';
                        events(end).stim_type = num2str(stimCodes(s));
                        events(end).button_press = 'n/a';

                        %Button press
                        respIdx = find(keyTimes >= stimOnset & keyTimes < stimEnd, 1);
                        if ~isempty(respIdx)
                            pressedKey = keyNames{respIdx};
                            switch pressedKey
                                case 'h', pressedKey = 'g';
                                case 'j', pressedKey = 'f';
                                case 'm', pressedKey = 'v';
                                case 'n', pressedKey = 'b';
                            end
                            events(end+1).onset = (keyTimes(respIdx) - pulseZero)/1000;
                            events(end).duration = 0;
                            events(end).trial_type = 'button pressed';
                            events(end).stim_type = 'n/a';
                            events(end).button_press = pressedKey;
                        end
                    end
                end
                i = j;
            else
                i = i + 1;
            end
        end

        %Calculating duration for block
        [~, idx] = sort([events.onset]);
        events = events(idx);

        blockIdx = find(strcmp({events.trial_type}, 'rest block') | strcmp({events.trial_type}, 'active block'));
        for b = 1:length(blockIdx)
            curIdx = blockIdx(b);
            if b < length(blockIdx)
                events(curIdx).duration = events(blockIdx(b+1)).onset - events(curIdx).onset;
            else
                events(curIdx).duration = (pulseTimes(end) - pulseZero)/1000 - events(curIdx).onset;
            end
        end

        %Save tsv file
        outPath = fullfile(outputDir, sprintf('sub-%s_task-VISMOT_acq-%s_events.tsv', subjID, trType));
        fid = fopen(outPath, 'w');
        fprintf(fid, 'onset\tduration\ttrial_type\tstim_type\tbutton_press\n');
        for k = 1:length(events)
            stimVal = events(k).stim_type;
            if isempty(stimVal), stimVal = 'n/a'; end
            buttonVal = events(k).button_press;
            if isempty(buttonVal), buttonVal = 'n/a'; end
            fprintf(fid, '%.3f\t%.3f\t%s\t%s\t%s\n', ...
                events(k).onset, events(k).duration, events(k).trial_type, ...
                stimVal, buttonVal);
        end
        fclose(fid);
    end
end
