%% EVENTS FILE FOR TASK VOB

function EVENT_VOB
    %Define path for input and output dir
    inputDir = 'Y:\projekty\ELICIT\derivatives\spm_processed_dataset\EVENTS\VOB';       
    outputDir = 'Y:\projekty\ELICIT\derivatives\spm_processed_dataset\EVENTS\VOB'; 

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %Find all logs which need to reformat
    logFiles = dir(fullfile(inputDir, 'ELICIT_VOB*.txt')); 

    %For loop across all log files
    for f = 1:length(logFiles)
        logPath = fullfile(logFiles(f).folder, logFiles(f).name);

        %Find ID of subject
        subjID = regexp(logFiles(f).name, 'ELICIT_VOBdistr_2_v[23]-(\d+)-\d+', 'tokens');   
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

        if length(pulseTimes) < 2
            warning('Two PulseTime are required in file %s.', logFiles(f).name);
            continue;
        end

        pulseZero = pulseTimes(1);
        TR = pulseTimes(2) - pulseTimes(1);

        if abs(TR - 1800) <= 5
            trType = 'tr1800';
        elseif abs(TR - 800) <= 5
            trType = 'tr800';
        else
            warning('Unknown TR');
            trType = sprintf('tr%d', TR);
        end
    
        %Name of output tsv file
        outPath = fullfile(outputDir, sprintf('sub-%s_task-VOB_acq-%s_events.tsv', subjID, trType));

        %Define the structure
        events = struct('onset', {}, 'duration', {}, 'trial_type', {}, ...
                    'stimulus', {}, 'response_time', {}, 'correct', {});

        %Work with logs
        i = 1;
        while i <= length(lines)
            if contains(lines{i}, '*** LogFrame Start ***')
                stimCode = '';
                stimText = '';
                stimOnset = NaN;
                stimDuration = NaN;
                responseTime = NaN;
                keyPressed = 'n/a';


                j = i + 1;
                while j <= length(lines) && ~contains(lines{j}, '*** LogFrame End ***')
                    line = strtrim(lines{j});
                    if startsWith(line, 'Stimulus: ')
                        stimCode = strtrim(extractAfter(line, 'Stimulus: '));
                    elseif startsWith(line, 'Stimulus_text: ')
                        stimText = strtrim(extractAfter(line, 'Stimulus_text: '));
                    elseif startsWith(line, 'Stimulus_Onset: ')
                        stimOnset = str2double(extractAfter(line, 'Stimulus_Onset: '));
                    elseif startsWith(line, 'Stimulus_Duration: ')
                        stimDuration = str2double(extractAfter(line, 'Stimulus_Duration: '));
                    elseif startsWith(line, 'Key_')
                        keyPressed = strtrim(extractAfter(line, ': '));
                    elseif startsWith(line, 'Time_')
                        responseTime = str2double(extractAfter(line, ': ')) - stimOnset;
                    end
                    j = j + 1;
                end

                %Define the type of stimulus
                if ~isempty(stimCode)
                    switch stimCode
                        case 'T'
                            trialType = 'target';
                        case 'F'
                            trialType = 'frequent';
                        case 'D'
                            trialType = 'distractor';
                        otherwise
                            trialType = 'unknown';
                    end
                

                    %Define accuracy
                    if strcmp(trialType, 'target')
                        if strcmp(keyPressed, 'n/a')
                            acc = 0;
                        else
                            acc = 1;
                        end
                    elseif ismember(trialType, {'frequent', 'distractor'})
                        if strcmp(keyPressed, 'n/a')
                            acc = 1;
                        else
                            acc = 0;
                        end
                    else
                        acc = NaN;
                    end

                    events(end+1).onset = (stimOnset - pulseZero) / 1000;
                    events(end).duration = stimDuration / 1000;
                    events(end).trial_type = trialType;
                    events(end).stimulus = stimText;
                    if isnan(responseTime)
                        events(end).response_time = 'n/a';
                    else
                        events(end).response_time = sprintf('%.3f', responseTime / 1000);
                    end
                    events(end).accuracy = acc;
                end
                i = j;
            else
                i = i + 1;
            end
        end

        %Save tsv file
        fid = fopen(outPath, 'w');
        fprintf(fid, 'onset\tduration\ttrial_type\tstimulus\tresponse_time\tcorrect\n');
        for k = 1:length(events)
            fprintf(fid, '%.3f\t%.3f\t%s\t%s\t%s\t%d\n', ...
                events(k).onset, events(k).duration, events(k).trial_type, events(k).stimulus, ...
                events(k).response_time, events(k).accuracy);
        end
        fclose(fid);
    end
end
