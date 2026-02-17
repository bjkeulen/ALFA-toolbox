% Author: B.J. Keulen
% Date: 23-01-2026
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to save an Excel file containing a data log with info about 
% which data is stored within individual JSON files, and a log containing 
% info about the number of total, valid and missing datapoints of 
% aggregated Timeline and Events data.
% 
% INPUT
%   dataLog         =   table with JSON filenames and info about which data 
%                       is stored within individual JSON files
%   dataTimeline    =   struct containing aggregated and processed Timeline
%                       data
%   dataEvents      =   struct containing aggregated and processed Events
%                       data
%   savepath        =   string with path into which to save data logs
%   savename        =   string with name under which to save data logs
%   rootName        =   string with name of folderset
%
% OUTPUT
%   none

function saveLogs(dataLog, dataTimeline, dataEvents, savepath, savename, rootName)

    % Create table for Timeline and Events logs
    aggrLog = table('Size',[2,7], 'VariableTypes',{'string','int64','int64','int64','int64','datetime','datetime'}, ...
                          'VariableNames',{'Datatype','n_files','n_total','n_valid','n_missing','start','end'});
    aggrLog{:,'Datatype'} = {'Timeline';'Events'};

    % Collect Timeline log
    if ~isempty(dataTimeline.Data.LfpPower)
        aggrLog{1,'n_files'} = length(dataTimeline.Info);
        aggrLog{1,'n_total'} = height(dataTimeline.Data);
        aggrLog{1,'n_valid'} = height(dataTimeline.Data(~any(isnan(dataTimeline.Data.LfpPower),2),:));
        aggrLog{1,'n_missing'} = height(dataTimeline.Data(any(isnan(dataTimeline.Data.LfpPower),2),:));

        dtTL = dataTimeline.Data.DateTime;
        dtTL.TimeZone = '';
        aggrLog{1,'start'} = dtTL(1);
        aggrLog{1,'end'} = dtTL(end);
    end

    % Collect Events log
    if ~isempty([dataEvents.Data.EventName])
        aggrLog{2,'n_files'} = length(dataEvents.Info);
        aggrLog{2,'n_total'} = height(dataEvents.Data);
        for i = 1:height(dataEvents.Data)
            if isempty(dataEvents.Data.Frequency{i})
                aggrLog{2,'n_missing'} = aggrLog{2,'n_missing'} + 1;
            else
                aggrLog{2,'n_valid'} = aggrLog{2,'n_valid'} + 1;
            end
        end

        dtSE = dataEvents.Data.DateTime;
        dtSE.TimeZone = '';
        aggrLog{2,'start'} = dtSE(1);
        aggrLog{2,'end'} = dtSE(end);
    end

    % Set savepath for logs and create folder if folderset processing
    if isempty(rootName)
        savepathLogs = savepath;
    else
        savepathLogs = [savepath filesep 'Logs'];
        if ~exist(savepathLogs, 'dir')
           mkdir(savepathLogs)
        end
    end

    % Save logs
    writetable(dataLog, [savepathLogs filesep 'log_' savename '.xlsx'], 'Sheet', 'Single JSON')
    writetable(aggrLog, [savepathLogs filesep 'log_' savename '.xlsx'], 'Sheet', 'Aggregated data')

    % Display results in command window
    types = ["Setup" "Survey" "Identifier" "Streaming" "Timeline" "Events"];
    if isempty(rootName)
        fprintf('------------------\nProcessing completed of %s\n  Single JSON files:\n', savename)
    else
        fprintf('------------------\nProcessing completed of %s in folderset %s\n  Single JSON files:\n', savename, rootName)
    end
    for f = 1:height(dataLog)
        if strcmp(dataLog{f,'File_check'}, 'Pass')
            if sum(dataLog{f,3:end}) == 0
                fprintf('    %s: No LFP data\n', dataLog{f,'Filename'})
            else
                fprintf('    %s: %s\n', dataLog{f,'Filename'}, strjoin(types(logical(dataLog{f,3:end})),', '))
            end
        else
            fprintf('    %s: %s\n', dataLog{f,'Filename'}, dataLog{f,'File_check'})
        end
    end
    if ~isempty(dataTimeline.Data.LfpPower) || ~isempty([dataEvents.Data.EventName])
        fprintf('  Aggregated data:\n')
        if ~isempty(dataTimeline.Data.LfpPower)
            fprintf('    Timeline: %i data points from %s to %s\n', aggrLog{1,'n_total'}, aggrLog{1,"start"}, aggrLog{1,'end'})
        end
        if ~isempty([dataEvents.Data.EventName])
            fprintf('    Events: %i recordings from %s to %s\n', aggrLog{2,'n_total'}, aggrLog{2,"start"}, aggrLog{2,'end'})
        end
    else
        fprintf('  No aggregated data\n')
    end
    fprintf('Data saved in %s\n------------------\n', savepath)
end