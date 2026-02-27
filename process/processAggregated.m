% Author: B.J. Keulen
% Date: 14-01-2026
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for further processing of aggregated Events and Timeline data.
% Duplicates are removed, missing values and parameters are added (Timeline
% only) and datetime values are time zone corrected.
% 
% INPUT
%   dataStruct      =   struct with aggregated dataTimeline or dataEvents
%                       data
%   dataTL          =   dataTimeline. Only required if dataEvents is given
%                       as input for dataStruct
%   fields          =   cell with names of fields in dataTimeline and 
%                       dataEvents which contain extracted LFP data. Only 
%                       required if datatimeline is given as input for 
%                       dataStruct
%   dataset         =   str defining type of dataset. Either single file 
%                       (0), a single folder (1) or a set of folders (2)
%   fileData        =   struct containing the path and name of the folder
%                       containing the selected JSON files, the file names
%                       and number of folders and files, all if applicable
%   tZone           =   str defining Timezone to use for correcting for UTC
%                       offset and daylight saving time
%   savepath        =   str defining path to which output structures should
%                       be saved
%   filename        =   str with name of file. Only required when using
%                       single file processing (dataset = 0)
%   f               =   int indicating index of folder. Only required when
%                       using folderset processing (dataset = 2)
%
% OUTPUT
%   dataStruct      =   dataTimeline or dataEvents after further processing

function dataStruct = processAggregated(dataStruct, dataTL, fields, dataset, fileData, tZone, savepath, filename, f)

    % Further processing Timeline data
    if strcmp(dataStruct.DataType, 'LFPTrendLogs')
        dataStruct.Data = sortrows(struct2table(dataStruct.Data),'DateTime');     % Convert to table and sort on datetime       
        dataStruct = checkDuplicates(dataStruct, fields);                         % Remove duplicates
        dataStruct = addNaN(dataStruct);                                          % Insert NaN for missing values
        dataStruct = addParameters(dataStruct, tZone);                            % Get sensing and stimulation parameters and assign to data

    % Further processing Events data
    elseif strcmp(dataStruct.DataType, 'LfpFrequencySnapshotEvents') 
        dataStruct.Data = dataStruct.Data(2:end);                                 % Remove first empty row
        if length(dataStruct.Data) > 1                                            % Convert to table and sort on datetime
            dataStruct.Data = sortrows(struct2table(dataStruct.Data),'DateTime');     
        end
        dataStruct = checkDuplicates(dataStruct);                                 % Remove duplicates
        
    else
        error('Incorrect input structure. Only dataTimeline and dataEvents are allowed.')
    end

    % Change timezone to correct for UTC offset
    dataStruct.Data.DateTime.TimeZone = 'UTC';
    dataStruct.Data.DateTime.TimeZone = tZone;
    dataStruct.TimeZone = tZone;

    % Save data
    if strcmp(dataStruct.DataType, 'LFPTrendLogs')

        % Duplicate struct with different name for saving
        dataTimeline = dataStruct;

        % Save Timeline data
        if dataset == 0
            save([savepath filesep filename(1:end-5) '_Timeline.mat'], 'dataTimeline')
        elseif dataset == 1
            save([savepath filesep fileData.rootName '_Timeline.mat'], 'dataTimeline')
        elseif dataset == 2
            if ~exist([savepath filesep 'Timeline'], 'dir')
               mkdir([savepath filesep 'Timeline'])
            end
            save([savepath filesep 'Timeline' filesep fileData.folders(f).name '_Timeline.mat'], 'dataTimeline')
        end

    elseif strcmp(dataStruct.DataType, 'LfpFrequencySnapshotEvents') 

        % Add stimulation amplitude, pulsewidth and frequency to Events
        dataStruct = addStimInfo(dataStruct, dataTL);

        % Duplicate struct with different name for saving
        dataEvents = dataStruct;

        % Save Events data
        if dataset == 0
            save([savepath filesep filename(1:end-5) '_Events.mat'], 'dataEvents')
        elseif dataset == 1
            save([savepath filesep fileData.rootName '_Events.mat'], 'dataEvents')
        elseif dataset == 2
            if ~exist([savepath filesep 'Events'], 'dir')
               mkdir([savepath filesep 'Events'])
            end
            save([savepath filesep 'Events' filesep fileData.folders(f).name '_Events.mat'], 'dataEvents')
        end
    end
end
