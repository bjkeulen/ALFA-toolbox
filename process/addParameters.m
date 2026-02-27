% Author: B.J. Keulen
% Date: 17-09-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Add group, stimulation frequency, pulsewidth, sensing channel, sensing 
% frequency, stimulation limits and LFP thresholds to dataTimeline.
%
% INPUT
%   dataTimeline    =   struct with Timeline data
%
% OUTPUT
%   dataTimeline    =   struct with Timeline data

function dataTimeline = addParameters(dataTimeline, tZone)

    % Retrieve sessions and group changes from info and EventLogs
    [sessions, changes] = getHistory(dataTimeline);

    % Assign groups to data based on group changes
    if isempty(changes)
        dataTimeline.Data{dataTimeline.Data{:,"DateTime"} < sessions{1,"SessionStart"},"ActiveGroup"} = sessions{1,"InitialActive"};
        for c = 1:height(sessions)-1
            dataTimeline.Data{isbetween(dataTimeline.Data{:,"DateTime"}, sessions{c,"SessionEnd"}, ...
                              sessions{c+1,"SessionStart"}, "openright"),"ActiveGroup"} = sessions{c,"FinalActive"};
        end
    else
        dataTimeline.Data{dataTimeline.Data{:,"DateTime"} < changes{1,"DateTime"},"ActiveGroup"} = changes{1,"OldGroup"};
        dataTimeline.Data{dataTimeline.Data{:,"DateTime"} >= changes{end,"DateTime"},"ActiveGroup"} = changes{end,"NewGroup"};
        for c = 1:height(changes)-1
            dataTimeline.Data{isbetween(dataTimeline.Data{:,"DateTime"}, changes{c,"DateTime"}, ...
                              changes{c+1,"DateTime"}, "openright"),"ActiveGroup"} = changes{c,"NewGroup"};
        end
    end
    
    % Assign group info to data using active groups and session datetimes
    grID = str2double(replace(dataTimeline.Data{:,'ActiveGroup'},{'GROUP_','A','B','C','D'},{'','1','2','3','4'}));
    for h = 1:height(sessions)
        if h == 1
            idx = dataTimeline.Data{:,"DateTime"} < sessions{1,"SessionStart"};
        else
            idx = isbetween(dataTimeline.Data{:,"DateTime"}, sessions{h-1,"SessionEnd"}, sessions{h,"SessionStart"}, "open");
        end
        if ~isempty(sessions.InitialSensingChannel{h}) && ~isempty(find(idx,1))
            dataTimeline.Data{idx,'SensingFrequency'} = sessions.InitialSensingFrequency{h}(grID(idx),:);
            dataTimeline.Data{idx,'SensingChannel'} = sessions.InitialSensingChannel{h}(grID(idx),:);
            dataTimeline.Data{idx,'LowerLfpThreshold'} = sessions.InitialLowerLfpThreshold{h}(grID(idx),:);
            dataTimeline.Data{idx,'UpperLfpThreshold'} = sessions.InitialUpperLfpThreshold{h}(grID(idx),:);
            dataTimeline.Data{idx,'LowerStimulationLimit'} = sessions.InitialLowerAmplitudeLimit{h}(grID(idx),:);
            dataTimeline.Data{idx,'UpperStimulationLimit'} = sessions.InitialUpperAmplitudeLimit{h}(grID(idx),:);
            dataTimeline.Data{idx,'StimulationFrequency'} = sessions.InitialStimulationFrequency{h}(grID(idx),:);
            dataTimeline.Data{idx,'PulseWidth'} = sessions.InitialPulseWidth{h}(grID(idx),:);
        end
    end

    % Change timezone
    sessions.SessionStart.TimeZone = 'UTC';
    sessions.SessionStart.TimeZone = tZone;
    sessions.SessionEnd.TimeZone = 'UTC';
    sessions.SessionEnd.TimeZone = tZone;
    changes.DateTime.TimeZone = 'UTC';
    changes.DateTime.TimeZone = tZone;

    % Add history to dataTimeline
    dataTimeline.History.Sessions = sessions;
    dataTimeline.History.Changes = changes;

end