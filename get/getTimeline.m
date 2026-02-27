% Author: B.J. Keulen
% Date: 29-03-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for the extraction of BrainSense Timeline data from a single
% json file. This data will be written into the larger dataTimeline struct
% in the main script.
% 
% Part of the structure of this code is adapted from the perceive.m file
% of the perceive toolbox by Wolf-Julian Neumann, Tomas Sieger and Gerd 
% Tinkhauser: https://github.com/neuromodulation/perceive.
% 
% INPUT
%   data            =   js.DiagnosticData.LFPTrendLogs, with
%                       js the raw data from the json file
%   fieldsTL        =   cell with names of fields in dataTimeline
%
% OUTPUT
%   trendlogs       =   struct with all timeline data within json file

function timeline = getTimeline(data, fieldsTL)

    % Initialise export structure
    timeline = struct;
    for n = 1:length(fieldsTL)
        timeline.(fieldsTL{n}) = [];
    end

    % Initiate arrays
    lfpL=[]; stimL=[]; dtL=datetime([],[],[]);
    lfpR=[]; stimR=[]; dtR=datetime([],[],[]);

    % Check if data from left hemisphere is present
    if isfield(data.LFPTrendLogs,'HemisphereLocationDef_Left')

        % Get data and runs
        left = data.LFPTrendLogs.HemisphereLocationDef_Left;
        days = fieldnames(left);

        for d=1:length(days)

            % Get data for each run
            clfp = [left.(days{d}).LFP];
            cstim = [left.(days{d}).AmplitudeInMilliAmps];
            cdt = datetime({left.(days{d}).DateTime},'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
            [cdt,i] = sort(cdt);

            % Append data of each run to one array
            lfpL = [lfpL,clfp(i)];
            stimL = [stimL,cstim(i)];
            dtL = [dtL,cdt];

        end
    end

    % Check if data from right hemisphere is present
    if isfield(data.LFPTrendLogs,'HemisphereLocationDef_Right')

        % Get data and runs
        right = data.LFPTrendLogs.HemisphereLocationDef_Right;
        days = fieldnames(right);

        for d=1:length(days)

            % Get data for each run
            clfp = [right.(days{d}).LFP];
            cstim = [right.(days{d}).AmplitudeInMilliAmps];
            cdt = datetime({right.(days{d}).DateTime},'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''');
            [cdt,i] = sort(cdt);
            
            % Append data of each run to one array
            lfpR = [lfpR,clfp(i)];
            stimR = [stimR,cstim(i)];
            dtR = [dtR,cdt];

        end
    end

    % Check presence of left/right datetime data and get all unique values
    if isempty(dtL)
        time = sort(dtR)';
    elseif isempty(dtR)
        time = sort(dtL)';
    else
        time=sort([dtL,setdiff(dtR,dtL)])';
    end

    % Assign LFP values of each side to one array corresponding to datetime
    lfp = nan(length(time),2);
    stim = nan(length(time),2);

    for t = 1:length(time)
        if ismember(time(t),dtL)
            i = find(dtL == time(t));
            lfp(t,1) = lfpL(i);
            stim(t,1) = stimL(i);
        else
            lfp(t,1) = nan;
            stim(t,1) = nan;
        end

        if ismember(time(t),dtR)
            i = find(dtR == time(t));
            lfp(t,2) = lfpR(i);
            stim(t,2) = stimR(i);
        else
            lfp(t,2) = nan;
            stim(t,2) = nan;
        end
    end

    % Add data to structure
    timeline.DateTime = time;
    timeline.ActiveGroup = strings(length(time),1);
    timeline.SensingChannel = cell(length(time),2);
    timeline.SensingFrequency = nan(length(time),2);
    timeline.LfpPower = lfp;
    timeline.LowerLfpThreshold = nan(length(time),2);
    timeline.UpperLfpThreshold = nan(length(time),2);
    timeline.StimulationAmplitude = stim;
    timeline.LowerStimulationLimit = nan(length(time),2);
    timeline.UpperStimulationLimit = nan(length(time),2);
    timeline.PulseWidth = nan(length(time),2);
    timeline.StimulationFrequency = nan(length(time),2);
    
end