% Author: B.J. Keulen
% Date: 13-09-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function for assigning stimulation parameters as stored in dataTimeline
% to dataEvents. Data are assigned if the nearest known value was less than
% 15 minutes from the time of event. If no Timeline data were recorded
% near the time of events, the stimulation parameters remain unknown.
% 
% INPUT
%   dataEvents      =   struct with Events data
%   dataTimeline    =   struct with Timeline data
%
% OUTPUT
%   dataEvents      =   struct with Events data

function dataEvents = addStimInfo(dataEvents, dataTimeline)

    % Check if dataTimeline were recorded
    if ~isempty(dataTimeline.Data.LfpPower)

        % Loop over events
        for e = 1:height(dataEvents.Data)
    
            % Find index of Timeline time array closest to Event time
            [delta, idx] = min(abs(dataTimeline.Data.DateTime - dataEvents.Data{e,'DateTime'}));
    
            % Add stimulation amplitude, pulsewidth and frequency if difference is less than 15 min
            if delta < duration('00:15:00')
                dataEvents.Data{e,'StimulationAmplitude'} = dataTimeline.Data.StimulationAmplitude(idx,:);
                dataEvents.Data{e,'PulseWidth'} = dataTimeline.Data.PulseWidth(idx,:);
                dataEvents.Data{e,'StimulationFrequency'} = dataTimeline.Data.StimulationFrequency(idx,:);
            end
        end
    end
end