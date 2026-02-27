% Author: B.J. Keulen
% Date: 29-03-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for the extraction of BrainSense Events data of a single json
% file. This data will be written into the larger dataEvents struct
% in the main script.
%
% Part of the structure of this code is adapted from the perceive.m file
% of the perceive toolbox by Wolf-Julian Neumann, Tomas Sieger and Gerd 
% Tinkhauser: https://github.com/neuromodulation/perceive.
% 
% INPUT
%   data            =   js.DiagnosticData.LfpFrequencySnapshotEvents, with
%                       js the raw data from the json file
%   fieldsSE        =   cell with names of fields in dataEvents which
%                       contain extracted event data
%
% OUTPUT
%   events          =   1xn struct with n the number of events in the json
%                       file 

function events = getEvents(data, fieldsSE)

    % Initialise export structure
    events = struct;
    events.EventName = {};
    for n = 1:length(fieldsSE)
        events.(fieldsSE{n}) = [];
    end
    e = 1;

    % Loop over events
    for c = 1:length(data)
        try 
            eventdata = data{c};
        catch
            eventdata = data(c);
        end

        % Get data from event
        id = eventdata.EventID;
        eventname = {eventdata.EventName};
        time = datetime(eventdata.DateTime(1:end-1),'InputFormat','yyyy-MM-dd''T''HH:mm:ss');

        % Create empty array for channel and PSD values
        channel = {{}, {}};
        psd = {{},{}};

        if isfield(eventdata,'LfpFrequencySnapshotEvents')

            % Get left lead contacts used for measurement
            if isfield(eventdata.LfpFrequencySnapshotEvents,'HemisphereLocationDef_Left')
                chanL = split(eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.SenseID,'.');
                channel(:,1) = chanL(end);
            end
            
            % Get right lead contacts used for measurement
            if isfield(eventdata.LfpFrequencySnapshotEvents,'HemisphereLocationDef_Right')
                chanR = split(eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Right.SenseID,'.');
                channel(:,2) = chanR(end);
            end

            % Get frequencies and bin data
            if isfield(eventdata.LfpFrequencySnapshotEvents,'HemisphereLocationDef_Left') && isfield(eventdata.LfpFrequencySnapshotEvents,'HemisphereLocationDef_Right')
                freq = eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.Frequency;
                psd(1) = {eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.FFTBinData};
                psd(2) = {eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Right.FFTBinData};

            elseif isfield(eventdata.LfpFrequencySnapshotEvents,'HemisphereLocationDef_Left')
                freq = eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.Frequency;
                psd(1) = {eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Left.FFTBinData};

            elseif isfield(eventdata.LfpFrequencySnapshotEvents,'HemisphereLocationDef_Right')
                freq = eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Right.Frequency;
                psd(2) = {eventdata.LfpFrequencySnapshotEvents.HemisphereLocationDef_Right.FFTBinData};

            else
               freq = [];
            end

        else
            freq = {};
        end

        % Add data to structure
        events(e).DateTime = time;
        events(e).EventName = eventname;
        events(e).EventID = id;
        events(e).SensingChannel = channel;
        events(e).Frequency = freq;
        events(e).Magnitude = psd;
        events(e).StimulationAmplitude = [nan, nan];
        events(e).StimulationFrequency = [nan, nan];
        events(e).PulseWidth = [nan, nan];
        e = e + 1;

    end
end