% Author: B.J. Keulen
% Date: 14-02-2025
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for the extraction of BrainSense Electrode Identifier data from
% json files, recorded after the software update released in January 2025,
% which introduced the Electrode Identifier. The data is saved into a
% folder within the location of the json file.
% 
% INPUT
%   info            =   struct containing information about json file
%   js              =   struct with data from json file
%   idx             =   index of recording, used to separate Survey from
%                       Identifier
%   savepath        =   string with path into which to save dataIdentifier
%   savename_json   =   string of name under which to save dataIdentifier
%
% OUTPUT
%   dataIdentifier  =   struct containing info and Electrode Identifier 
%                       data

function dataIdentifier = getIdentifier(info, js, idx, savepath, savenameJSON)

    % Initialise export structure
    dataIdentifier = [];
    dataIdentifier.DataType = 'BrainSenseSurveys/BrainSenseSurveysTimeDomain';
    dataIdentifier.Mode = js.BrainSenseSurveysTimeDomain{idx}.SurveyMode;
    dataIdentifier.Info = info;
    sides = {'Left', 'Right'};
    modes = {'Ring', 'Segment'};
    fields = {'Run','DateTime', 'Channel','Time','LFP','Frequency','PSD','Artifact'};
    for s = 1:2
        for m = 1:2
            dataIdentifier.Data.(sides{s}).(modes{m}) = struct;
        end
    end

    % Define channels
    chans = cell(1,2);
    chans{1} = {'ZERO';'ONE';'TWO';'THREE'};
    chans{2} = {'ONE_A';'ONE_B';'ONE_C';'TWO_A';'TWO_B';'TWO_C'};

    % Loop over recordings
    for i = 1:size(js.BrainSenseSurveysTimeDomain{idx}.ElectrodeIdentifier,1)

        % Initialise structure for each recording
        for f = 1:length(fields)
            identifier.(fields{f}) = [];
        end

        % Get time domain data
        try 
            TD = js.BrainSenseSurveysTimeDomain{idx}.ElectrodeIdentifier{i};
        catch
            TD = js.BrainSenseSurveysTimeDomain{idx}.ElectrodeIdentifier(i);
        end

        % Check if data is present
        if isempty(TD.TimeDomainDatainMicroVolts)
            continue
        end

        % Get sampling frequency, if not filled in yet
        if ~isfield(dataIdentifier,'SamplingFrequency')
            dataIdentifier.SamplingFrequency = TD.SampleRateInHz;
        end

        % Retrieve characteristics of recording
        side = [TD.Hemisphere(1) lower(TD.Hemisphere(2:end))];
        if contains(TD.Channel,{'_A','_B','_C'})
            mode = 'Segment';
        else
            mode = 'Ring';
        end
        parts = split(TD.Channel,'_');
        chan = strjoin(parts(~contains(parts,'ELECTRODE')),"_");

        % Get LFP data and construct time array
        identifier.DateTime = char(datetime(strrep(TD.FirstPacketDateTime(1:end-5),'T',' ')) + hours(str2double(info.ProgrammerUtcOffset(1:3))) ...
                                   + minutes(str2double(info.ProgrammerUtcOffset([1, 5:6]))));
        identifier.UtcOffset = info.ProgrammerUtcOffset;
        identifier.Channel = chan;
        identifier.LFP = TD.TimeDomainDatainMicroVolts;
        identifier.Time = (0 : 1/TD.SampleRateInHz : (length(identifier.LFP)-1)/TD.SampleRateInHz)';

        % Retrieve corresponding frequency domain data if present
        try
            FD = js.BrainSenseSurveys{idx}.ElectrodeIdentifier(all([contains({js.BrainSenseSurveys{idx}.ElectrodeIdentifier.Hemisphere},[upper(side(1)),side(2:end)]); contains({js.BrainSenseSurveys{idx}.ElectrodeIdentifier.SensingElectrodes},chan)]));
        catch
            for m = 1:length(js.BrainSenseSurveys{idx}.ElectrodeIdentifier)
                if all([contains(js.BrainSenseSurveys{idx}.ElectrodeIdentifier{m}.Hemisphere,[upper(side(1)),side(2:end)]); contains(js.BrainSenseSurveys{idx}.ElectrodeIdentifier{m}.SensingElectrodes,chan)])
                    FD = js.BrainSenseSurveys{idx}.ElectrodeIdentifier{m};
                end
            end
        end

        % If present, fill data into struct
        if ~isempty(FD)
            identifier.Frequency = FD.LFPFrequencyinHertz;
            identifier.Magnitude = FD.LFPMagnitudeinMicroVoltPeak;
            if strcmp(FD.ArtifactStatus,'ARTIFACT_NOT_PRESENT')
                identifier.Artifact = 0;
            else
                identifier.Artifact = 1;
            end
        end

        % Add data to main structure
        if isempty(fieldnames(dataIdentifier.Data.(side).(mode)))
            dataIdentifier.Data.(side).(mode) = identifier;
        else
            dataIdentifier.Data.(side).(mode)(end+1) = identifier;
        end
    end

    % Check for runs and empty fields
    nRuns = [];
    for s = 1:2
        for m = 1:2

            % Check for number of runs and assign value to each recording
            if ~isempty(fieldnames(dataIdentifier.Data.(sides{s}).(modes{m})))
                for c = 1:length(chans{m})
                    idx = find(contains({dataIdentifier.Data.(sides{s}).(modes{m}).Channel}, chans{m}{c}));
                    nRuns = [nRuns length(idx)];
                    for r = 1:length(idx)
                        dataIdentifier.Data.(sides{s}).(modes{m})(idx(r)).Run = r;
                    end
                end

            % Remove field for mode if empty
            else
                dataIdentifier.Data.(sides{s}) = rmfield(dataIdentifier.Data.(sides{s}), modes{m});
            end
        end

        % Remove field for hemisphere if empty
        if isempty(fieldnames(dataIdentifier.Data.(sides{s})))
            dataIdentifier.Data = rmfield(dataIdentifier.Data, sides{s});
        end
    end

    % Save number of runs
    dataIdentifier.nRuns = max(nRuns);

    % Save data
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    save([savepath filesep savenameJSON '_Identifier'], 'dataIdentifier')

end