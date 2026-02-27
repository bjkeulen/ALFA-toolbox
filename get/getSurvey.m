% Author: B.J. Keulen
% Date: 26-07-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for the extraction of BrainSense Electrode Survey data from json
% files, recorded after the software update released in January 2025, which
% introduced the Electrode Identifier. The data is saved into a folder 
% within the location of the json file.
% 
% INPUT
%   info            =   struct containing information about json file
%   js              =   struct with data from json file
%   idx             =   index of recording, used to separate Survey from
%                       Identifier
%   savepath        =   string with path into which to save dataSurvey
%   savename_json   =   string of name under which to save dataSurvey
%
% OUTPUT
%   dataSurvey      =   struct containing info and Electrode Survey data

function dataSurvey = getSurvey(info, js, idx, savepath, savenameJSON)

    % Initialise export structure
    dataSurvey = [];
    if isfield(js, 'BrainSenseSurveysTimeDomain')
        dataSurvey.DataType = 'BrainSenseSurveys/BrainSenseSurveysTimeDomain';
        dataSurvey.Mode = js.BrainSenseSurveysTimeDomain{idx}.SurveyMode;
        allTime = js.BrainSenseSurveysTimeDomain{idx}.ElectrodeSurvey;
        allFrequency = js.BrainSenseSurveys{idx}.ElectrodeSurvey;
    else
        dataSurvey.DataType = 'LFPMontage/LfpMontageTimeDomain';
        dataSurvey.Mode = 'ElectrodeSurvey';
        allTime = js.LfpMontageTimeDomain;
        allFrequency = js.LFPMontage;
    end
    dataSurvey.Info = info;
    sides = {'Left', 'Right'};
    modes = {'Ring', 'Segment'};
    fields = {'Run','DateTime', 'Channel','Time','LFP','Frequency','PSD','Artifact'};
    for s = 1:2
        for m = 1:2
            dataSurvey.Data.(sides{s}).(modes{m}) = struct;
        end
    end

    % Define channels
    chans = cell(1,2);
    chans{1} = {'ZERO_AND_THREE';'ONE_AND_THREE';'ZERO_AND_TWO';'ONE_AND_TWO';'ZERO_AND_ONE';'TWO_AND_THREE'};
    chans{2} = {'ONE_A_AND_ONE_B';'ONE_B_AND_ONE_C';'ONE_A_AND_ONE_C';'TWO_A_AND_TWO_B';'TWO_B_AND_TWO_C'; ...
                'TWO_A_AND_TWO_C';'ONE_A_AND_TWO_A';'ONE_B_AND_TWO_B';'ONE_C_AND_TWO_C'};

    % Loop over recordings
    for i = 1:size(allTime,1)

        % Initialise structure for each recording
        for f = 1:length(fields)
            survey.(fields{f}) = [];
        end

        % Get time domain data
        try 
            TD = allTime{i};
        catch
            TD = allTime(i);
        end

        % Check if data is present
        try
            if isempty(TD.TimeDomainDatainMicroVolts)
                continue
            end
        catch
            if isempty(TD.TimeDomainData)
                continue
            end
        end

        % Get sampling frequency, if not filled in yet
        if ~isfield(dataSurvey,'SamplingFrequency')
            dataSurvey.SamplingFrequency = TD.SampleRateInHz;
        end

        % Retrieve characteristics of recording
        if isfield(TD, 'Hemisphere')
            side = [TD.Hemisphere(1) lower(TD.Hemisphere(2:end))];
        else
            if contains(TD.Channel, 'LEFT')
                side = 'Left';
            else
                side = 'Right';
            end
        end
        if contains(TD.Channel,{'_A_','_B_','_C_'})
            mode = 'Segment';
        else
            mode = 'Ring';
        end
        parts = split(TD.Channel,'_');
        chan = strjoin(parts(~contains(parts,{'RING','SEGMENT','LEFT','RIGHT'})),"_");

        % Get LFP data and construct time array
        survey.DateTime = char(datetime(strrep(TD.FirstPacketDateTime(1:end-5),'T',' ')) + hours(str2double(info.ProgrammerUtcOffset(1:3))) ...
                               + minutes(str2double(info.ProgrammerUtcOffset([1, 5:6]))));
        survey.UtcOffset = info.ProgrammerUtcOffset;
        survey.Channel = chan;
        try
            survey.LFP = TD.TimeDomainDatainMicroVolts;
        catch
            survey.LFP = TD.TimeDomainData;
        end
        survey.Time = (0 : 1/TD.SampleRateInHz : (length(survey.LFP)-1)/TD.SampleRateInHz)';

        % Retrieve corresponding frequency domain data if present
        try
            FD = allFrequency(all([contains({allFrequency.Hemisphere},[upper(side(1)),side(2:end)]); contains({allFrequency.SensingElectrodes},chan)]));
        catch
            for m = 1:length(allFrequency)
                if all([contains(allFrequency{m}.Hemisphere,[upper(side(1)),side(2:end)]); contains(allFrequency{m}.SensingElectrodes,chan)])
                    FD = allFrequency{m};
                end
            end
        end

        % If present, fill data into struct
        if ~isempty(FD)
            try
                survey.Frequency = FD.LFPFrequencyinHertz;
                survey.Magnitude = FD.LFPMagnitudeinMicroVoltPeak;
            catch
                survey.Frequency = FD.LFPFrequency;
                survey.Magnitude = FD.LFPMagnitude;
            end
            if strcmp(FD.ArtifactStatus,'ArtifactStatusDef.ARTIFACT_NOT_PRESENT')
                survey.Artifact = 0;
            else
                survey.Artifact = 1;
            end
        end

        % Add data to main structure
        if isempty(fieldnames(dataSurvey.Data.(side).(mode)))
            dataSurvey.Data.(side).(mode) = survey;
        else
            dataSurvey.Data.(side).(mode)(end+1) = survey;
        end
    end

    % Check for runs and empty fields
    nRuns = [];
    for s = 1:2
        for m = 1:2

            % Check for number of runs and assign value to each recording
            if ~isempty(fieldnames(dataSurvey.Data.(sides{s}).(modes{m})))
                for c = 1:length(chans{m})
                    iRec = find(contains({dataSurvey.Data.(sides{s}).(modes{m}).Channel}, chans{m}{c}));
                    nRuns = [nRuns length(iRec)];
                    for r = 1:length(iRec)
                        dataSurvey.Data.(sides{s}).(modes{m})(iRec(r)).Run = r;
                    end
                end

            % Remove field for mode if empty
            else
                dataSurvey.Data.(sides{s}) = rmfield(dataSurvey.Data.(sides{s}), modes{m});
            end
        end

        % Remove field for hemisphere if empty
        if isempty(fieldnames(dataSurvey.Data.(sides{s})))
            dataSurvey.Data = rmfield(dataSurvey.Data, sides{s});
        end
    end

    % Save number of runs
    dataSurvey.nRuns = max(nRuns);

    % Save data
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    save([savepath filesep savenameJSON '_Survey'], 'dataSurvey')

end