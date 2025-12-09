% Author: B.J. Keulen
% Date: 05-11-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function for the extraction of BrainSense Setup data from json files.
% The data is saved into a folder within the location of the json file.
% 
% INPUT
%   info            =   struct containing information about json file
%   js              =   struct with data from json file
%   savepath        =   string with path into which to save dataSetup
%   savename_json   =   string of name under which to save dataSetup
%
% OUTPUT
%   dataSetup       =   struct containing info and data saved in
%                       MostRecentInSessonSignalCheck and SenseChannelTests

function dataSetup = getSetup(info, js, savepath, savename_json)

    % Initialise export structure
    dataSetup = [];
    dataSetup.DataType = 'MostRecentInSessionSignalCheck/SenseChannelTests';
    dataSetup.Info = info;
    dataSetup.Data = struct;
    sides = {'Left', 'Right'};
    fields = {'Run','DateTime','Channel','Time','LFP','Frequency','PSD','Artifact'};
    for s = 1:2
        for m = 1:2
            dataSetup.Data.(sides{s}) = struct;
        end
    end

    % Define channels
    chans = {'ZERO_THREE';'ONE_THREE';'ZERO_TWO'};
    
    % Loop over recordings
    for i = 1:size(js.SenseChannelTests,1)

        % Initialise structure for each recording
        setup = struct;
        for f = 1:length(fields)
            setup.(fields{f}) = [];
        end

        % Extract Setup data
        try 
            TD = js.SenseChannelTests{i};
        catch
            TD = js.SenseChannelTests(i);
        end

        % Get sampling frequency, if not filled in yet
        if ~isfield(dataSetup, 'SamplingFrequency')
            dataSetup.SamplingFrequency = TD.SampleRateInHz;
        end

        % Retrieve characteristics of recording
        parts = split(TD.Channel,'_');
        chan = strjoin(parts(1:end-1),'_');
        side = [parts{end}(1) lower(parts{end}(2:end))];

        % Get LFP data and construct time array
        setup.DateTime = char(datetime(strrep(TD.FirstPacketDateTime(1:end-5),'T',' ')) + hours(str2double(info.ProgrammerUtcOffset(1:3))) ...
                              + minutes(str2double(info.ProgrammerUtcOffset([1, 5:6]))));
        setup.UtcOffset = info.ProgrammerUtcOffset;
        setup.Channel = chan;
        setup.LFP = TD.TimeDomainData;
        setup.Time = (0 : 1/TD.SampleRateInHz : (length(setup.LFP)-1)/TD.SampleRateInHz)';

        % Retrieve corresponding frequency domain data if present
        try
            FD = js.MostRecentInSessionSignalCheck(contains({js.MostRecentInSessionSignalCheck.Channel},TD.Channel));
        catch
            if ~isempty(js.MostRecentInSessionSignalCheck)
                for m = 1:length(js.MostRecentInSessionSignalCheck)
                    if contains({js.MostRecentInSessionSignalCheck{m}.Channel},TD.Channel)
                        FD = js.MostRecentInSessionSignalCheck{m};
                    end
                end
            else
                FD = [];
            end
        end

        % If present, fill data into struct
        if ~isempty(FD)
            setup.Frequency = FD.SignalFrequencies;
            setup.Magnitude = FD.SignalPsdValues;
            if strcmp(FD.ArtifactStatus,'ArtifactStatusDef.ARTIFACT_NOT_PRESENT')
                setup.Artifact = 0;
            else
                setup.Artifact = 1;
            end
        end

        % Add data to main structure
        if isempty(fieldnames(dataSetup.Data.(side)))
            dataSetup.Data.(side) = setup;
        else
            dataSetup.Data.(side)(end+1) = setup;
        end
    end

    % Check for runs and empty fields
    nRuns = [];
    for s = 1:2

        % Check for number of runs and assign value to each recording
        if ~isempty(fieldnames(dataSetup.Data.(sides{s})))
            for c = 1:length(chans)
                idx = find(contains({dataSetup.Data.(sides{s}).Channel}, chans{c}));
                nRuns = [nRuns length(idx)];
                for r = 1:length(idx)
                    dataSetup.Data.(sides{s})(idx(r)).Run = r;
                end
            end

        % Remove field for hemisphere if empty
        else
            dataSetup.Data.(sides{s}) = rmfield(dataSetup.Data, sides{s});
        end
    end

    % Save number of runs
    dataSetup.nRuns = max(nRuns);

    % Save data
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    savename = [savename_json '_Setup'];
    save([savepath filesep savename '.mat'], 'dataSetup')

end