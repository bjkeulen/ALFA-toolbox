% Author: B.J. Keulen, M.J. Stam & J.T. Boonstra
% Date: 29-03-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for the extraction of BrainSenseTimeDomain Streaming data and
% reconstruction of the corresponding time vector. The LFP data is notch
% filtered, ECG-suppressed if applicable, and saved. The output struct
% contains for each recording within JSON the info, BrainSenseTimeDomain 
% Streaming time array and LFP data in raw, notch filtered, ECG suppressed 
% and notch filtered + ECG suppressed form.
%
% Part of the structure of this code is adapted from the perceive.m file
% of the perceive toolbox by Wolf-Julian Neumann, Tomas Sieger and Gerd 
% Tinkhauser: https://github.com/neuromodulation/perceive.
% 
% INPUT
%   info            =   struct containing information about json file
%   js              =   struct with data from json file
%   ecgFilter       =   int indicating type of ECG suppression
%   linenoise       =   1x2 double with frequency and bandwidth for 
%                       line-noise filter
%   savepath        =   string with path into which to save dataStreaming
%   savename_json   =   string of name under which to save dataStreaming
%   plotData        =   int indicating whether to plot data (1) or not (0)
%   showFig         =   int indicating whether to show figures (1) or not
%                       (0)
%
% OUTPUT
%   dataStreaming   =   1xn struct containing Streaming data with n rows
%                       corresponding to each recording within JSON

function dataStreaming = getStreaming(info, js, linenoise, ecgMethod, rTime, savepath, savenameJSON, plotData, showFig)

    % Extract required fields
    TDdata = js.BrainSenseTimeDomain;
    BSdata = js.BrainSenseLfp;

    % Check number of recordings in file
    FirstPacketDateTime = strrep(strrep({TDdata(:).FirstPacketDateTime},'T',' '),'Z','');

    % ----- below ---- is adjusted !! MJS

    % runs = unique(FirstPacketDateTime);

    ChannelRaw = {TDdata(:).Channel};

    % Obtain LEFT/RIGHT after last "_" in ChannelRaw:
    out  = regexp(ChannelRaw, '[^_]+$', 'match');
    side = cellfun(@(x) x{1}, out, 'UniformOutput', false);   % Nx1 cell: 'LEFT'/'RIGHT'

    % Group on timestamp
    [uTimes, ~, g] = unique(FirstPacketDateTime, 'stable');

    runCount = 0;

    for k = 1:numel(uTimes)
        idx = find(g == k);              % all rows with same timestamp (in order)
        s   = side(idx);                 % hemispheres at this timestamp

        hasL = any(strcmpi(s, 'LEFT'));
        hasR = any(strcmpi(s, 'RIGHT'));

         if hasL && hasR
             % both LEFT and RIGHT available at timestamp: put together in one run
             runCount = runCount + 1;
             runs{1,runCount} = uTimes{k};
         else
             % only LEFT or only RIGHT: each occurrence is another run
             for j = 1:numel(idx)
                 runCount = runCount + 1;
                 runs{1,runCount} = uTimes{k};
             end
         end
    end

    % ----- above ---- is adjusted !! MJS

    usedCount = containers.Map('KeyType','char','ValueType','double');
    
    % Loop over recordings within file
    for c = 1:length(runs)
        i=perceive_ci(runs{c}, FirstPacketDateTime);

        % ----- below ---- is adjusted !! (in order to collect all runs, even if they are from the same channel on the same timestamp) - MJS

        % --- fix duplicated hits: split or merge depending on channels ---
        if numel(i) > 1
            % ChannelRaw needed for these hits
            ch = {TDdata(i).Channel};  % cell array of chars

            % determine "side" of label after "_":
            tmp  = regexp(ch, '[^_]+$', 'match');
            lab  = cellfun(@(x) x{1}, tmp, 'UniformOutput', false);  % e.g. 'LEFT'/'RIGHT' or other last token

            % CASE A: multiple different labels (eg LEFT+RIGHT) -> put together in one run (keep "i" like it is):
            if numel(unique(lower(lab))) > 1
                % keep i as vector (merge)
            else
                % CASE B: same label (eg RIGHT+RIGHT) -> treat as different runs:
                % select one index from i in an alternating manner, depending on how often this run has occurred
                key = runs{c};  % timestamp string als key

                if ~isKey(usedCount, key)
                    usedCount(key) = 0;
                end
                usedCount(key) = usedCount(key) + 1;

                pick = mod(usedCount(key)-1, numel(i)) + 1;
                i = i(pick);  % make "i" scalar
            end
        end
        % --- end fix ---

        % Save FirstPacketDateTime
        info.FirstPacketDateTime = datetime(FirstPacketDateTime{i(1)});
    
        % Collect ticks and packet sizes of recording
        ticksins = str2num(TDdata(i(1)).TicksInMses)/1000;
        packetsizes = str2num(TDdata(i(1)).GlobalPacketSizes);
    
        % Collect LFP data and check for sample size between channels
        try
            lfp=[TDdata(i).TimeDomainData]';
        catch
            for xi=1:length(i)
                sl(xi)=length(TDdata(i(xi)).TimeDomainData);
            end
            smin=min(sl);
            lfp=[];
            for xi = 1:length(xi)
                lfp(xi,:) = TDdata(i(xi)).TimeDomainData(1:smin);
            end
            warning('Sample size differed between channels. Check session affiliation.')
        end
    
        % Create time array (method from Percept PC whitepaper)
        fs = TDdata.SampleRateInHz;
        seg_missing = 0;
    
        TDtime = ticksins(end) - (packetsizes(end)-1)/fs : 1/fs : ticksins(end);
        for j = length(ticksins)-1:-1:1
            if diff([ticksins(j), ticksins(j+1)]) > (packetsizes(j)+1)/fs
                prev_packet = ticksins(j) - (packetsizes(j)-1)/fs : 1/fs : ticksins(j);
                seg_missing = 1;
            else
                prev_packet = TDtime(1) - packetsizes(j)/fs : 1/fs : TDtime(1) - 1/fs;
            end
            TDtime = [prev_packet, TDtime];
        end

        % Start time array at t=0
        TDtime = TDtime - TDtime(1);

        % Add missing values as NaN
        if seg_missing
            seg_idx = find(round(diff(TDtime),4)-1/fs);
            missing_time = [];
            for k = 1:length(seg_idx)
                missing_time = [missing_time TDtime(seg_idx(k))+1/fs : 1/fs : TDtime(seg_idx(k)+1)-1/fs];
            end
            
            [TDtime, idx_sort] = sort([TDtime missing_time]);
            lfp = [lfp nan(size(lfp,1),length(missing_time))];
            lfp = lfp(:,idx_sort);
        end

        % Match LFP data to stimulation info from BrainSenseLfp
        [~,k] = min(abs(datetime(strrep(strrep({BSdata(:).FirstPacketDateTime},'T',' '),'Z','')) - runs{c}));

        % Retrieve time array and check for missing data
        BStime = [BSdata(k).LfpData.TicksInMs]';
        iM = find(diff(BStime)>mode(diff(BStime)));
        if ~isempty(iM)
            nanTime = [];
            for j = 1:length(iM)
                nanTime = [nanTime BStime(iM(j)) + mode(diff(BStime)) : mode(diff(BStime)) : BStime(iM(j)+1) - mode(diff(BStime))];
            end
            [BStime, iT] = sort([BStime; nanTime']);
        end
        
        % Allign with TD time
        if abs(ticksins(1) - BSdata(k).LfpData(1).TicksInMs/1000) < 1
            BStime = BStime/1000 - ticksins(1);
        elseif isfield(TDdata(i(1)), 'FirstPacketDateTimeOffsetInSeconds')
            BStime = (BStime-BSdata(k).LfpData(1).TicksInMs)/1000 + diff([TDdata(i(1)).FirstPacketDateTimeOffsetInSeconds ...
                      BSdata(k).FirstPacketDateTimeOffsetInSeconds]);
        else
            BStime = (BStime-BSdata(k).LfpData(1).TicksInMs)/1000 + seconds(diff([datetime(strrep(strrep(TDdata(i(1)).FirstPacketDateTime,'T',' '),'Z','')) ...
                      datetime(strrep(strrep(BSdata(k).FirstPacketDateTime,'T',' '),'Z',''))]));
        end

        % Retrieve stimulation amplitudes and LFP power values
        stim = {};
        stim.Time = BStime;
        stim.Amplitude = nan(length(BSdata(k).LfpData),2);
        
        power = {};
        power.Time = BStime;
        power.Value = nan(length(BSdata(k).LfpData),2);

        if isfield(BSdata(k).LfpData, 'Left')
            for s = 1:length(BSdata(k).LfpData)
                if isfield(BSdata(k).LfpData(s).Left, 'mA')
                    stim.Amplitude(s,1) = BSdata(k).LfpData(s).Left.mA;
                end
                if isfield(BSdata(k).LfpData(s).Left, 'LFP')  && contains([TDdata(i(:)).Channel], 'LEFT')
                    power.Value(s,1) = BSdata(k).LfpData(s).Left.LFP;
                end
            end
        end
        if isfield(BSdata(k).LfpData, 'Right')
            for s = 1:length(BSdata(k).LfpData)
                if isfield(BSdata(k).LfpData(s).Right, 'mA')
                    stim.Amplitude(s,2) = BSdata(k).LfpData(s).Right.mA;
                end
                if isfield(BSdata(k).LfpData(s).Right, 'LFP') && contains([TDdata(i(:)).Channel], 'RIGHT')
                    power.Value(s,2) = BSdata(k).LfpData(s).Right.LFP;
                end
            end
        end
        if ~isempty(iM)
            stim.Amplitude = [stim.Amplitude; nan(length(nanTime),2)];
            stim.Amplitude = stim.Amplitude(iT,:);
            power.Value = [power.Value; nan(length(nanTime),2)];
            power.Value = power.Value(iT,:);
        end
        power.Value = power.Value(:,~all(isnan(power.Value)));

        % Retrieve settings
        settings = {};
        settings.StimulationFrequency = nan(1,2);
        settings.PulseWidth = nan(1,2);
        settings.SensingChannel = {'NOT_CONFIGURED' 'NOT_CONFIGURED'};
        settings.SensingFrequency = nan(1,2);
        settings.SamplingFrequencyTimeDomain = fs;
        settings.SamplingFrequencyPower = BSdata(k).SampleRateInHz;
        settings.LowerLfpThreshold = nan(1,2);
        settings.UpperLfpThreshold = nan(1,2);
        settings.LowerStimulationLimit = nan(1,2);
        settings.UpperStimulationLimit = nan(1,2);
        settings.AdaptiveMode = {'NOT_CONFIGURED' 'NOT_CONFIGURED'};
        settings.AdaptiveStatus = {'NOT_CONFIGURED' 'NOT_CONFIGURED'};

        if isfield(BSdata(k).TherapySnapshot, 'Left')
            settings.StimulationFrequency(1) = BSdata(k).TherapySnapshot.Left.RateInHertz;
            settings.PulseWidth(1) = BSdata(k).TherapySnapshot.Left.PulseWidthInMicroSecond;
            settings.SensingFrequency(1) = BSdata(k).TherapySnapshot.Left.FrequencyInHertz;
            settings.LowerLfpThreshold(1) = BSdata(k).TherapySnapshot.Left.LowerLfpThreshold;
            settings.UpperLfpThreshold(1) = BSdata(k).TherapySnapshot.Left.UpperLfpThreshold;
            settings.LowerStimulationLimit(1) = BSdata(k).TherapySnapshot.Left.LowerLimitInMilliAmps;
            settings.UpperStimulationLimit(1) = BSdata(k).TherapySnapshot.Left.UpperLimitInMilliAmps;
            chan = split(BSdata(k).TherapySnapshot.Left.SensingChannel,'.');
            settings.SensingChannel{1} = chan{2};

            if isfield(BSdata(k).TherapySnapshot.Left, 'StreamingAdaptiveMode')
                aMode = split(BSdata(k).TherapySnapshot.Left.StreamingAdaptiveMode,'.');
                settings.AdaptiveMode{1} = aMode{2};
            end
            if isfield(BSdata(k).TherapySnapshot.Left, 'AdaptiveTherapyStatus')
                status = split(BSdata(k).TherapySnapshot.Left.AdaptiveTherapyStatus,'.');
                settings.AdaptiveStatus{1} = status{2};
            end
        end
        if isfield(BSdata(k).TherapySnapshot, 'Right')
            settings.StimulationFrequency(2) = BSdata(k).TherapySnapshot.Right.RateInHertz;
            settings.PulseWidth(2) = BSdata(k).TherapySnapshot.Right.PulseWidthInMicroSecond;
            settings.SensingFrequency(2) = BSdata(k).TherapySnapshot.Right.FrequencyInHertz;
            settings.LowerLfpThreshold(2) = BSdata(k).TherapySnapshot.Right.LowerLfpThreshold;
            settings.UpperLfpThreshold(2) = BSdata(k).TherapySnapshot.Right.UpperLfpThreshold;
            settings.LowerStimulationLimit(2) = BSdata(k).TherapySnapshot.Right.LowerLimitInMilliAmps;
            settings.UpperStimulationLimit(2) = BSdata(k).TherapySnapshot.Right.UpperLimitInMilliAmps;
            chan = split(BSdata(k).TherapySnapshot.Right.SensingChannel,'.');
            settings.SensingChannel{2} = chan{2};

            if isfield(BSdata(k).TherapySnapshot.Right, 'StreamingAdaptiveMode')
                aMode = split(BSdata(k).TherapySnapshot.Right.StreamingAdaptiveMode,'.');
                settings.AdaptiveMode{2} = aMode{2};
            end
            if isfield(BSdata(k).TherapySnapshot.Right, 'AdaptiveTherapyStatus')
                status = split(BSdata(k).TherapySnapshot.Right.AdaptiveTherapyStatus,'.');
                settings.AdaptiveStatus{2} = status{2};
            end
        end
    
        % Add info and data to one struct
        dataRec = {};
        dataRec.DataType = 'BrainSenseTimeDomain/BrainSenseLfp';
        dataRec.Run = c;
        dataRec.DateTime = char(info.FirstPacketDateTime + hours(str2double(info.ProgrammerUtcOffset(1:3))) ...
                                + minutes(str2double(info.ProgrammerUtcOffset([1, 5:6]))));
        dataRec.UtcOffset = info.ProgrammerUtcOffset;
        dataRec.Info = info;
        dataRec.Settings = settings;
        dataRec.Data = {};
        dataRec.Data.LfpPower = power;
        dataRec.Data.Stimulation = stim;
        dataRec.Data.Channel = {TDdata(i(:)).Channel};
        dataRec.Data.Time = TDtime';
        dataRec.Data.LFP = lfp';
        dataRec.RawPower = BSdata(k);
        dataRec.RawTimeDomain = TDdata(i);
        dataRec.EcgMethod = ecgMethod;

        % Apply Linenoise notch and ECG suppression
        if ecgMethod == 0
            dataRec = filterLinenoise(dataRec, linenoise);
            dataRec.Data.LFP_Ecg = NaN(size(dataRec.Data.LFP,1), size(dataRec.Data.LFP,2));
            dataRec.Data.LFP_LinenoiseEcg = NaN(size(dataRec.Data.LFP,1), size(dataRec.Data.LFP,2));
        elseif ecgMethod == 1
            dataRec = filterLinenoise(dataRec, linenoise);
            dataRec = filterECG(dataRec, rTime);
        elseif ecgMethod == 2
            dataRec = filterLinenoise(dataRec, linenoise);
            dataRec = filterECGfull(dataRec, rTime, savepath, savenameJSON, c);
        end

        % Plot data if applicable
        if plotData
            plotStreaming(dataRec, length(runs), ecgMethod, [savepath filesep 'Figures'], [savenameJSON '_Streaming_rec' num2str(c)], showFig);
        end

        % Save all recordings into one struct, which can be given as output if needed
        dataStreaming(c) = dataRec;
    end

    % Save data
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    save([savepath filesep savenameJSON '_Streaming'], 'dataStreaming')

end
