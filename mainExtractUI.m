% Author: B.J. Keulen
% Date: 05-03-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
% 
% GitHub Repository: https://github.com/bjkeulen/ALFA-toolbox.
%
% Main script of the Amsterdam Local Field potential Analysis (ALFA) 
% toolbox created by B.J. Keulen and M.J. Stam. The ALFA toolbox can be 
% used for extracting, processing, rearranging and visualizing data from 
% Medtronic PerceptTM JSON reports containing local field potential (LFP) 
% data. General settings are given through a user interface which appears
% when running this script.
%
% Part of the structure of this code is adapted from the perceive.m file
% of the perceive toolbox by Wolf-Julian Neumann, Tomas Sieger and Gerd 
% Tinkhauser: https://github.com/neuromodulation/perceive.
%
% Other parts of the code within this file are inspired by the Percept 
% toolbox by Yohann Thenaisie and Bart Keulen:
% https://github.com/YohannThenaisie/PerceptToolbox.
%
% This is an open research tool that is not intended for clinical purposes.

%% Settings
datafields = {'BrainSenseTimeDomain','DiagnosticData','LfpMontageTimeDomain','BrainSenseSurveysTimeDomain','SenseChannelTests'};
fieldsTL = {'DateTime','ActiveGroup','SensingChannel','SensingFrequency','LfpPower','LowerLfpThreshold','UpperLfpThreshold','StimulationAmplitude','LowerStimulationLimit','UpperStimulationLimit','StimulationFrequency','PulseWidth'};
fieldsSE = {'DateTime','EventName','EventID','SensingChannel','Frequency','Magnitude','StimulationAmplitude','StimulationFrequency','PulseWidth'};

% Check if underlying folders are added to path, otherwise do so
if ~exist('mainPath','var')
    toolboxPath = pwd;
    addpath(genpath(toolboxPath))
end

% Call user interface for settings
settings = userInterface();
if isempty(fieldnames(settings))
    error('Analysis canceled.')
end
            
% Get files and directories
[fileData, savepath] = getFileData(settings.dataset);
if settings.dataset == 2
    settings.showFig = 0;
end

%% Loop over folders
for f = 1:fileData.nfolders

    % Get files within directory in case of folderset processing
    if settings.dataset == 2
        fileData.files = dir([fileData.rootPath filesep fileData.folders(f).name filesep '**' filesep '*.json']);
        if isempty(fileData.files)
            continue
        else
            fileTable = struct2table(fileData.files);
            if height(fileTable) == 1
                [~,idx,~] = unique({fileTable.name});
            else
                [~,idx,~] = unique(fileTable.name);
            end
            fileData.files = fileData.files(idx,:);
            fileData.nfiles = length(fileData.files);
        end
    end

    % Get empty output structures for Timeline and Events data
    dataTimeline = getStructures('Timeline', fieldsTL);
    dataEvents = getStructures('Events', fieldsSE);
    dataLog = table('Size',[fileData.nfiles,8], 'VariableTypes',{'string','string','int64','int64','int64','int64','int64','int64'}, ...
                     'VariableNames',{'Filename','File_check','Setup_n_runs','Survey_n_runs','Identifier_n_runs','Streaming_n_runs','Timeline_present','Events_present'});

    %% Loop over files
    for i = 1:fileData.nfiles

        % Retrieve file and folder; define name to save data under
        if settings.dataset == 0
            folder = fileData.rootPath;
            filename = fileData.files;
            savenameJSON = filename(1:end-5);
        elseif settings.dataset == 1
            folder = fileData.files(i).folder;
            filename = fileData.files(i).name;
            savenameJSON = [fileData.rootName '_json' num2str(i) '_' filename(1:end-5)];
        elseif settings.dataset == 2
            folder = fileData.files(i).folder;
            filename = fileData.files(i).name;
            savenameJSON = [fileData.folders(f).name '_json' num2str(i) '_' filename(1:end-5)];
        end

        % Add filename to log
        dataLog{i,"Filename"} = string(filename);

        % Load JSON, check format, and collect general info
        try
            js = jsondecode(fileread([folder filesep filename]));
            try
                datetime(strrep(js.SessionDate(1:end-1),'T',' '));
            catch
                warning('Could not parse SessionDate for file: %s', filename);
                dataLog{i,'File_check'} = "Incorrect format SessionDate";
                continue
            end
            try
                datetime(strrep(js.SessionEndDate(1:end-1),'T',' '));
            catch
                warning('Could not parse SessionEndDate for file: %s', filename);
                dataLog{i,'File_check'} = "Incorrect format SessionEndDate";
                dataLog{i,3:end} = -1;
                continue
            end
            dataLog{i,'File_check'} = "Pass";
        catch
            warning('Could not load file: %s', filename);
            dataLog{i,'File_check'} = "Could not load file";
            dataLog{i,3:end} = -1;
            continue
        end
        info = getInfo(folder, filename, js);
        dataTimeline.Info = vertcat(dataTimeline.Info, info);
    
        % Check data types in file and loop over data types present
        for b = 1:length(datafields)
            if isfield(js,datafields{b})
                data = js.(datafields{b});
                if isempty(data)
                    continue
                end
        
                switch datafields{b}

                    % Setup
                    case 'SenseChannelTests'

                        % Folderset processing
                        if settings.dataset == 2
                            savepath_setup = [savepath filesep 'Setup'];
                            if ~exist([savepath filesep 'Setup'], 'dir')
                               mkdir([savepath filesep 'Setup'])
                            end
                            dataSetup = getSetup(info, js, [savepath filesep 'Setup'], savenameJSON);
                            if settings.plotData
                                plotSetup(dataSetup, [savepath filesep 'Setup' filesep 'Figures'], savenameJSON, settings.showFig);
                            end

                        % Single file and folder processing
                        else
                            dataSetup = getSetup(info, js, savepath, savenameJSON);
                            if settings.plotData
                                plotSetup(dataSetup, [savepath filesep 'Figures'], savenameJSON, settings.showFig);
                            end
                        end

                        % Add number of runs to log
                        dataLog{i,'Setup_n_runs'} = dataSetup.nRuns;

                    % Survey old format
                    case 'LfpMontageTimeDomain'

                        % Only if new format is not present
                        if ~isfield(js, 'BrainSenseSurveysTimeDomain')

                            % Folderset processing
                            if settings.dataset == 2
                                if ~exist([savepath filesep 'Survey'], 'dir')
                                   mkdir([savepath filesep 'Survey'])
                                end
                                dataSurvey = getSurvey(info, js, [], [savepath filesep 'Survey'], savenameJSON);
                                if settings.plotData
                                    plotSurveys(dataSurvey, [savepath filesep 'Survey' filesep 'Figures'], savenameJSON, settings.showFig);
                                end

                            % Single file and folder processing
                            else
                                dataSurvey = getSurvey(info, js, [], savepath, savenameJSON);
                                if settings.plotData
                                    plotSurveys(dataSurvey, [savepath filesep 'Figures'], savenameJSON, settings.showFig);
                                end
                            end

                            % Add number of runs to log
                            dataLog{i,'Survey_n_runs'} = dataSurvey.nRuns;
                        end

                    % Survey and Identifier
                    case 'BrainSenseSurveysTimeDomain'

                        % Loop over fields and check type of Survey
                        for idx = 1:length(js.BrainSenseSurveysTimeDomain)
                            if isfield(js.BrainSenseSurveysTimeDomain{idx}, 'ElectrodeSurvey')
                                if ~isempty(vertcat(js.BrainSenseSurveysTimeDomain{idx}.ElectrodeSurvey.TimeDomainDatainMicroVolts))

                                    % Folderset processing
                                    if settings.dataset == 2
                                        savepath_survey = [savepath filesep 'Survey'];
                                        if ~exist([savepath filesep 'Survey'], 'dir')
                                           mkdir([savepath filesep 'Survey'])
                                        end
                                        dataSurvey = getSurvey(info, js, idx, [savepath filesep 'Survey'], savenameJSON);
                                        if settings.plotData
                                            plotSurveys(dataSurvey, [savepath filesep 'Survey' filesep 'Figures'], savenameJSON, settings.showFig);
                                        end

                                    % Single file and folder processing
                                    else
                                        dataSurvey = getSurvey(info, js, idx, savepath, savenameJSON);
                                        if settings.plotData
                                            plotSurveys(dataSurvey, [savepath filesep 'Figures'], savenameJSON, settings.showFig);
                                        end
                                    end

                                    % Add number of runs to log
                                    dataLog{i,'Survey_n_runs'} = dataSurvey.nRuns;
                                end

                            elseif isfield(js.BrainSenseSurveysTimeDomain{idx}, 'ElectrodeIdentifier')
                                if ~isempty(vertcat(js.BrainSenseSurveysTimeDomain{idx}.ElectrodeIdentifier.TimeDomainDatainMicroVolts))

                                    % Folderset processing
                                    if settings.dataset == 2
                                        savepath_identifier = [savepath filesep 'Identifier'];
                                        if ~exist([savepath filesep 'Identifier'], 'dir')
                                           mkdir([savepath filesep 'Identifier'])
                                        end
                                        dataIdentifier = getIdentifier(info, js, idx, [savepath filesep 'Identifier'], savenameJSON);
                                        if settings.plotData
                                            plotSurveys(dataIdentifier, [savepath filesep 'Identifier' filesep 'Figures'], savenameJSON, settings.showFig);
                                        end

                                    % Single file and folder processing
                                    else
                                        dataIdentifier = getIdentifier(info, js, idx, savepath, savenameJSON);
                                        if settings.plotData
                                            plotSurveys(dataIdentifier, [savepath filesep 'Figures'], savenameJSON, settings.showFig);
                                        end
                                    end

                                    % Add number of runs to log
                                    dataLog{i,'Identifier_n_runs'} = dataIdentifier.nRuns;
                                end
                            end
                        end

                    % Streaming
                    case 'BrainSenseTimeDomain'

                        % Folderset processing
                        if settings.dataset == 2
                            if ~exist([savepath filesep 'Streaming'], 'dir')
                               mkdir([savepath filesep 'Streaming'])
                            end
                            dataStreaming = getStreaming(info, js, settings.linenoise, settings.ecgMethod, settings.rWindow, [savepath filesep 'Streaming'], savenameJSON, settings.plotData, settings.showFig);

                        % Single file and folder processing
                        else
                            dataStreaming = getStreaming(info, js, settings.linenoise, settings.ecgMethod, settings.rWindow, savepath, savenameJSON, settings.plotData, settings.showFig);
                        end

                        % Add number of runs to log
                        dataLog{i,'Streaming_n_runs'} = length(dataStreaming);

                    case 'DiagnosticData'

                        % Timeline
                        if isfield(data,'LFPTrendLogs')

                            % Retrieve data from recording and concatenate with all data
                            timeline = getTimeline(data, fieldsTL);
                            for n = 1:length(fieldsTL)
                                dataTimeline.Data.(fieldsTL{n}) = vertcat(dataTimeline.Data.(fieldsTL{n}), timeline.(fieldsTL{n}));
                            end
                            dataLog{i,'Timeline_present'} = 1;
                        end

                        % Events
                        if isfield(data,'LfpFrequencySnapshotEvents')

                            % Retrieve data from recording and concatenate with all data
                            events = getEvents(data.LfpFrequencySnapshotEvents, fieldsSE);
                            dataEvents.Info = vertcat(dataEvents.Info, info);
                            dataEvents.Data = horzcat(dataEvents.Data, events);
                            dataLog{i,'Events_present'} = 1;
                        end
                end
            end
        end
    end

    %% Process, save and plot aggregated Timeline data
    if ~isempty(dataTimeline.Data.LfpPower)
        dataTimeline = processAggregated(dataTimeline, {}, fieldsTL, settings.dataset, fileData, settings.tZone, savepath, filename, f);
        if settings.plotData == 1
            if settings.dataset == 0
                plotTimeline(dataTimeline, [savepath filesep 'Figures'], filename(1:end-5), settings.showFig)
            elseif settings.dataset == 1
                plotTimeline(dataTimeline, [savepath filesep 'Figures'], fileData.rootName, settings.showFig)
            elseif settings.dataset == 2
                plotTimeline(dataTimeline, [savepath filesep 'Timeline' filesep 'Figures'], fileData.folders(f).name, settings.showFig)
            end
        end
    end

    %% Process, save and plot aggregated Events data
    if ~isempty([dataEvents.Data.EventName])
        dataEvents = processAggregated(dataEvents, dataTimeline, {}, settings.dataset, fileData, settings.tZone, savepath, filename, f);
        if settings.plotData == 1
            if settings.dataset == 0
                plotEvents(dataEvents, [savepath filesep 'Figures'], filename(1:end-5), settings.showFig)
            elseif settings.dataset == 1
                plotEvents(dataEvents, [savepath filesep 'Figures'], fileData.rootName, settings.showFig)
            elseif settings.dataset == 2
                plotEvents(dataEvents, [savepath filesep 'Events' filesep 'Figures'], fileData.folders(f).name, settings.showFig)
            end
        end
    end

    %% Complete and save datalog
    if settings.dataset == 0
        saveLogs(dataLog, dataTimeline, dataEvents, savepath, filename(1:end-5), [])
    elseif settings.dataset == 1
        saveLogs(dataLog, dataTimeline, dataEvents, savepath, fileData.rootName, [])
    elseif settings.dataset == 2
        saveLogs(dataLog, dataTimeline, dataEvents, savepath, fileData.folders(f).name, fileData.rootName)
    end
end