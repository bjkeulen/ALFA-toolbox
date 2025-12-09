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

%% Call user interface for settings
addpath(genpath(cd))
settings = userInterface();
if isempty(fieldnames(settings))
    error('Analysis canceled.')
end
            
%% Get files and directories

% Get single .json file of patient
if settings.dataset == 0
    [fname, fldr] = uigetfile('*.json', 'Select a single JSON file');

    if ischar(fname) && ischar(fldr)
        if strcmp(fname(end-4:end),'.json')
            savepath = [fldr filesep fname(1:end-5) '_unpacked'];
            nfldr = 1;
            nfiles = 1;
        else
            error('No JSON file selected.')
        end
    else
        error('No file selected.')
    end

% Get folder and all .json files within folder
elseif settings.dataset == 1
    rootFolder = uigetdir('', 'Select a folder');

    if ischar(rootFolder)
        files = dir([rootFolder '\**\*.json']);

        % Check if files are found, remove duplicates if applicable
        if isempty(files)
            error('No JSON files found in selected folder.')
        else
            if length(files) > 1
                fileTable = struct2table(files);
                [~,idx,~] = unique(fileTable.name);
                files = files(idx,:);
            end
            nfiles = length(files);
            nfldr = 1;

            [~,rootName] = fileparts(rootFolder);
            savepath = [rootFolder filesep 'folder_' rootName '_unpacked'];
        end
    else
        error('No folder selected.')
    end

% Get folder and subfolders
elseif settings.dataset == 2
    rootFolder = uigetdir('', 'Select a folderset');
    settings.showFig = 0;

    if ischar(rootFolder)
        fldrs = dir(rootFolder);
        fldrs = fldrs([fldrs.isdir],:);
        fldrs = fldrs(~matches({fldrs.name},[".",".."]),:);
        nfldr = length(fldrs);

        % Check if files are found, remove duplicates if applicable
        if isempty(fldrs)
            error('No subfolders found. Select a folder containing subfolders')
        else
            filecount = [];
            for c = 1:length(fldrs)
                filecount = [filecount; dir([rootFolder filesep fldrs(c).name '\**\*.json'])];
            end

            if isempty(filecount)
                error('No JSON files found in subfolders of selected folderset.')
            else
                [~,rootName] = fileparts(rootFolder);
                savepath = [rootFolder filesep 'folderset_' rootName '_unpacked'];
            end
        end
    else
        error('No folder selected.')
    end
end

%% Settings and initialisation of structures

% Define which data to extract from json
datafields = {'BrainSenseTimeDomain','DiagnosticData','LfpMontageTimeDomain','BrainSenseSurveysTimeDomain','SenseChannelTests'};

% Create folder to save data and plots into, if needed
if ~exist(savepath, 'dir')
   mkdir(savepath)
end

% Define fieldnames for Events and Timeline
fieldsTL = {'DateTime','ActiveGroup','SensingChannel','SensingFrequency','LfpPower','LowerLfpThreshold','UpperLfpThreshold','StimulationAmplitude','LowerStimulationLimit','UpperStimulationLimit','StimulationFrequency','PulseWidth'};
fieldsSE = {'DateTime','EventName','EventID','SensingChannel','Frequency','Magnitude','StimulationAmplitude','StimulationFrequency','PulseWidth'};

%% Loop over folders
for f = 1:nfldr

    % Get files within directory in case of settings.dataset of folders
    if settings.dataset == 2
        files = dir([rootFolder '\' fldrs(f).name '\**\*.json']);
        if isempty(files)
            continue
        else
            fileTable = struct2table(files);
            if height(fileTable) == 1
                [~,idx,~] = unique({fileTable.name});
            else
                [~,idx,~] = unique(fileTable.name);
            end
            files = files(idx,:);
            nfiles = length(files);
        end
    end

    % Initialise structure for Timeline data
    dataTimeline = struct;
    dataTimeline.DataType = 'LFPTrendLogs';
    dataTimeline.Info = [];
    dataTimeline.Labels = {'LEFT','RIGHT'};
    for n = 1:length(fieldsTL)
        dataTimeline.Data.(fieldsTL{n}) = [];
    end
    dataTimeline.History = struct;
    
    % Initialise structure for Events data
    dataEvents = struct;
    dataEvents.DataType = 'LfpFrequencySnapshotEvents';
    dataEvents.Info = [];
    dataEvents.Labels = {'LEFT', 'RIGHT'};
    dataEvents.Data = struct;
    for n = 1:length(fieldsSE)
        dataEvents.Data.(fieldsSE{n}) = [];
    end

    %% Loop over files
    for i = 1:nfiles

        % Retrieve file and folder; define name to save data under
        if settings.dataset == 0
            savename_json = fname(1:end-5);
        elseif settings.dataset == 1
            fldr = files(i).folder;
            fname = files(i).name;
            savename_json = [rootName '_json' num2str(i) '_' fname(1:end-5)];
        elseif settings.dataset == 2
            fldr = files(i).folder;
            fname = files(i).name;
            savename_json = [fldrs(f).name '_json' num2str(i) '_' fname(1:end-5)];
        end

        % Load JSON
        js = jsondecode(fileread([fldr filesep fname]));

        % Collect relevant information
        info = {};
        info.OriginalFolder = fldr;
        info.OriginalFile = fname;
        info.SessionStartDateUtc = datetime(strrep(js.SessionDate(1:end-1),'T',' '));
        info.SessionEndDateUtc = datetime(strrep(js.SessionEndDate(1:end-1),'T',' '));
        info.ProgrammerUtcOffset = js.ProgrammerUtcOffset;
        info.GroupsInitial = js.Groups.Initial;
        info.GroupsFinal = js.Groups.Final;
        if isfield(js, 'GroupHistory')
            info.GroupHistory = js.GroupHistory;
        else
            info.GroupHistory = {};
        end
        if isfield(js, 'DiagnosticData')
            info.EventLogs = js.DiagnosticData.EventLogs;
        else
            info.EventLogs = {};
        end
    
        % Always add info to dataTimeline, for first sensing info
        dataTimeline.Info = vertcat(dataTimeline.Info, info);
    
        % Check data types in file and loop over data types present
        for b = 1:length(datafields)
            if isfield(js,datafields{b})
                data = js.(datafields{b});
                if isempty(data)
                    continue
                end
        
                switch datafields{b}
                        
                    % % Streaming
                    case 'BrainSenseTimeDomain'

                        % Define savepath
                        if settings.dataset == 2
                            savepath_streaming = [savepath filesep 'Streaming'];
                            if ~exist(savepath_streaming, 'dir')
                               mkdir(savepath_streaming)
                            end
                        else
                            savepath_streaming = savepath;
                        end

                        % Get data and plot if applicable
                        dataStreaming = getStreaming(info, js, settings.linenoise, settings.ecgMethod, settings.rWindow, savepath_streaming, savename_json, settings.plotData, settings.showFig);

                    % Setup
                    case 'SenseChannelTests'

                        % Define savepath
                        if settings.dataset == 2
                            savepath_setup = [savepath filesep 'Setup'];
                            if ~exist(savepath_setup, 'dir')
                               mkdir(savepath_setup)
                            end
                        else
                            savepath_setup = savepath;
                        end

                        % Get data and plot if applicable
                        dataSetup = getSetup(info, js, savepath_setup, savename_json);
                        if settings.plotData
                            plotSetup(dataSetup, [savepath_setup filesep 'Figures'], savename_json, settings.showFig);
                        end

                    % Survey old format
                    case 'LfpMontageTimeDomain'

                        % Only if new format is not present
                        if ~isfield(js, 'BrainSenseSurveysTimeDomain')

                            % Define savepath
                            if settings.dataset == 2
                                savepath_survey = [savepath filesep 'Survey'];
                                if ~exist(savepath_survey, 'dir')
                                   mkdir(savepath_survey)
                                end
                            else
                                savepath_survey = savepath;
                            end

                            % Get data and plot if applicable
                            dataSurvey = getSurvey(info, js, [], savepath_survey, savename_json);
                            if settings.plotData
                                plotSurveys(dataSurvey, [savepath_survey filesep 'Figures'], savename_json, settings.showFig);
                            end
                        end

                    % Survey and Identifier
                    case 'BrainSenseSurveysTimeDomain'

                        % Loop over fields and check type of Survey
                        for idx = 1:length(js.BrainSenseSurveysTimeDomain)
                            if isfield(js.BrainSenseSurveysTimeDomain{idx}, 'ElectrodeSurvey')
                                if ~isempty(vertcat(js.BrainSenseSurveysTimeDomain{idx}.ElectrodeSurvey.TimeDomainDatainMicroVolts))

                                    % Define savepath
                                    if settings.dataset == 2
                                        savepath_survey = [savepath filesep 'Survey'];
                                        if ~exist(savepath_survey, 'dir')
                                           mkdir(savepath_survey)
                                        end
                                    else
                                        savepath_survey = savepath;
                                    end

                                    % Get data and plot if applicable
                                    dataSurvey = getSurvey(info, js, idx, savepath_survey, savename_json);
                                    if settings.plotData
                                        plotSurveys(dataSurvey, [savepath_survey filesep 'Figures'], savename_json, settings.showFig);
                                    end
                                end

                            elseif isfield(js.BrainSenseSurveysTimeDomain{idx}, 'ElectrodeIdentifier')
                                if ~isempty(vertcat(js.BrainSenseSurveysTimeDomain{idx}.ElectrodeIdentifier.TimeDomainDatainMicroVolts))

                                    % Define savepath
                                    if settings.dataset == 2
                                        savepath_identifier = [savepath filesep 'Identifier'];
                                        if ~exist(savepath_identifier, 'dir')
                                           mkdir(savepath_identifier)
                                        end
                                    else
                                        savepath_identifier = savepath;
                                    end

                                    % Get data and plot if applicable
                                    dataIdentifier = getIdentifier(info, js, idx, savepath_identifier, savename_json);
                                    if settings.plotData
                                        plotSurveys(dataIdentifier, [savepath_identifier filesep 'Figures'], savename_json, settings.showFig);
                                    end
                                end
                            end
                        end

                    case 'DiagnosticData'

                        % Timeline
                        if isfield(data,'LFPTrendLogs')

                            % Retrieve data from recording and concatenate with all data
                            timeline = getTimeline(data, fieldsTL);
                            for n = 1:length(fieldsTL)
                                dataTimeline.Data.(fieldsTL{n}) = vertcat(dataTimeline.Data.(fieldsTL{n}), timeline.(fieldsTL{n}));
                            end
                        end

                        % Events
                        if isfield(data,'LfpFrequencySnapshotEvents')

                            % Retrieve data from recording and concatenate with all data
                            events = getEvents(data.LfpFrequencySnapshotEvents, fieldsSE);
                            dataEvents.Info = vertcat(dataEvents.Info, info);
                            dataEvents.Data = horzcat(dataEvents.Data, events);
                        end
                end
            end
        end
    end

    %% Further processing Timeline data
    if ~isempty(dataTimeline.Data.LfpPower)

        % Define savename
        if settings.dataset == 0
            savename = fname(1:end-5);
        elseif settings.dataset == 1
            savename = rootName;
        elseif settings.dataset == 2
            savename = fldrs(f).name;
        end
        
        % Further processing
        dataTimeline.Data = sortrows(struct2table(dataTimeline.Data),'DateTime');     % Convert to table and sort on datetime       
        dataTimeline = checkDuplicates(dataTimeline, fieldsTL);                       % Remove duplicates
        dataTimeline = addNaN(dataTimeline);                                          % Insert NaN for missing values
        dataTimeline = addParameters(dataTimeline, tZone);                            % Get sensing and stimulation parameters and assign to data

        % Change timezone to correct for UTC offset
        dataTimeline.Data.DateTime.TimeZone = 'UTC';
        dataTimeline.Data.DateTime.TimeZone = tZone;
        dataTimeline.TimeZone = tZone;
    
        % Save data
        if settings.dataset == 0
            save([savepath filesep savename '_Timeline.mat'], 'dataTimeline')
        elseif settings.dataset == 1
            save([savepath filesep savename '_Timeline.mat'], 'dataTimeline')
        elseif settings.dataset == 2
            savepath_TL = [savepath filesep 'Timeline'];
            if ~exist(savepath_TL, 'dir')
               mkdir(savepath_TL)
            end
            save([savepath_TL filesep savename '_Timeline.mat'], 'dataTimeline')
        end
    
        % Plot data
        if plotData == 1
            if settings.dataset == 0
                plotTimeline(dataTimeline, [savepath filesep 'Figures'], savename, showFig)
            elseif settings.dataset == 1
                plotTimeline(dataTimeline, [savepath filesep 'Figures'], savename, showFig)
            elseif settings.dataset == 2
                plotTimeline(dataTimeline, [savepath_TL filesep 'Figures'], savename, showFig)
            end
        end
    end
    
    %% Further processing Events data
    if ~isempty([dataEvents.Data.EventName])

        % Further processing
        dataEvents.Data = dataEvents.Data(2:end);                         % Remove first empty row
        if length(dataEvents.Data) > 1                                    % Convert to table and sort on datetime
            dataEvents.Data = sortrows(struct2table(dataEvents.Data),'DateTime');     
        end
        dataEvents = checkDuplicates(dataEvents);                         % Remove duplicates

        % Change timezone to correct for UTC offset
        dataEvents.Data.DateTime.TimeZone = 'UTC';
        dataEvents.Data.DateTime.TimeZone = tZone;
        dataEvents.TimeZone = tZone;

        % Add stimulation amplitude, pulsewidth and frequency
        dataEvents = addStimInfo(dataEvents, dataTimeline);

        % Save data
        if settings.dataset == 0
            save([savepath filesep fname(1:end-5) '_Events.mat'], 'dataEvents')
        elseif settings.dataset == 1
            save([savepath filesep rootName '_Events.mat'], 'dataEvents')
        elseif settings.dataset == 2
            savepath_SE = [savepath filesep 'Events'];
            if ~exist(savepath_SE, 'dir')
               mkdir(savepath_SE)
            end
            save([savepath_SE filesep fldrs(f).name '_Events.mat'], 'dataEvents')
        end

        % Plot data
        if plotData == 1
            if settings.dataset == 0
                plotEvents(dataEvents, [savepath filesep 'Figures'], fname(1:end-5), showFig)
            elseif settings.dataset == 1
                plotEvents(dataEvents, [savepath filesep 'Figures'], rootName, showFig)
            elseif settings.dataset == 2
                plotEvents(dataEvents, [savepath_SE filesep 'Figures'], fldrs(f).name, showFig)
            end
        end
    end
end