% Author: B.J. Keulen
% Date: 21-10-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
% 
% GitHub Repository: https://github.com/bjkeulen/ALFA-toolbox.
%
% Main script of the Amsterdam Local Field potential Analysis (ALFA) 
% toolbox created by B.J. Keulen and M.J. Stam. This main script can be 
% used for plotting outputs from the ALFA toolbox (see mainExtract.m and 
% mainExtractUI.m). One or more files can be selected. Possible inputs are 
% MAT files containing either dataStreaming, dataTimeline, dataEvents,
% dataSurvey, dataIdentifier or dataSetup. Figures are created and saved 
% into a subfolder of the folder containing the MAT file(s).
%
% This is an open research tool that is not intended for clinical purposes.

%% Check if underlying folders are added to path, otherwise do so
if ~exist('mainPath','var')
    mainPath = pwd;
    addpath(genpath(mainPath))
end

%% Settings
folder = 1;      % Select either a single file (0) or a folder (1)
showFig = 1;     % Show figures (1) or not (0)

%% Clear existing output variables if needed
if exist('dataStreaming', 'var')
    clear dataStreaming
end
if exist('dataSurvey', 'var')
    clear dataSurvey
end
if exist('dataIdentifier', 'var')
    clear dataIdentifier
end
if exist('dataSetup', 'var')
    clear dataSetup
end
if exist('dataTimeline', 'var')
    clear dataTimeline
end
if exist('dataEvents', 'var')
    clear dataEvents
end

%% Get files and directories

% Get single .json file of patient
if folder == 0
    [file, rootFolder] = uigetfile('*.mat', 'Select a single MAT file');

    if ischar(file) && ischar(rootFolder)
        files = {};
        files.name = file;
        files.folder = rootFolder;
        savepath = [rootFolder filesep 'Figures'];
    else
        error('No file selected. Try again.')
    end

% Get folder and all .json files within folder
else
    rootFolder = uigetdir('', 'Select a folder containing MAT files');

    if ischar(rootFolder)
        files = dir([rootFolder filesep '**' filesep '*.mat']);

        % Check if files are found, remove duplicates if applicable
        if isempty(files)
            error('No files found in selected folder. Select a folder containing MAT files.')
        elseif length(files) > 1
            fileTable = struct2table(files);
            [~,idx,~] = unique(fileTable.name);
            files = files(idx,:);
            savepath = [rootFolder filesep 'Figures'];
        end
    else
        error('No folder selected. Try again.')
    end
end

%% Loop over files and load data
for i = 1:size(files,1)
    load([files(i).folder filesep files(i).name]);
    fsave = split(files(i).name,"_");
    fsave = strjoin(fsave(1:end-1),"_");

    % Streaming
    if exist('dataStreaming', 'var')
        for r = 1:length(dataStreaming)
            try 
                dataRec = dataStreaming(r);
            catch
                error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
            end

            % Loop over recordings in data
            if strcmp(dataRec.DataType, 'BrainSenseTimeDomain/BrainSenseLfp')
                plotStreaming(dataRec, length(dataStreaming), dataRec.EcgMethod, savepath, [fsave '_Streaming_rec' num2str(r)], showFig)
                clear dataRec
            else
                error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
            end
        end
        clear dataStreaming
        
    % Survey
    elseif exist('dataSurvey', 'var')
        if strcmp(dataSurvey.DataType, 'LFPMontage/LfpMontageTimeDomain') || ...
                  strcmp(dataSurvey.DataType, 'BrainSenseSurveys/BrainSenseSurveysTimeDomain')
            plotSurveys(dataSurvey, savepath, fsave, showFig)
            clear dataSurvey
        else
            error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
        end

    % Identifier
    elseif exist('dataIdentifier', 'var')
        if strcmp(dataIdentifier.DataType, 'BrainSenseSurveys/BrainSenseSurveysTimeDomain')
            plotSurveys(dataIdentifier, savepath, fsave, showFig)
            clear dataIdentifier
        else
            error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
        end
       
    % Setup
    elseif exist('dataSetup', 'var')
        if strcmp(dataSetup.DataType, 'MostRecentInSessionSignalCheck/SenseChannelTests')
            plotSetup(dataSetup, savepath, fsave, showFig)
            clear dataSetup
        else
            error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
        end

    % Timeline
    elseif exist('dataTimeline', 'var')
        if strcmp(dataTimeline.DataType, 'LFPTrendLogs')
            plotTimeline(dataTimeline, savepath, fsave, showFig)
            clear dataTimeline
        else
            error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
        end

    % Events
    elseif exist('dataEvents', 'var')
        if strcmp(dataEvents.DataType, 'LfpFrequencySnapshotEvents')
            plotEvents(dataEvents, savepath, fsave, showFig)
            clear dataEvents
        else
            error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
        end

    else
        error('The file %s is of incorrect format. Only outputs from the ALFA toolbox are accepted.', [files(i).folder filesep fsave])
    end
end

