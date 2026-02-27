% Author: B.J. Keulen
% Date: 13-01-2026
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
%
% Function for the collection of JSON files and their directories, and to
% define and create the path to which output structured will be saved.
% 
% INPUT
%   dataset         =   str defining type of dataset. Either single file 
%                       (0), a single folder (1) or a set of folders (2)
%
% OUTPUT
%   fileData        =   struct containing the path and name of the folder
%                       containing the selected JSON files, the file names
%                       and number of folders and files, all if applicable
%   savepath        =   str defining path to which output structures should
%                       be saved

function [fileData, savepath] = getFileData(dataset)

    % Define output variable
    fileData.rootPath = '';
    fileData.rootName = '';
    fileData.folders = '';
    fileData.files = '';
    fileData.nfolders = nan;
    fileData.nfiles = nan;

    % Get single .json file
    if dataset == 0
        [fileData.files, fileData.rootPath] = uigetfile('*.json', 'Select a single JSON file');
    
        if ischar(fileData.files) && ischar(fileData.rootPath)
            fileData.rootPath = fileData.rootPath(1:end-1);
            [~,fileData.rootName] = fileparts(fileData.rootPath);
            if strcmp(fileData.files(end-4:end),'.json')
                savepath = [fileData.rootPath filesep fileData.files(1:end-5) '_unpacked'];
                fileData.nfolders = 1;
                fileData.nfiles = 1;
            else
                error('No JSON file selected.')
            end
        else
            error('No file selected.')
        end
    
    % Get folder and all .json files within folder
    elseif dataset == 1
        fileData.rootPath = uigetdir('', 'Select a folder containing JSON files');
    
        if ischar(fileData.rootPath)
            fileData.files = dir([fileData.rootPath filesep '**' filesep '*.json']);
    
            % Check if files are found, remove duplicates if applicable
            if isempty(fileData.files)
                error('No JSON files found in selected folder.')
            else
                if length(fileData.files) > 1
                    fileTable = struct2table(fileData.files);
                    [~,idx,~] = unique(fileTable.name);
                    fileData.files = fileData.files(idx,:);
                end
                fileData.nfiles = length(fileData.files);
                fileData.nfolders = 1;
    
                [~,fileData.rootName] = fileparts(fileData.rootPath);
                savepath = [fileData.rootPath filesep 'folder_' fileData.rootName '_unpacked'];
            end
        else
            error('No folder selected.')
        end
    
    % Get folder and subfolders
    elseif dataset == 2
        fileData.rootPath = uigetdir('', 'Select a folder containing subfolders with JSON files');
    
        if ischar(fileData.rootPath)
            fileData.folders = dir(fileData.rootPath);
            fileData.folders = fileData.folders([fileData.folders.isdir],:);
            fileData.folders = fileData.folders(~matches({fileData.folders.name},[".",".."]),:);
            fileData.nfolders = length(fileData.folders);
    
            % Check if files are found, remove duplicates if applicable
            if isempty(fileData.folders)
                error('No subfolders found. Select a folder containing subfolders')
            else
                filecount = [];
                for c = 1:length(fileData.folders)
                    filecount = [filecount; dir([fileData.rootPath filesep fileData.folders(c).name filesep '**' filesep '*.json'])];
                end
    
                if isempty(filecount)
                    error('No JSON files found in subfolders of selected folderset.')
                else
                    [~,fileData.rootName] = fileparts(fileData.rootPath);
                    savepath = [fileData.rootPath filesep 'folderset_' fileData.rootName '_unpacked'];
                end
            end
        else
            error('No folder selected.')
        end
    else
        error('Incorrect input for dataset. Options are 0, 1 or 2.')
    end

    % Create folder to save data and plots into, if needed
    if ~exist(savepath, 'dir')
       mkdir(savepath)
    end
end