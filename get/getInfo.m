% Author: B.J. Keulen
% Date: 13-01-2026
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to collect general info on JSON file into one struct. This
% struct will later be added to all output structures.
% 
% INPUT
%   folder          =   str with path to folder
%   filename        =   str with filename
%   js              =   JSON file loaded into MATLAB
%
% OUTPUT
%   info            =   struct with general info on JSON file

function info = getInfo(folder, filename, js)

    info = {};
    info.OriginalFolder = folder;
    info.OriginalFile = filename;
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
end
