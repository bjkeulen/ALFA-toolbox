% Author: B.J. Keulen
% Date: 13-01-2026
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Function to initialise the output structures for Timeline and Events data
% 
% INPUT
%   datatype        =   str with type of data. Either 'Timeline' or 'Events'
%
% OUTPUT
%   structure       =   empty output structure for either Timeline or
%                       Events data

function structure = getStructures(datatype, fields)

    structure = struct;
    if strcmp(datatype,'Timeline')
        structure.DataType = 'LFPTrendLogs';
    elseif strcmp(datatype,'Events')
        structure.DataType = 'LfpFrequencySnapshotEvents';
    else
        error("Incorrect input for datatype. Options are 'Timeline' or 'Events'")
    end
    structure.Info = [];
    structure.Labels = {'LEFT','RIGHT'};
    structure.Data = struct;
    if strcmp(datatype,'Timeline')
        for n = 1:length(fields)
            structure.Data.(fields{n}) = [];
        end
        structure.History = struct;
    else
        for n = 1:length(fields)
            structure.Data.(fields{n}) = [];
        end
    end
end