% Authors: M.J. Stam & B.J. Keulen
% Date: 29-03-2024
%
% Copyright 2026 B.J. Keulen, M.J. Stam & J.T. Boonstra
% SPDX-License-Identifier: Apache-2.0
% 
% Look for duplicates in Timeline or Events. The same data may be present
% in multiple json files, which is why this is needed. Timeline data points
% spaced more than 4 minutes apart will be removed. 
%
% INPUT
%   data            =   struct with either Timeline or Events data
%   fieldsTL        =   cell with names of fields in dataTimeline which
%                       contain extracted LFP data. Only needed if data
%                       is dataTimeline
%
% OUTPUT
%   data            =   struct with either Timeline or Events data

function data = checkDuplicates(data, fieldsTL)

    arguments
        data;
        fieldsTL = {};
    end

    % Define max difference between two recordings
    max_diff = minutes(4);

    if data.DataType == "LFPTrendLogs"
        if ~isempty(fieldsTL)

            % Go through data and check differences
            tdiff = diff(data.Data{:,'DateTime'});
            idx_diff = logical([0; tdiff >= max_diff]);
            data.Data = data.Data(idx_diff,:);

        end

    elseif data.DataType == "LfpFrequencySnapshotEvents"

        % Retrieve unique events
        events = unique(data.Data{:,'EventName'});

        % Go through full data and check differences
        for e = 1:length(events)
            eventData = data.Data(strcmp(data.Data{:,'EventName'}, events{e}),:);

            eventAll = eventData(~cellfun(@isempty, eventData.Frequency),:);
            diff_full = diff(eventAll.DateTime);
            eventAll([false; diff_full < max_diff],:) = [];

            eventEmpty = eventData(cellfun(@isempty, eventData.Frequency),:);
            for i = 1:height(eventEmpty)
                if ~any((eventAll.DateTime - eventEmpty{i,'DateTime'}) < max_diff)
                    eventAll = [eventAll; eventEmpty(i,:)];
                end
            end

            if e == 1
                allData = eventAll;
            else
                allData = [allData; eventAll];
            end
        end

        % Remove Events with less than 100 values
        idx = [];
        for e = 1:height(allData)
            if isempty(allData{e,'Frequency'}{:}) || length(allData{e,'Frequency'}{:}) == 100
                idx = [idx, e];
            end
        end
        allData = allData(idx,:);

        % Sort on time and replace
        data.Data = sortrows(allData,'DateTime');

    end
end