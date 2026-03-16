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
            iDiff = logical([0; tdiff >= max_diff]);
            data.Data = data.Data(iDiff,:);

        end

    elseif data.DataType == "LfpFrequencySnapshotEvents"

        % Retrieve unique events
        events = unique(data.Data{:,'EventName'});

        % Go through full data and check differences
        for e = 1:length(events)

            % Get all events of the same type, split based on data or not
            eventData = data.Data(strcmp(data.Data{:,'EventName'}, events{e}),:);
            eventFull = eventData(~cellfun(@isempty, eventData.Frequency),:);
            eventEmpty = eventData(cellfun(@isempty, eventData.Frequency),:);

            % Delete absolute duplicates of events with data based on datetime
            [~, iU, ~] = unique(eventFull(:,'DateTime'));
            eventFull = eventFull(iU,:);

            % Find indices of events with data less than 10 minutes apart and if
            % PSD magnitude is equal, consider as duplicate and remove the latter
            % NOTE: added 13-03-2026; after publication of user documentation
            iDiff = (find(diff(eventFull.DateTime) < minutes(10)))+1;
            for i = length(iDiff):-1:1
                if isempty(eventFull{iDiff(i),'Magnitude'}{1})
                    if isequal(eventFull{iDiff(i),'Magnitude'}{2},eventFull{iDiff(i)-1,'Magnitude'}{2})
                        eventFull(iDiff(i),:) = [];
                    end
                else
                    if isequal(eventFull{iDiff(i),'Magnitude'}{1},eventFull{iDiff(i)-1,'Magnitude'}{1})
                        eventFull(iDiff(i),:) = [];
                    end
                end
            end

            % Add empty events if more than 4 minutes apart from other events
            for i = 1:height(eventEmpty)
                if ~any((eventFull.DateTime - eventEmpty{i,'DateTime'}) < max_diff)
                    eventFull = [eventFull; eventEmpty(i,:)];
                end
            end

            % Create or add to table with all events
            if e == 1
                allData = eventFull;
            else
                allData = [allData; eventFull];
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