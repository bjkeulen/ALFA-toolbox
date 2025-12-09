% Author: B.J. Keulen
% Date: 14-05-2024
%
% Copyright 2025 B.J. Keulen and M.J. Stam
% SPDX-License-Identifier: Apache-2.0
%
% Impute NaN values for missing values in dataTimeline, defined as gaps in
% time of more than 20 minutes. NaN will then be imputed each 10 minutes 
% after the most recent time point with data until the next time point with
% data. Note that this may result in a time difference of less than 10 
% minutes between the last NaN of each missing segment and the next time 
% point with data.
%
% INPUT
%   dataTimeline    =   struct with Timeline data as table
%
% OUTPUT
%   dataTimeline    =   dataTimeline with NaN imputed for missing values


function dataTimeline = addNaN(dataTimeline)

    % Calculate timesteps and find gaps >15 min
    timesteps = diff(dataTimeline.Data{:,'DateTime'});
    istart = find(timesteps>minutes(15));
    
    % Check for lost data
    nantime = [];
    if ~isempty(istart)
        for i = 1:length(istart)
            nantime = [nantime dataTimeline.Data{istart(i),'DateTime'} + minutes(10) : minutes(10) : dataTimeline.Data{istart(i)+1,'DateTime'} - minutes(10)];
        end
    end

    % Create table with NaN or empty values
    nantable = table(nantime(:),'VariableNames',{'DateTime'});
    for f = 2:length(dataTimeline.Data.Properties.VariableNames)
        if strcmp(dataTimeline.Data.Properties.VariableNames{f},'ActiveGroup')
            nantable{:,'ActiveGroup'} = strings(length(nantime),1);
        elseif any(strcmp(dataTimeline.Data.Properties.VariableNames{f},'SensingChannel'))
            nantable{:,dataTimeline.Data.Properties.VariableNames{f}} = cell(length(nantime),2);
        else
            nantable{:,dataTimeline.Data.Properties.VariableNames{f}} = nan(length(nantime),2);
        end
    end

    % Join table and sort on time
    if ~isempty(nantable)
        dataTimeline.Data = sortrows([dataTimeline.Data; nantable],1);
    end

end
