% This function plots the input trend data vs its temperature. It creates a
% trend line for each year.
%
% Inputs: 
%       trend_data          An array [mx1] that contains the trend data
%
%       temperature_data    An array [mx1] that contains the temperature
%                           data
%
%       event_dates         A cell array [mx1] containing corresponding
%                           dates for each event (data point). Dates should
%                           already be sorted by year.
%
%       events_per_group    The number of events each group contains. This
%                           will start from the most recent and work
%                           backwards to the early events. If there are any
%                           left that are not enough to create their own
%                           group, they will be added to the first group.
%
%       Title               A string to be used as the title.
%
%       y_label_txt         A string to be used as the y axis label.
%
% Outputs:
%       None


function Spinoff_Plot_Temperature_by_event(trend_data,temperature_data,event_dates,events_per_group,Title,y_label_txt)

years = zeros(size(event_dates));
for i = 1:length(event_dates)
    years(i) = str2double(event_dates{i}(1:4));
end

group_ranges = length(trend_data):-events_per_group:1;
if group_ranges(end) < events_per_group && length(group_ranges) > 1
    group_ranges = group_ranges(end-1:-1:1);
else
    group_ranges = group_ranges(end:-1:1);
end
group_ranges = [1,group_ranges(1:end-1)+1;group_ranges];

n_groups = size(group_ranges,2);

%% Generate Temperature plot
fid = figure('position',[100,100,1100,800]);
axid = axes('parent',fid);
title(axid,Title)
if iscell(Title)
    set(fid,'name',['Spinoff: Temperature Scatter plot ' Title{1}])
else
    set(fid,'name',['Spinoff: Temperature Scatter plot ' Title])
end
xlabel(axid,'Temperature (C)')
ylabel(axid,y_label_txt)
grid(axid,'on')
hold(axid,'all')

h_lines = zeros(2,n_groups);
legend_txt = cell(2,n_groups);
for i = 1:n_groups
    xvals = temperature_data(group_ranges(1,i):group_ranges(2,i));
    yvals = trend_data(group_ranges(1,i):group_ranges(2,i));
    
    if length(unique(yvals)) > 1
        line_poly = polyfit(xvals,yvals,1);
    else
        line_poly = [0,yvals(1)];
    end
    
    h_lines(1,i) = plot(axid,xvals,yvals,'*');
    colour = get(h_lines(1,i),'color');
    
    h_lines(2,i) = plot(axid,xvals,polyval(line_poly,xvals),'color',colour);
    legend_txt(:,i) = {sprintf('Group %i data (%i events %i-%i)',i,length(xvals),years(group_ranges(1,i)),years(group_ranges(2,i)));...
        sprintf('Group %i linear fit: %.5f x %+.2f',i,line_poly)};
end

legend(h_lines(:),legend_txt(:),'location','eastoutside')

end
