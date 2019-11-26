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
%       Title               A string to be used as the title.
%
%       y_label_txt         A string to be used as the y axis label.
%
% Outputs:
%       None


function Spinoff_Plot_Temperature_by_year(trend_data,temperature_data,event_dates,Title,y_label_txt)

years = zeros(size(event_dates));
for i = 1:length(event_dates)
    years(i) = str2double(event_dates{i}(1:4));
end

year_splits = [0,find(diff(years)),length(years)];

n_years = length(year_splits) - 1;

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

h_lines = zeros(2,n_years);
legend_txt = cell(2,n_years);
for i = 1:n_years
    xvals = temperature_data(year_splits(i)+1:year_splits(i+1));
    yvals = trend_data(year_splits(i)+1:year_splits(i+1));
    
    if length(unique(yvals)) > 1
        line_poly = polyfit(xvals,yvals,1);
    else
        line_poly = [0,yvals(1)];
    end
    
    h_lines(1,i) = plot(axid,xvals,yvals,'*');
    colour = get(h_lines(1,i),'color');
    
    h_lines(2,i) = plot(axid,xvals,polyval(line_poly,xvals),'color',colour);
    legend_txt(:,i) = {sprintf('%i data',years(year_splits(i+1)));...
        sprintf('%i linear fit: %.5f x %+.2f',years(year_splits(i+1)),line_poly)};
end

legend(h_lines(:),legend_txt(:),'location','eastoutside')

end
