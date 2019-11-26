% This function plots the input trend data and its temperature compensated
% version using the input compensation. It artificially brings the trend
% values back to zero degrees C.
%
% Inputs: 
%       trend_data          An array [mx1] that contains the trend data
%
%       temperature_data    An array [mx1] that contains the temperature
%                           data
%
%       correction          A value corresponding to the correlation
%                           between the trend value and temperature (value
%                           is in units per degree C). This function will
%                           use that slope value to bring the trend values
%                           back to 0 degrees C.
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


function Spinoff_Plot_Temperature_Compensated_Manually(trend_data,temperature_data,event_dates,correction,correction_units,correction_temperature,Title,y_label_txt,legend_txt)

years = zeros(size(event_dates));
xvals = zeros(length(event_dates),2); % For year separators
corrected_vals = zeros(size(trend_data));
for i = 1:length(event_dates)
    corrected_vals(i) = trend_data(i) + correction * (correction_temperature - temperature_data(i));
    years(i) = str2double(event_dates{i}(1:4));
    if ~any(years(i)==xvals(:,2))
        xvals(i,:) = [i-0.5,years(i)];
    end
end
xvals = xvals(xvals(:,1)~=0,:);

%% Generate Temperature plot
fid = figure('position',[100,100,1100,800]);
axid = axes('parent',fid);
set(axid,'position',[0.06,0.20,0.92,0.72])
title(axid,Title)
if iscell(Title)
    set(fid,'name',['Spinoff: Temperature Compensated Trend ' Title{1}])
else
    set(fid,'name',['Spinoff: Temperature Compensated Trend ' Title])
end
ylabel(axid,y_label_txt)
set(axid,'ygrid','on')
hold(axid,'all')

plot(axid,1:length(trend_data),trend_data,'-ob');
plot(axid,1:length(corrected_vals),corrected_vals,'-*r');

legend(axid,legend_txt,['Temperature Corrected ' legend_txt sprintf(' (%f %s)',correction,correction_units)])

set(axid,'xtick',1:length(trend_data))
set(axid,'xlim',[0.5,length(trend_data)+0.5])
set(axid,'xticklabel',[])

ylims = get(axid,'ylim');
xticks = get(axid,'XTick');

ypos = 1.01*ylims(1) - 0.01*ylims(2);
text(xticks,ypos*ones(size(xticks)),event_dates,'Parent',axid,...
    'HorizontalAlignment','right','rotation',90,'FontSize',8);

for i = 1:size(xvals,1)
    plot(axid,[xvals(i,1) xvals(i,1)],ylims,':','color','black');
end
text((xvals(:,1)+[xvals(2:end,1);length(years)+0.5])/2,...
    (ylims(1)+0.05*(ylims(2)-ylims(1)))*ones(length(xvals),1),...
    num2str(xvals(:,2)),'rotation',90,'Parent',axid);

set(axid,'ylim',ylims)
end
