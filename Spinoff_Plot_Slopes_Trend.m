% This function plots the slope of the input data.
%
% Inputs: 
%       data                A cell array [mx1] that contains data for each
%                           individual event. Each cell should contain an
%                           [nx2] array, where the first column is the x
%                           axis data and the second column is the y axis
%                           data.
%
%       event_dates         A cell array [mx1] containing corresponding
%                           dates for each event (cell) in "data"
%
%       Title               A string to be used as the title.
%
%       y_label_txt         A string to be used as the y axis label.
%
% Outputs:
%       None

function Spinoff_Plot_Slopes_Trend(data,event_dates,Title,y_label_txt)

slopes = zeros(size(data));
for i = 1:length(data)
    xvals = data{i}(:,1);
    yvals = data{i}(:,2);
    linear_fit_poly = polyfit(xvals,yvals,1);
    slopes(i) = linear_fit_poly(1);
end

%% Generate Slopes plot
fid = figure('position',[100,100,1100,800]);
axid = axes('parent',fid);
set(axid,'position',[0.06,0.20,0.92,0.72])
title(axid,Title)
if iscell(Title)
    set(fid,'name',['Spinoff: Slopes ' Title{1}])
else
    set(fid,'name',['Spinoff: Slopes ' Title])
end

ylabel(axid,y_label_txt)
hold(axid,'all')

plot(axid,1:length(slopes),slopes,'-*b');

set(axid,'xtick',1:length(data))
set(axid,'xlim',[0.5,length(data)+0.5])
set(axid,'xticklabel',[])

ylims = get(axid,'ylim');
xticks = get(axid,'XTick');

ypos = 1.01*ylims(1) - 0.01*ylims(2);
text(xticks,ypos*ones(size(xticks)),event_dates,'Parent',axid,...
    'HorizontalAlignment','right','rotation',90,'FontSize',8);

years = zeros(size(event_dates));
xvals = zeros(length(event_dates),2);
for i = 1:length(event_dates)
    years(i) = str2double(event_dates{i}(1:4));
    if ~any(years(i)==xvals(:,2))
        xvals(i,:) = [i-0.5,years(i)];
    end
end
xvals = xvals(xvals(:,1)~=0,:);
for i = 1:size(xvals,1)
    plot(axid,[xvals(i,1) xvals(i,1)],ylims,':','color','black');
end
text((xvals(:,1)+[xvals(2:end,1);length(years)+0.5])/2,...
    (ylims(1)+0.05*(ylims(2)-ylims(1)))*ones(length(xvals),1),...
    num2str(xvals(:,2)),'rotation',90,'Parent',axid);
set(axid,'ygrid','on')

end
