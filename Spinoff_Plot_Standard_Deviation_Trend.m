% This function plots the standard deviation of the input data. It can also
% plot a baseline of the standard deviation values using the average values
% and a 3-sigma bound (the standard deviation of the standard deviations).
%
% Inputs: 
%       data                A cell array [mx1] that contains data for each
%                           individual event. Each cell should contain an
%                           [nx1] array, which is the data segment to take
%                           the standard deviation of.
%
%       event_dates         A cell array [mx1] containing corresponding
%                           dates for each event (cell) in "data"
%
%       Title               A string to be used as the title.
%
%       y_label_txt         A string to be used as the y axis label.
%
%       bool_baseline       A boolean representing whether or not to
%                           include baselines in the standard deviation
%                           plots.
%
% Outputs:
%       None

function Spinoff_Plot_Standard_Deviation_Trend(data,event_dates,Title,y_label_txt,bool_baseline)

std_dev = zeros(size(data));

for i = 1:length(data)
   %Applies std to each variable point in array
    std_dev(i) = std(data{i});
end
outliers = true;
inliers = true(1,length(std_dev));
while outliers
    avg_val = mean(std_dev(inliers));
    sig_val = 3*std(std_dev(inliers));
    new_inliers = std_dev < avg_val + sig_val & std_dev > avg_val - sig_val;
    %If new inliers == old inliers, it becomes an outlier
    if all(new_inliers == inliers)
        outliers = false;
    else
        inliers = new_inliers;
    end
end

%% Generate Std Dev Function plot
fid = figure('position',[100,100,1100,800]);
axid = axes('parent',fid);
set(axid,'position',[0.06,0.20,0.92,0.72])
title(axid,Title)
if iscell(Title)
    set(fid,'name',['Spinoff: Smoothing Function ' Title{1}])
else
    set(fid,'name',['Spinoff: Smoothing Function ' Title])
end

ylabel(axid,y_label_txt)
hold(axid,'all')

hlines(1) = plot(axid,find(inliers),std_dev(inliers),'-*b');
if bool_baseline
    hlines(2) = plot(axid,[1,length(std_dev)],[avg_val,avg_val],'color',[0,0.5,0]);
    hlines(3) = plot(axid,[1,length(std_dev)],[avg_val+sig_val,avg_val+sig_val],'-m');
    hlines(4) = plot(axid,[1,length(std_dev)],[avg_val-sig_val,avg_val-sig_val],'-c');
    if any(~inliers)
        hlines(5) = plot(axid,find(~inliers),std_dev(~inliers),'*r');
    end
    plot(axid,1,0,'w') % Ensures 0 is included in the y axis
end

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

if length(hlines) == 5 && bool_baseline
    legend(hlines,'Average values',sprintf('Overall Average: %.3f',avg_val),...
        sprintf('Upper 3%s bound: %.3f','\sigma',avg_val+sig_val),...
        sprintf('Lower 3%s bound: %.3f','\sigma',avg_val-sig_val),...
        'Outliers','location','best')
elseif bool_baseline
    legend(hlines,'Average values',sprintf('Overall Average: %.3f',avg_val),...
        sprintf('Upper 3%s bound: %.3f','\sigma',avg_val+sig_val),...
        sprintf('Lower 3%s bound: %.3f','\sigma',avg_val-sig_val),...
        'location','best')
end

end
