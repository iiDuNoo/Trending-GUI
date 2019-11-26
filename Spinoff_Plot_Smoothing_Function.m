% This function approximates data sets, attempting to capture important patterns while eliminating noise in the trends 
% that distracts the user
%
% Inputs: 
%       data                A cell array [mx1] that contains data for each
%                           individual event. Each cell should contain an
%                           [nx1] array.
%
%       event_dates         A cell array [mx1] containing corresponding
%                           dates for each event (cell) in "data"
%
%       Title               A string to be used as the title.
%
%       y_label_txt         A string to be used as the y axis label.
%
%
% Outputs:
%       None

function Spinoff_Plot_Smoothing_Function(data,event_dates,Title,y_label_txt,column_int)
%setting low
max=-10000000;
min=100000000;

sm_data= zeros(size(data));

%%Pop-Up Menu
dlg_title = 'Type of Moving Average';
sum=0;
answer = inputdlg('Enter the type of Moving Average (i.e. 2 = two-point M.A., 3 = three-point M.A.)',dlg_title);
answer=str2double(answer);
sorttype = inputdlg('Do you want to use max(1) or min(2) values?','Type of Data');
sorttype=str2double(sorttype);
if (sorttype==1)%max option
    for i=1:length(data)%Extracts cell
        for x=1:(length(data{i}(:,column_int)))%looks at all items in column 2 (often varies on number of points)
            if data{i}(x,column_int)>max%If the data point is larger than any other in the column, it overwrites
               max = data{i}(x,column_int);%Allocates largest points to Max
            end
        end
        sm_data(i) = max;%sm_data is filled with maximum values 
        max=-100000000;%Resets
        
    end 

elseif (sorttype ==2) %min option
    for i=1:length(data)%Extracts cell
        for x=1:(length(data{i}(:,column_int)))%looks at all items in column 2 (often varies on number of points)
            if data{i}(x,column_int)<min%If the data point is larger than any other in the column, it overwrites
               min = data{i}(x,column_int);%Allocates largest points to Max
            end
        end
        sm_data(i) = min;%sm_data is filled with maximum values 
        min=100000000;%Resets
    end

else %Invalid input
    disp('3')
    msgbox('Please select an option that is either 1 or 2')
    return
end
    
filtered_data=zeros(size(data)-answer);
if (answer>0)%for any case
    
    offset=0;
    for y=1:((length(sm_data)-((2*answer))))%Goes through data set
        for p=1:(2*answer)+1
           sum= sum+(sm_data(p+offset)); %adds x amount on each side plus middle 
        end

       filtered_data(y) = sum/((2*answer)+1);%y+count follows the index of the middle point
        sum=0;
        offset=offset+1;
        if (p+offset)>length(data)
            break; 
        end
    end     
    
elseif (answer>(0.5*length(sm_data)))
    msgbox('Please select a smaller Moving Average Value')
    return
else 
    msgbox('Please enter a valid option')
    return
end
%adjust size
for g=1:answer
    sm_data(g)=[];
    sm_data(length(sm_data)-g)=[];
    event_dates(g)=[];
    event_dates(length(sm_data)-g)=[];
end
%% Generate Smoothing Function plot
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
for i=1:length(filtered_data)
    %Makes array for aestetic purposes
     yvals(i) = filtered_data(i);  
end

%Plots max values
plot (1:length(filtered_data),yvals,1:length(filtered_data),sm_data)
legend ('smoothed data','original data')
if length(data)<30
    grid on
end
set(axid,'xtick',1:length(data)-answer)
set(axid,'xlim',[0.5,length(data)+0.5-answer])
set(axid,'xticklabel',event_dates)
end
