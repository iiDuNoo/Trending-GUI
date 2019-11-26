% This function plots the Power Spectral Density (PSD) of all the motor
% rates which are in the variable data. It will extract the longest
% segment containing 10 Hz data. It also applies a hamming window to the
% data.
%
% Inputs: 
%       data                A cell array [mx1] that contains data for each
%                           individual event. Each cell should contain an
%                           [nx3] array, where the first column is time
%                           data, the second column is position data, and
%                           the third column is motor rate data.
%
%       event_dates         A cell array [mx1] containing corresponding
%                           dates for each event (cell) in "data"
%
%       Title               A string to be used as the first line of the
%                           title. The second and third line are
%                           automatically assigned.
%
% Outputs:
%       None

function Spinoff_Plot_FFT(data,event_dates,Title)
%% Extract the longest segment containing 10 Hz data

percent_included = 80; % Specifies the minimum percent of events to use in FFT analysis

data_lengths = zeros(size(data));
time_col = 1;
pos_col = 2;
rate_col = 3;

for i = 1:length(data)
    offsets = rem(data{i}(:,time_col),0.1)*2*pi;
    xs = cos(offsets);
    ys = sin(offsets);
    offset = atan(sum(ys)/sum(xs))/2/pi; % Average of a cyclic quantity
    data{i}(:,time_col) = data{i}(:,time_col) - offset;
    rounded_time = round(data{i}(:,time_col)*10); % converted to deciseconds
    good_data_bool = true(length(data{i}(:,time_col)),1);

    for j = min(rounded_time):max(rounded_time) % 10 hz data should increment by integer deciseconds (j)
        potential_data_bool = rounded_time == j;
        if sum(potential_data_bool) > 1 % Multiple data found at that time
            % Select the one closest to the expected time
            [~,min_i] = min(abs(data{i}(:,time_col)-j/10));
            if ~potential_data_bool(min_i)
                error('Error: Bug in code')
            else
                good_data_bool(potential_data_bool) = false;
                good_data_bool(min_i) = true;
            end
        end
    end
    data{i} = data{i}(good_data_bool,:); % Remove all points that do not fit in 10Hz data
    
    rounded_time = round(data{i}(:,time_col)*10);
    bad_jumps_i = ...
        [0,find(diff(rounded_time(:)') ~= 1),length(rounded_time)];
    [big_10hz,big_i] = max(diff(bad_jumps_i));
    start_i = bad_jumps_i(big_i) + 1;
    end_i = bad_jumps_i(big_i + 1);
    data{i} = data{i}(start_i:end_i,:);
    rounded_time = round(data{i}(:,time_col)*10);
    
    if any(diff(rounded_time(:)') ~= 1)
        error('Error: Converting to 10 Hz failed')
    end
    if length(bad_jumps_i) > 2
        warning(['Warning: File contains gap in 10 Hz data. Largest '...
            'continuous segment will be taken for FFT: %i/%i data '...
            'points (%s)'],big_10hz,length(good_data_bool),event_dates{i})
    end
    data_lengths(i) = size(data{i},1);
end

sorted_lengths = sort(data_lengths);
pct_pts_i = floor(length(data_lengths)*(1-percent_included/100));
if ~pct_pts_i
    pct_pts_i = 1;
end
pct_pts = sorted_lengths(pct_pts_i);
min_pts = min(data_lengths);
max_pts = max(data_lengths);
if min_pts>max_pts/2
    num_pts = min_pts;
else
    num_pts = pct_pts;
end
data = data(data_lengths>=num_pts);
fprintf(['\nData length range: %i to %i. The first %i data points will '...
    'be taken for all files (%.1f%% of the files will be used)\n'],...
    min_pts,max_pts,num_pts,length(data)/length(data_lengths)*100)

%% Generate motor rate vs position plot (with position offset)
fid = figure('position',[100,100,1100,800]);
axid = axes('parent',fid);
title(axid,{Title; '3D Motor Rate vs Position - Segments used for FFT'})
xlabel(axid,'Position (Offset to start at zero)')
ylabel(axid,'Event Date')
zlabel(axid,'Absolute Motor Rate')
hold(axid,'all')

for i = 1:length(data)
    data{i} = data{i}(1:num_pts,:);
    data{i}(:,pos_col) = data{i}(:,pos_col) - data{i}(1,pos_col); % Apply position offset.
    data{i}(:,rate_col) = abs(data{i}(:,rate_col)); % Absolute rates only.
    
    xvals = data{i}(:,pos_col)';
    yvals = i*ones(size(xvals));
    zvals = data{i}(:,rate_col)';
    surface([xvals;xvals],[yvals;yvals],[zvals;zvals],[zvals;zvals],...
        'facecol','no','edgecol','interp','linew',2,'parent',axid)
end


set(axid,'ytick',1:length(data))
set(axid,'ylim',[0.5,length(data)+0.5])
set(axid,'yticklabel',event_dates)
view(axid,[-25,40])

%% Perform FFT Analysis and generate PSD plot
fid = figure('position',[100,100,1100,800]);
axid = axes('parent',fid);
% Extract motor rates and determine their PSDs using fft
for i = 1:length(data)
    motor_rate_data = data{i}(:,rate_col);
    avg_rate = mean(motor_rate_data);
    fft_signal = ...
        (motor_rate_data-avg_rate).*hamming(size(motor_rate_data,1));
    fft_data = fft(fft_signal,num_pts);
    freq_conv=2*pi/(abs(avg_rate));
    
    frequencies = (0:(num_pts/2)-1)/(num_pts*0.1)*freq_conv;
        
    power_spec_sig = (fft_data.*conj(fft_data)/num_pts)/freq_conv;
    PSD = power_spec_sig(1:length(frequencies));
    xvals = frequencies;
    yvals = i*ones(size(frequencies));
    zvals = PSD';
    surface([xvals;xvals],[yvals;yvals],[zvals;zvals],[zvals;zvals],...
        'facecol','no','edgecol','interp','linew',2,'parent',axid)
end

view(axid,[-25 40])
%X-axis label on graph
xlabel(axid,'Frequency (cycles/mtr rev)')
%Y-axis label on graph
ylabel(axid,'Event')
%Z-axis label on graph
zlabel(axid,'Power (amp^2/freq)')
title(axid,{Title;'Power Spectrum Density vs Event';...
    ['Number of Points Used for FFT: ' num2str(num_pts)]})
if iscell(Title)
    set(fid,'name',['Spinoff: FFT ' Title{1}])
else
    set(fid,'name',['Spinoff: FFT ' Title])
end

set(axid, 'Ytick', 1:length(data))
set(axid, 'Yticklabel',event_dates,'fontsize',8)
set(axid,'Ylim',[0.5,length(data)+0.5])

end
