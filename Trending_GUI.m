function varargout = Trending_GUI(varargin)
%% TRENDING_GUI MATLAB code for Trending_GUI.fig
%      TRENDING_GUI, by itself, creates a new TRENDING_GUI or raises the existing
%      singleton*.
%
%      H = TRENDING_GUI returns the handle to a new TRENDING_GUI or the handle to
%      the existing singleton*.
%
%      TRENDING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRENDING_GUI.M with the given input arguments.
%
%      TRENDING_GUI('Property','Value',...) creates a new TRENDING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Trending_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Trending_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Trending_GUI

% Last Modified by GUIDE v2.7 20-Aug-2015 15:01:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Trending_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Trending_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Trending_GUI is made visible.
function Trending_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
%% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Trending_GUI (see VARARGIN)

% Choose default command line output for Trending_GUI
handles.output = hObject;

% Set initial values for the GUI
set(hObject,'Name',['Trend-based, Holistic Observation, Manipulation '...
    '& Analysis Software (Trending GUI) V2.7'])

set(handles.System_selection,'String',{...
    'SSRMS LEE',...
    'SPDM OTCM',...
    'SPDM LEE',...
    })
handles.loading.processing_functions = {...
    @process_SSRMS_LEE_data,...
    @process_OTCM_data,...
    @process_SPDM_LEE_data,...
    };
%Menu for Spinoffs
set(handles.Spinoff_selection,'String',{...
    'FFT',...
    'Standard Deviation Trend',...
    'Trend of the Slopes',...
    'Difference of 2 Trend Lines',...
    'Temperature vs Event',...
    'Trend vs Temperature Scatter Plot',...
    'Temperature Compensated Trend Line',...
    'Smoothing Function'})
    

% Update handles structure
handles.loading.folders = struct('files',{},'name',{});
handles.loading.import_dir = ...
    '\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI';
handles.loading.session_dir = ['\\vsgroups\Groups\SE\'...
    'Trending\Raw_data_and_Trending_GUI\Saved_GUI_Sessions'];

handles.data = [];

handles.figdata = [];

% Initialize and set default values
handles.settings.raw_data_visible = true;
handles.settings.enhanced_view = false;
handles.settings.grids_visible = true;
handles.settings.baselines_visible = false;
handles.settings.animations_enabled = false;
handles.settings.current_view = [-25,40];
handles.settings.reversed_y_axis = false;
handles.settings.colour_scheme = 'Z gradient';
handles.settings.marker_style = 'none';
handles.settings.marker_size = 6;
handles.settings.line_style = '-';
handles.settings.line_width = 0.5;

handles = initialize_auto_sessions(handles);

guidata(hObject, handles);

% UIWAIT makes Trending_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Trending_GUI_OutputFcn(hObject, eventdata, handles) 
%% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Load_automatic_button.
%Loads added session (see Add_auto_session_menu_Callback
function Load_auto_session_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Load_automatic_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.Automatic_Sessions_list,'value');
session_names = get(handles.Automatic_Sessions_list,'string');
session = ['\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI'...
    '\Automatic_Updating_Scripts\Automatic_Sessions\'...
    session_names{val} '.mat'];

Load_session_menu_Callback([], session, handles)



% --- Executes on button press in Add_folder_button.
%Imports data from text files in the selected folders
function Add_folder_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Add_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder = uigetdir([handles.loading.import_dir,'\..\'],...
    'Select Folder Containing Data in Text Files');

if any(strcmp(folder,{handles.loading.folders.name}))
    msgbox('Folder has already been selected')
    return
end

if folder
    handles.loading.import_dir = folder;
    folder_i = length(handles.loading.folders) + 1;
    handles.loading.folders(folder_i).name = folder;
end

handles = update_import_status(handles);

guidata(handles.output,handles)



% --- Executes on button press in Remove_folder_button.
%Removes added folders from the query
function Remove_folder_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Remove_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vals = get(handles.Folder_list,'value');
handles.loading.folders(vals) = [];

handles = update_import_status(handles);

guidata(handles.output,handles)



% --- Executes on button press in System_selection.
%Can select which system the data is from (SPDM LEE/UTCM or SSRMS LEE)
function System_selection_dropdown_Callback(hObject, eventdata, handles)
%% hObject    handle to System_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = update_import_status(handles);

guidata(handles.output,handles)


% --- Executes on button press in Manual_import_button.
%Processes imported folders
function Import_process_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Manual_import_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bool_new_process_fcn = false;
fcn_i = get(handles.System_selection,'value');
if ~isfield(handles.loading,'active_processing_fcn')
    handles.loading.active_processing_fcn = ...
        handles.loading.processing_functions{fcn_i}([],1);
    bool_new_process_fcn = true;
end

if ~isequal(eventdata,'auto')
    handles.loading.processing_fcn_handle = ...
        handles.loading.processing_functions{fcn_i};
end
processing_fcn = handles.loading.processing_fcn_handle;

% Make sure processing function is the same. If not, ask if user wants to
% import and process all available data
if ~isequal(processing_fcn([],1),handles.loading.active_processing_fcn)
    if ~isequal(eventdata,'auto') % Supress user interaction for auto updating
        response = questdlg(['Warning: Either the system selected or '...
            'the processing function version has changed. All raw data '...
            'must be imported to continue'],['Different Processing '...
            'Function Found'],'Continue','Cancel','Continue');
        if ~strcmp(response,'Continue')
            return
        end
    end
    handles.data = [];
    handles.loading.active_processing_fcn = processing_fcn([],1);
    bool_new_process_fcn = true;
end

num_files = 0;
num_new = 0;
for i = 1:length(handles.loading.folders)
    num_new = num_new + ...
        sum((~[handles.loading.folders(i).files.isprocessed]&...
        ~[handles.loading.folders(i).files.failed_processing]));
    num_files = num_files + length(handles.loading.folders(i).files);
    if bool_new_process_fcn
        [handles.loading.folders(i).files.isprocessed] = deal(false);
        [handles.loading.folders(i).files.failed_processing] = deal(false);
    end
end
% Import and processs all data that has not already been processed.
num_processed = 0;
all_files = cell(num_files,1);
files_i = 0;
set(handles.Import_status,'foregroundcolor','blue')
for i = 1:length(handles.loading.folders)
    for j = 1:length(handles.loading.folders(i).files)
        files_i = files_i + 1;
        all_files(files_i) = {[handles.loading.folders(i).name,'\',...
            handles.loading.folders(i).files(j).name]};
        if ~handles.loading.folders(i).files(j).isimported && ...
                ~handles.loading.folders(i).files(j).isprocessed && ...
                ~handles.loading.folders(i).files(j).failed_processing
            try
                handles.loading.folders(i).files(j).rawdata = ...
                    importdata([handles.loading.folders(i).name,'\',...
                    handles.loading.folders(i).files(j).name]);
            catch err
                msgbox(sprintf('An error occured while importing:\n%s',...
                    [handles.loading.folders(i).name,'\',...
                    handles.loading.folders(i).files(j).name]))
                disable_filtering(handles);
                rethrow(err)
            end
        end
        if ~handles.loading.folders(i).files(j).isprocessed && ...
                ~handles.loading.folders(i).files(j).failed_processing
            try
                num_processed = num_processed + 1;
                set(handles.Import_status,'string',sprintf(['Processing'...
                    ' %i/%i (%.1f%% complete)'],num_processed,num_new,...
                    100*num_processed/num_new))
                drawnow
                fprintf('\nProcessing %i/%i: %s...\n...%s\n',...
                    num_processed,num_new,[handles.loading.folders(...
                    i).name,'\'],handles.loading.folders(i).files(j).name)
                data = processing_fcn(...
                    handles.loading.folders(i).files(j).rawdata,0);
                data.eventsources = {[handles.loading.folders(i).name,...
                    '\',handles.loading.folders(i).files(j).name]};
            catch err
                msgbox(sprintf('An error occured while processing:\n%s',...
                    [handles.loading.folders(i).name,'\',...
                    handles.loading.folders(i).files(j).name]))
                disable_filtering(handles);
                rethrow(err)
            end
            if isempty(data.parameterdata)
                % Handle the no good data case
                handles.loading.folders(i).files(j).failed_processing = ...
                    true;
            elseif isempty(handles.data)
                handles.data = data;
            else
                % Append data
                new_i = length(handles.data.eventdates)+1;
                handles.data.filterdata(new_i) = data.filterdata;
                handles.data.parameterdata(new_i) = data.parameterdata;
                handles.data.eventdates(new_i) = data.eventdates;
                handles.data.temperature(new_i) = data.temperature;
                handles.data.eventsources(new_i) = data.eventsources;
            end
        end
    end
end
% Sort data
if ~isempty(handles.data)
    conversion_vals = [10^12,10^9,10^7,10^5,10^3,1];
    datevals = zeros(size(handles.data.eventdates));
    for event_i = 1:length(datevals)
        if any(strcmp(handles.data.eventsources{event_i},all_files))
            datevals(event_i) = sum(sscanf(...
                handles.data.eventdates{event_i},'%d%*c',6)' .* ...
                conversion_vals);
        else
            datevals(event_i) = 0; % Events to be removed
        end
    end
    [sorted_vals,sort_I] = sort(datevals);
    sort_I = sort_I(~~(sorted_vals));
    handles.data.filterdata = handles.data.filterdata(sort_I);
    handles.data.parameterdata = handles.data.parameterdata(sort_I);
    handles.data.eventdates = handles.data.eventdates(sort_I);
    handles.data.temperature = handles.data.temperature(sort_I);
    handles.data.eventsources = handles.data.eventsources(sort_I);
end

handles = update_import_status(handles);
if bool_new_process_fcn && ~isempty(handles.data)
    handles = reset_filters(handles);
end
handles = disable_3d_plotting(handles);
guidata(handles.output,handles)


% --- Executes on selection change in Data_limit_selection.
%Modifies the options in the Select Filter Value drop down menu
function Filter_selection_Callback(hObject, eventdata, handles)
%% hObject    handle to Data_limit_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Data_limit_selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Data_limit_selection
filter_i = get(handles.Filter_selection,'value');
set(handles.Filter_value_Selection,'value',1)
set(handles.Filter_value_Selection,'String',...
    handles.data.filters(filter_i).options)


% --- Executes on button press in Apply_filters_button.
%Modifies the selected custom filter values
function Add_filter_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Apply_filters_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filter_i = get(handles.Filter_selection,'value');
update_val = get(handles.Filter_value_Selection,'value');
selected_val = handles.data.selected_filter_vals(filter_i);
active_val = handles.data.active_filter_vals(filter_i);

if update_val == selected_val
    msgbox('Filter has already been selected')
else
    handles.data.selected_filter_vals(filter_i) = update_val;
    cur_status = get(handles.Filter_status,'string');
    msg = [handles.data.filters(filter_i).name,': ',...
        handles.data.filters(filter_i).options{...
        handles.data.active_filter_vals(filter_i)}];
    if update_val ~= active_val
        msg = [msg,' (changed to ',...
            handles.data.filters(filter_i).options{update_val},')'];
        cur_status{filter_i + 1} = msg;
        set(handles.Filter_status,'string',cur_status)
        set(handles.Filter_status,'ForegroundColor','red')
    elseif all(handles.data.selected_filter_vals == ...
            handles.data.active_filter_vals)
        cur_status{filter_i + 1} = msg;
        set(handles.Filter_status,'string',cur_status)
        set(handles.Filter_status,'ForegroundColor',[0,0.5,0])
    end
end

guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Apply_filters_button.
%Applies the selected custom filters, data restictions, and data limits, in
%that orders
function Apply_filters_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Apply_filters_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.active_filter_vals = handles.data.selected_filter_vals;
handles.data.active_restrictions = handles.data.selected_restrictions;
handles.data.active_data_limits = handles.data.selected_data_limits;

handles.data.filtered_data = handles.data.parameterdata;

handles.data.full_params = {handles.data.parameters.string};
for i = 1:length(handles.data.full_params)
    dependance = handles.data.parameters(i).unitdep;
    if dependance
        unit_i = handles.data.filters(dependance).values{handles.data.active_filter_vals(dependance)};
        unit = handles.data.parameters(i).unit{unit_i};
    else
        unit = handles.data.parameters(i).unit;
    end
    handles.data.full_params{i} = ...
        [handles.data.full_params{i}, ' (' unit ')'];
end
set(handles.Restriction_data_selection,'String',handles.data.full_params);
set(handles.Data_limit_selection,'String',handles.data.full_params);

for i = 1:length(handles.data.parameterdata)
    bool_good_data = true(size(handles.data.filtered_data{i},1),1);
    
    % Apply custom filters
    for j = 1:length(handles.data.filters)
        bool_good_data = bool_good_data & ...
            ismember(handles.data.filterdata{i}(:,j),...
            handles.data.filters(j).values{handles.data.active_filter_vals(j)});
    end
    handles.data.filtered_data{i} = handles.data.filtered_data{i}(bool_good_data,:);
    
    bool_good_data = true(size(handles.data.filtered_data{i},1),1);
    % Apply event restrictions
    for j = 1:size(handles.data.active_restrictions,1)
        bool_contains = true(size(handles.data.filtered_data{i},1),1);
        if ~isempty(handles.data.active_restrictions{j,3}) % Min restriction
            param_i = handles.data.active_restrictions{j,2};
            bool_contains = bool_contains & ...
                (handles.data.filtered_data{i}(:,param_i) > handles.data.active_restrictions{j,3});
        end
        if ~isempty(handles.data.active_restrictions{j,4}) % Max restriction
            param_i = handles.data.active_restrictions{j,2};
            bool_contains = bool_contains & ...
                (handles.data.filtered_data{i}(:,param_i) < handles.data.active_restrictions{j,4});
        end
        bool_good_data = bool_good_data & ...
            (any(bool_contains) == handles.data.active_restrictions{j,1}); % Good depends on must contain vs must not contain
    end
    handles.data.filtered_data{i} = handles.data.filtered_data{i}(bool_good_data,:);
    
    % Apply custom data filter
    for j = 1:size(handles.data.active_data_limits,1)
        bool_contains = true(size(handles.data.filtered_data{i},1),1);
        if ~isempty(handles.data.active_data_limits{j,2}) % Min restriction
            param_i = handles.data.active_data_limits{j,1};
            bool_contains = bool_contains & ...
                (handles.data.filtered_data{i}(:,param_i) > handles.data.active_data_limits{j,2});
        end
        if ~isempty(handles.data.active_data_limits{j,3}) % Max restriction
            param_i = handles.data.active_data_limits{j,1};
            bool_contains = bool_contains & ...
                (handles.data.filtered_data{i}(:,param_i) < handles.data.active_data_limits{j,3});
        end
        handles.data.filtered_data{i}(~bool_contains,:) = NaN; % Don't display this data
    end
    if any(all(isnan(handles.data.filtered_data{i}))) || (size(handles.data.filtered_data{i},1) < 2)
        handles.data.filtered_data{i} = [];
    end
end

status_str = cell(length(handles.data.filters),1);
for i = 1:length(handles.data.filters)
    status_str{i} = [handles.data.filters(i).name,': ',...
        handles.data.filters(i).options{handles.data.active_filter_vals(i)}];
end

if isfield(handles.data,'filtered_dates')
    prev_dates = handles.data.filtered_dates;
else
    prev_dates = {};
end

no_data = cellfun(@isempty,handles.data.filtered_data);
handles.data.filtered_data = handles.data.filtered_data(~no_data);
handles.data.filtered_dates = handles.data.eventdates(~no_data);
handles.data.filtered_temperatures = handles.data.temperature(~no_data);
handles.data.omitted_dates = handles.data.eventdates(no_data);

if all(no_data) % No data found
    set(handles.Filter_status,'string',[{'No data found for:'};status_str])
    set(handles.Filter_status,'ForegroundColor','red')
    handles = disable_3d_plotting(handles);
else
    event_stat = {sprintf('Data found for %i of %i events:',...
        sum(~no_data),length(handles.data.eventdates))};
    set(handles.Filter_status,'string',[event_stat;status_str;])
    set(handles.Filter_status,'ForegroundColor',[0,0.5,0])
    
    if ~isequal(prev_dates,handles.data.filtered_dates)
        event_visible_strs = strcat(cellstr(num2str((1:length(handles.data.filtered_dates))')),{': '},handles.data.filtered_dates');
        set(handles.Event_table,'Data',[{'All Events',true};...
            event_visible_strs,...
            num2cell(true(size(handles.data.filtered_dates')))])
    end
    
    handles = enable_3D_plotting(handles);
    
    handles = update_figure_handles(handles);
    handles = update_plot_data(handles);
    if ~isempty(handles.figdata)
        handles = enable_full_plotting(handles);
    end
end
handles = update_restrictions_list(handles);
handles = update_data_limit_list(handles);

guidata(handles.output,handles)



% --- Executes on button press in Export_button.
%Exports an excel summary of the selected trend lines
function Export_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = update_figure_handles(handles);

[xlsxfile,xlsxdir] = uiputfile('*.xlsx','Select name of Excel File being saved','Saved_trends.xlsx');

h = msgbox('Exporting summary. Please wait');
for i = find([handles.data.trends.is_visible])
    excel_sheet = {'Trend type:',handles.data.trends(i).type};
    excel_sheet(3:4,1) = {'Filters:';'Filter selections:'};
    excel_sheet(6,1) = {'Restrictions:'};
    excel_sheet(12,1) = {'Limits:'};
    
    % Filter info
    excel_sheet(3,2:length(handles.data.filters)+1) = {handles.data.filters.name};
    filter_vals = cell(1,length(handles.data.filters));
    for j = 1:length(handles.data.filters)
        filter_vals{j} = handles.data.filters(j).options{handles.data.trends(i).filters(j)};
    end
    excel_sheet(4,2:length(handles.data.filters)+1) = filter_vals;
    
    % Restriction info
    restrictions = cell(4,size(handles.data.trends(i).restrictions,1));
    for j = 1:size(restrictions,2)
        must_vals = {'Must not contain','Must contain'};
        restrictions(1,j) = must_vals(handles.data.trends(i).restrictions{j,1}+1);
        restrictions(2,j) = {[handles.data.parameters(handles.data.trends(i).restrictions{j,2}).string,' between:']};
        min_val = handles.data.trends(i).restrictions{j,3};
        if isempty(min_val)
            min_val = 'inf';
        else
            min_val = num2str(min_val);
        end
        max_val = handles.data.trends(i).restrictions{j,4};
        if isempty(max_val)
            max_val = 'inf';
        else
            max_val = num2str(max_val);
        end
        restrictions(3:4,j) = {min_val;max_val};
    end
    excel_sheet(7:10,1:size(restrictions,2)) = restrictions;
    
    % Limit info
    limits = cell(3,size(handles.data.trends(i).limits,1));
    for j = 1:size(limits,2)
        limits(1,j) = {[handles.data.parameters(handles.data.trends(i).limits{j,1}).string,':']};
        min_val = handles.data.trends(i).limits{j,2};
        if isempty(min_val)
            min_val = 'inf';
        else
            min_val = num2str(min_val);
        end
        max_val = handles.data.trends(i).limits{j,3};
        if isempty(max_val)
            max_val = 'inf';
        else
            max_val = num2str(max_val);
        end
        limits(2:3,j) = {min_val;max_val};
    end
    excel_sheet(13:15,1:size(limits,2)) = limits;
    
    % Trend line(s) values
    trend_info = [{'Event date:'};handles.data.trends(i).dates'];
    for j = 1:length(handles.figdata)
        col = j*3-1;
        param_x = handles.figdata(j).params(1);
        param_y = handles.figdata(j).params(2);

        dependance = handles.data.parameters(param_x).unitdep;
        if dependance
            unit_i = handles.data.filters(dependance).values{handles.data.active_filter_vals(param_x)};
            unit_x = handles.data.parameters(param_x).unit{unit_i};
        else
            unit_x = handles.data.parameters(param_x).unit;
        end
        dependance = handles.data.parameters(param_y).unitdep;
        if dependance
            unit_i = handles.data.filters(dependance).values{handles.data.active_filter_vals(param_y)};
            unit_y = handles.data.parameters(param_y).unit{unit_i};
        else
            unit_y = handles.data.parameters(param_y).unit;
        end
        
        header1 = [handles.data.trends(i).type(1:end-1),' ',handles.data.parameters(param_y).string,' (',unit_y,')'];
        if strcmp(handles.data.trends(i).type,'Averages')
            header2 = [handles.data.trends(i).type(1:end-1),' ',handles.data.parameters(param_x).string,' (',unit_x,')'];
        else
            header2 = [handles.data.parameters(param_x).string,' at ',...
                handles.data.trends(i).type(1:end-1),' ',handles.data.parameters(param_y).string,' (',unit_x,')'];
        end
        trend_info(1,col:col+1) = {header1,header2};
        trend_info(2:end,col:col+1) = num2cell(handles.data.trends(i).vals(param_y,:,[param_y,param_x]));
    end
    
    % Temperature info
    temperature_names = unique([handles.data.trends(i).temperatures.name],'stable');
    temp_names = reshape([handles.data.trends(i).temperatures.name],...
        length(handles.data.trends(i).temperatures(1).name),...
        length(handles.data.trends(i).temperatures))';
    
    temperature_data = reshape([handles.data.trends(i).temperatures.value],...
        length(handles.data.trends(i).temperatures(1).name),length(handles.data.trends(i).temperatures))';
    trend_info{1,size(trend_info,2)+1} = [];
    for j = 1:length(temperature_names)
        temp_bool = strcmp(temp_names,temperature_names{j});
        temps = NaN(size(temp_bool,1),1);
        for k = 1:size(temp_bool,1)
            temp_i = find(temp_bool(k,:),1,'first');
            if ~isempty(temp_i)
                temps(k) = temperature_data(k,temp_i);
            end
        end
        trend_info(:,size(trend_info,2)+1) = [temperature_names(j);num2cell(temps)];
    end
    excel_sheet(17:16+size(trend_info,1),1:size(trend_info,2)) = trend_info;
    
    % Write to excel file
    try
        xlswrite([xlsxdir,xlsxfile],excel_sheet,['Trend line ' num2str(i)])
    catch %err
        msgbox('Exporting summary failed. Make sure that the excel file is not already open')
%         rethrow(err)
    end
end
close(h)



% --- Executes on button press in Add_restriction_button.
%Adds the restriction to the selected restrictions list
function Add_restriction_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Add_restriction_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_must = 2 - get(handles.Contain_selection,'value'); % 0 == must not have, 1 == must have
new_param = get(handles.Restriction_data_selection,'value');
new_min = sscanf(get(handles.Min_restriction,'String'),'%f',1);
new_max = sscanf(get(handles.Max_restriction,'String'),'%f',1);
new_restriction = {new_must,new_param,new_min,new_max};

for i = 1:size(handles.data.selected_restrictions,1)
    if isequal(handles.data.selected_restrictions(i,:),new_restriction)
        msgbox('Restriction has already been selected')
        return
    end
end

i = size(handles.data.selected_restrictions,1)+1;

handles.data.selected_restrictions(i,:) = new_restriction;

handles = update_restrictions_list(handles);

guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Remove_restriction_button.
%Removes the restriction from the selected restriction list
function Remove_restriction_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Remove_restriction_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.data.selected_restrictions)
    msgbox('No restrictions to remove')
    return
end

i_to_remove = get(handles.Restriction_list,'value');

i_to_remove = ismember(1:size(handles.data.selected_restrictions,1), ...
    i_to_remove(i_to_remove <= size(handles.data.selected_restrictions,1)));

if ~any(i_to_remove)
    msgbox('Select a valid restriction to remove')
end

handles.data.selected_restrictions = ...
    handles.data.selected_restrictions(~i_to_remove,:);

handles = update_restrictions_list(handles);

guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Add_data_limit_button.
%Adds the data limit to the selected limits list
function Add_data_limit_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Add_data_limit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

new_param = get(handles.Data_limit_selection,'value');
new_min = sscanf(get(handles.Min_limit,'String'),'%f',1);
new_max = sscanf(get(handles.Max_limit,'String'),'%f',1);
new_data_filter = {new_param,new_min,new_max};

if any([handles.data.selected_data_limits{:,1}] == new_param)
    i = find([handles.data.selected_data_limits{:,1}] == new_param);
else
    i = size(handles.data.selected_data_limits,1)+1;
end

handles.data.selected_data_limits(i,:) = new_data_filter;

handles = update_data_limit_list(handles);

guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Remove_data_limit_button.
%Removes the data limit from the selected limits list
function Remove_data_limit_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Remove_data_limit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.data.selected_data_limits)
    msgbox('No data limits to remove')
    return
end

i_to_remove = get(handles.Limit_list,'value');

i_to_remove = ismember(1:size(handles.data.selected_data_limits,1), ...
    i_to_remove(i_to_remove <= size(handles.data.selected_data_limits,1)));

if ~any(i_to_remove)
    msgbox('Select a valid data limit to remove')
end

handles.data.selected_data_limits = ...
    handles.data.selected_data_limits(~i_to_remove,:);

handles = update_data_limit_list(handles);

guidata(get(hObject,'Parent'),handles)



% --- Executes on value change in Event_table.
function Event_table_Callback(hObject, eventdata, handles)
%% hObject    handle to Event_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

if eventdata.Indices(1) == 1
    prev_data = get(handles.Event_table,'Data');
    [prev_data{:,2}] = deal(prev_data{1,2});
    set(handles.Event_table,'Data',prev_data)
end
handles = update_plot_settings(handles);



% --- Executes on selection change in Event_table.
function Event_table_selection_Callback(hObject, eventdata, handles)
%% hObject    handle to Event_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'UserData',eventdata.Indices(:,1))



% --- Executes on value change in Deselect_Group_context.
%Disables the visibility of selected elements
function Deselect_Group_context_Callback(hObject, eventdata, handles)
%% hObject    handle to Deselect_Group_context (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

indices = get(handles.Event_table,'UserData');
if ~isempty(indices)
    prev_data = get(handles.Event_table,'Data');
    [prev_data{indices,2}] = deal(false);
    set(handles.Event_table,'Data',prev_data)
    handles = update_plot_settings(handles);
end



% --- Executes on value change in Select_Group_context.
%Enables the visibility of selected elements
function Select_Group_context_Callback(hObject, eventdata, handles)
%% hObject    handle to Select_Group_context (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

indices = get(handles.Event_table,'UserData');
if ~isempty(indices)
    prev_data = get(handles.Event_table,'Data');
    [prev_data{indices,2}] = deal(true);
    set(handles.Event_table,'Data',prev_data)
    handles = update_plot_settings(handles);
end

% --- Executes on button press in Generate_3d_button.
%Generates a new 3D plot 
function Generate_3d_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Generate_3d_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

selected_params = [get(handles.X_axis_selection,'value'),...
    get(handles.Z_axis_selection,'value')];
for i = 1:length(handles.figdata)
    if all(selected_params == handles.figdata(i).params)
        msgbox('Figure already exists')
        return
    end
end

% Create 3d figure
new_i = length(handles.figdata)+1;
handles.figdata(new_i).figID = figure('position',[55+15*new_i,145-15*new_i,1100,800]);
ax1 = axes('Position',[0 0 1 1],'Visible','off','parent',handles.figdata(end).figID);
set(ax1,'DeleteFcn',sprintf('set(get(%.13f,''children''),''DeleteFcn'','''')',ax1))
set(ax1,'UserData',sprintf(['text(0.85,0.01,char([71,101,110,101,114,'...
    '97,116,101,100,32,98,121,32,84,72,79,77,65,83]),''parent'',%.13f,'...
    '''fontsize'',6,''SelectionHighlight'',''off'',''DeleteFcn'',get(%.13f,''UserData''))'],ax1,ax1))
axes('Position',[0 0 1 1],'Visible','off','parent',handles.figdata(end).figID);
text(0.85,0.01,char([71,101,110,101,114,97,116,101,100,32,98,121,32,84,...
    72,79,77,65,83]),'parent',ax1,'fontsize',6,'SelectionHighlight',...
    'off','DeleteFcn',get(ax1,'UserData'))
handles.figdata(new_i).axesID = axes('parent',handles.figdata(end).figID);
handles.figdata(new_i).params = selected_params;
handles.figdata(new_i).h3dlines = [];
handles.figdata(new_i).h2dlines = [];
handles.figdata(new_i).enhanced_h_lines = [];
handles.figdata(new_i).enhanced_h_years = [];
handles.figdata(new_i).enhanced_h_y_txt = [];
handles.figdata(new_i).htrends = [];
handles.figdata(new_i).hbaselines = [];

view(handles.figdata(new_i).axesID,handles.settings.current_view)

handles = update_plot_data(handles);
handles = enable_full_plotting(handles);

guidata(handles.output,handles)



% --- Executes on button press in Generate_trend_button.
%Creates a new trend line element
function Generate_trend_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Generate_trend_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

trend_i = length(handles.data.trends) + 1;

handles.data.trends(trend_i).is_visible = true;

if get(handles.Max_radio,'value')
    handles.data.trends(trend_i).type = 'Maximums';
elseif get(handles.Avg_radio,'value')
    handles.data.trends(trend_i).type = 'Averages';
elseif get(handles.Min_radio,'value')
    handles.data.trends(trend_i).type = 'Minimums';
else
    error('Something went wrong')
end

on_off_vals = get(handles.Event_table,'data');
on_off_vals = [on_off_vals{2:end,2}];

handles.data.trends(trend_i).dates = handles.data.filtered_dates(on_off_vals);
handles.data.trends(trend_i).temperatures = handles.data.filtered_temperatures(on_off_vals);
handles.data.trends(trend_i).filters = handles.data.active_filter_vals;
handles.data.trends(trend_i).restrictions = handles.data.active_restrictions;
handles.data.trends(trend_i).limits = handles.data.active_data_limits;

handles.data.trends(trend_i).vals = ...
    zeros(length(handles.data.parameters),sum(on_off_vals),...
    length(handles.data.parameters));

for i = 1:size(handles.data.trends(trend_i).vals,1) % Parameter that finds the peak
    for j = 1:size(handles.data.trends(trend_i).vals,2) % Event index
        data_i = find(strcmp(handles.data.filtered_dates,handles.data.trends(trend_i).dates(j)));
        switch handles.data.trends(trend_i).type
            case 'Maximums'
                [~,max_i] = max(handles.data.filtered_data{data_i}(:,i));
                handles.data.trends(trend_i).vals(i,j,:) = handles.data.filtered_data{data_i}(max_i,:);
            case 'Averages'
                is_num = all(~isnan(handles.data.filtered_data{data_i}),2);
                handles.data.trends(trend_i).vals(i,j,:) = mean(handles.data.filtered_data{data_i}(is_num,:),1);
            case 'Minimums'
                [~,min_i] = min(handles.data.filtered_data{data_i}(:,i));
                handles.data.trends(trend_i).vals(i,j,:) = handles.data.filtered_data{data_i}(min_i,:);
            otherwise
                error('Something different went wrong')
        end
    end
end

for i = 1:length(handles.data.trends)-1
    if isequal(handles.data.trends(i).vals,handles.data.trends(trend_i).vals)
        msgbox('Trend already exists')
        return % Doesn't save changes to data since guidata() is not called
    end
end

handles = update_plot_data(handles);

guidata(handles.output,handles)


% --- Executes on button press in Toggle_trend_button.
%Toggles the visibility of seleced trend lines
function Show_Hide_trend_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Toggle_trend_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

hide_val = get(handles.Trend_list,'value');
if ~isempty(handles.data.trends)
    for i = hide_val
        handles.data.trends(i).is_visible = ~handles.data.trends(i).is_visible;
    end
end

handles = update_plot_data(handles);

guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Clear_trend_button.
%Deletes selected trend line element 
function Clear_trend_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Clear_trend_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

del_val = get(handles.Trend_list,'value');
handles.data.trends(del_val) = [];

set(handles.Trend_list,'value',[])

handles = update_plot_data(handles);

guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Baseline_check.
%Toggles the visibility of the plots baseline
function Baseline_check_Callback(hObject, eventdata, handles)
%% hObject    handle to Baseline_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Baseline_check
handles.settings.baselines_visible = get(handles.Baseline_check,'value');
handles = update_figure_handles(handles);
handles = update_plot_settings(handles);
guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Raw_data_check.
%Toggles viewing the raw data in the plot
function Raw_data_check_Callback(hObject, eventdata, handles)
%% hObject    handle to Raw_data_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Raw_data_check
handles.settings.raw_data_visible = get(handles.Raw_data_check,'value');
handles = update_figure_handles(handles);
handles = update_plot_settings(handles);
guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Baseline_check.
%Offers a broader view of the data by extending the Y-Axis
function Enhancedviewcheck_Callback(hObject, eventdata, handles)
%% hObject    handle to Baseline_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Baseline_check
handles.settings.enhanced_view = get(handles.Enhancedviewcheck,'value');
handles = update_figure_handles(handles);
handles = update_plot_settings(handles);
guidata(get(hObject,'Parent'),handles)


% --- Executes on button press in Baseline_check.
%Toggles the gridlines on the plot
function Grid_check_Callback(hObject, eventdata, handles)
%% hObject    handle to Baseline_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Baseline_check
handles.settings.grids_visible = get(handles.Grid_check,'value');
handles = update_figure_handles(handles);
handles = update_plot_settings(handles);
guidata(get(hObject,'Parent'),handles)

% --- Executes on button press in Spinoff_button.
%Generates selected Spinoff in dropdown menu
function Spinoff_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Spinoff_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

spinoff_val = get(handles.Spinoff_selection,'value');

on_off_vals = get(handles.Event_table,'data');
on_off_vals = [on_off_vals{2:end,2}];

spin_dates = handles.data.filtered_dates(on_off_vals);

title_1 = handles.data.processed_fcn;
lim_str = 'Limits: ';
for i = 1:size(handles.data.active_data_limits,1)
    min_txt = num2str(handles.data.active_data_limits{i,2});
    if isempty(min_txt)
        min_txt = 'inf';
    end
    max_txt = num2str(handles.data.active_data_limits{i,3});
    if isempty(max_txt)
        max_txt = 'inf';
    end
    lim_str = [lim_str,handles.data.parameters(...
        handles.data.active_data_limits{i,1}).string,...
        ':(',min_txt,',',max_txt,'), '];
end
if isempty(handles.data.active_data_limits)
    lim_str = [lim_str,'None'];
end

for j = 1:length(handles.data.filters)
    title_1 = [title_1(1:end),...
        handles.data.filters(j).titletxt{handles.data.active_filter_vals(j)}];
end

switch spinoff_val
    case 1 % FFT Spinoff selected
        time_col = find(strcmpi({handles.data.parameters.string},'Time'),1);
        pos_col = find(strcmpi({handles.data.parameters.string},'Position'),1);
        rate_col = find(strcmpi({handles.data.parameters.string},'Motor Rate'),1);
        
        if length([time_col,pos_col,rate_col]) ~= 3
            msgbox('Error: At least one of the following was not found: Time, Position, Motor Rate')
            return
        end
        
        FFT_data = handles.data.filtered_data(on_off_vals);
        for i = 1:length(FFT_data)
            is_num = all(~isnan(FFT_data{i}),2);
            FFT_data{i} = FFT_data{i}(is_num,[time_col,pos_col,rate_col]);
        end
        
        FFT_title = title_1;
        
        Spinoff_Plot_FFT(FFT_data,spin_dates,FFT_title)
        
    case 2 % Standard Deviation Spinoff selected
        if isempty(handles.figdata)
            msgbox({'Please generate a 3d plot before generating a standard deviation plot',...
                'The z-axis value will be used to generate the spin-off plot'})
        end
        for i = 1:length(handles.figdata)
            param_i = handles.figdata(i).params(2);
            
            std_title = {title_1;lim_str;...
                ['Standard Deviation of ' ...
                handles.data.parameters(param_i).string]};
            bool_baseline = get(handles.Baseline_check,'value');
            
            std_data = handles.data.filtered_data(on_off_vals);
            for j = 1:length(std_data)
                is_num = all(~isnan(std_data{j}),2);
                std_data{j} = std_data{j}(is_num,param_i);
            end
            y_label_txt = handles.data.full_params{param_i};
            
            Spinoff_Plot_Standard_Deviation_Trend(...
                std_data,spin_dates,std_title,y_label_txt,bool_baseline)
        end
        
    case 3 % Slopes Spinoff selected
        if isempty(handles.figdata)
            msgbox({'Please generate a 3d plot before generating a slopes plot',...
                'The axes selected will be used to generate the spin-off plot'})
        end
        for i = 1:length(handles.figdata)
            param_x = handles.figdata(i).params(1);
            param_y = handles.figdata(i).params(2);
            
            dependance = handles.data.parameters(param_x).unitdep;
            if dependance
                unit_i = handles.data.filters(dependance).values{handles.data.active_filter_vals(param_x)};
                unit_x = handles.data.parameters(param_x).unit{unit_i};
            else
                unit_x = handles.data.parameters(param_x).unit;
            end
            dependance = handles.data.parameters(param_y).unitdep;
            if dependance
                unit_i = handles.data.filters(dependance).values{handles.data.active_filter_vals(param_y)};
                unit_y = handles.data.parameters(param_y).unit{unit_i};
            else
                unit_y = handles.data.parameters(param_y).unit;
            end
            
            y_label_txt = ['Slope (' unit_y '/' unit_x ')'];
            
            slope_title = {title_1;lim_str;...
                ['Slopes of ' ...
                handles.data.parameters(param_x).string ' vs '...
                handles.data.parameters(param_y).string]};
            
            slope_data = handles.data.filtered_data(on_off_vals);
            for j = 1:length(slope_data)
                is_num = all(~isnan(slope_data{j}),2);
                slope_data{j} = slope_data{j}(is_num,[param_x,param_y]);
            end
            
            Spinoff_Plot_Slopes_Trend(...
                slope_data,spin_dates,slope_title,y_label_txt)
        end
        
    case 4 % Difference of 2 trendlines
        handles = update_figure_handles(handles);
        
        bool_good_conditions = true;
        if isempty(handles.data.trends) || isempty(handles.figdata)
            bool_good_conditions = false;
        elseif sum([handles.data.trends.is_visible]) < 2
            bool_good_conditions = false;
        else
            bool_potential_trends = [handles.data.trends.is_visible];
            for i = 1:length(bool_potential_trends)
                bool_matching_dates = false;
                for j = 1:length(bool_potential_trends)
                    if i ~= j && isequal(handles.data.trends(i).dates,handles.data.trends(j).dates)
                        bool_matching_dates = true;
                    end
                end
                bool_potential_trends(i) = bool_potential_trends(i) && bool_matching_dates;
            end
            bool_good_conditions = bool_good_conditions && any(bool_potential_trends);
        end
        
        if ~bool_good_conditions
            msgbox('Error: There must be at least 2 visible trends with matching dates to generate this plot')
            return
        end
        
        possible_trends_str = get(handles.Trend_list,'string');
        possible_trends_i = find(bool_potential_trends);
        minuend_str_cell = [{'Select Minuend (Minuend - Subtrahend = Difference)'};...
            possible_trends_str(bool_potential_trends)];
        minuend_str = sprintf('%s\n',minuend_str_cell{1});
        for i = 2:length(minuend_str_cell)
            minuend_str = sprintf('%s\n',[minuend_str,num2str(i-1),'. ', minuend_str_cell{i}]);
        end
        minuend_response = inputdlg(minuend_str,'Select Minuend');
        if ~cellfun(@isempty,minuend_response)
        else
            msgbox('Spinoff generation terminated by user')
            return
        end
        minuend_i = sscanf(minuend_response{1},'%i');
        if isempty(minuend_i)
            msgbox('No Minuend selected')
            return
        elseif ~any(minuend_i == 1:sum(bool_potential_trends))
            msgbox('Selection out of range')
            return
        end
        minuend_i = possible_trends_i(minuend_i); % Referencing actual trend now.
        
        bool_possible_subtrahend = bool_potential_trends & ~(1:length(bool_potential_trends) == minuend_i);
        for i = find(bool_possible_subtrahend)
            if ~isequal(handles.data.trends(i).dates,handles.data.trends(minuend_i).dates)
                bool_possible_subtrahend(i) = false;
            end
        end
        possible_subtrahend_i = find(bool_possible_subtrahend);
        if length(possible_subtrahend_i) == 1
            subtrahend_i = possible_subtrahend_i;
            uiwait(msgbox(sprintf('%s selected as subtrahend',possible_trends_str{subtrahend_i})))
        else
            subtrahend_str_cell = [{'Select Subtrahend (Minuend - Subtrahend = Difference)'};...
                possible_trends_str(bool_possible_subtrahend)];
            subtrahend_str = sprintf('%s\n',subtrahend_str_cell{1});
            for i = 2:length(subtrahend_str_cell)
                subtrahend_str = sprintf('%s\n',[subtrahend_str,num2str(i-1),'. ', subtrahend_str_cell{i}]);
            end
            subtrahend_response = inputdlg(subtrahend_str,'Select Subtrahend');
            if ~cellfun(@isempty,subtrahend_response)
            else
                msgbox('Spinoff generation terminated by user')
                return
            end
            subtrahend_i = sscanf(subtrahend_response{1},'%i');
            if isempty(subtrahend_i)
                msgbox('No Subtrahend selected')
                return
            elseif ~any(subtrahend_i == 1:sum(bool_possible_subtrahend))
                msgbox('Selection out of range')
                return
            end
            subtrahend_i = possible_subtrahend_i(subtrahend_i); % Referencing actual trend now.
        end
        
        for i = 1:length(handles.figdata)
            param_i = handles.figdata(i).params(2);
            
            minuend_data = handles.data.trends(minuend_i).vals(param_i,:,param_i);
            subtrahend_data = handles.data.trends(subtrahend_i).vals(param_i,:,param_i);
            
            difference_dates = handles.data.trends(minuend_i).dates;
            difference_title = {title_1;...
                [possible_trends_str{minuend_i} ' Minus '...
                possible_trends_str{subtrahend_i}]};
            y_label_txt = handles.data.full_params{param_i};
            
            Spinoff_Plot_Difference_of_2_Trends(minuend_data,subtrahend_data,difference_dates,difference_title,y_label_txt)
        end
    case 5 % Temperature vs Event
        temp_str = sprintf('%s\n','Select temperature to use');
        unique_temps = unique([handles.data.filtered_temperatures.name],'stable');
        for i = 1:length(unique_temps)
            temp_str = sprintf('%s\n',[temp_str,num2str(i),'. ',unique_temps{i}]);
        end
        temp_response = inputdlg(temp_str,'Select Temperature');
        if ~cellfun(@isempty,temp_response)
        else
            msgbox('Spinoff generation terminated by user')
            return
        end
        temp_i = sscanf(temp_response{1},'%i');
        if isempty(temp_i)
            msgbox('No Temperature selected')
            return
        elseif temp_i > length(unique_temps)
            msgbox('Selection out of range')
            return
        end
        
        selected_temp_name = unique_temps{temp_i};
        
        temp_names = reshape([handles.data.filtered_temperatures(:).name],...
            length(handles.data.filtered_temperatures(1).name),...
            length(handles.data.filtered_temperatures))';
        
        good_data_bool = strcmp(selected_temp_name,temp_names);
        good_data_bool = good_data_bool(on_off_vals,:);
        
        temperature_data = reshape([handles.data.filtered_temperatures.value],...
            length(handles.data.filtered_temperatures(1).name),length(handles.data.filtered_temperatures))';
        temperature_data = temperature_data(on_off_vals,:);
        temperature_data = temperature_data(good_data_bool);
        
        contains_temp = any(good_data_bool,2);
        
        trend_dates = handles.data.filtered_dates(on_off_vals);
        trend_dates = trend_dates(contains_temp);
        
        temperature_title = {title_1;...
            [selected_temp_name ' vs event']};
        y_label_txt = 'Temperature (C)';
        
        Spinoff_Plot_Input_Trend_Basic(temperature_data,trend_dates,temperature_title,y_label_txt)
        
    case 6 % Trend vs temperature scatter plot
        handles = update_figure_handles(handles);
        
        bool_good_conditions = true;
        if isempty(handles.data.trends) || isempty(handles.figdata)
            bool_good_conditions = false;
        elseif sum([handles.data.trends.is_visible]) < 1
            bool_good_conditions = false;
        end
        
        if ~bool_good_conditions
            msgbox('Error: There must be at least 1 visible trend to generate this plot')
            return
        end
        
        possible_trends_str = get(handles.Trend_list,'string');
        possible_trends_bool = [handles.data.trends.is_visible];
        possible_trends_i = find(possible_trends_bool);
        trend_str_cell = [{'Select Trend on which to perform temperature analysis'};...
            possible_trends_str(possible_trends_bool)];
        trend_str = sprintf('%s\n',trend_str_cell{1});
        for i = 2:length(trend_str_cell)
            trend_str = sprintf('%s\n',[trend_str,num2str(i-1),'. ', trend_str_cell{i}]);
        end
        if length(possible_trends_i) == 1
            trend_i = possible_trends_i;
            uiwait(msgbox(sprintf('Trend selected:\n%s',trend_str_cell{2})))
        else
            trend_response = inputdlg(trend_str,'Select Trend');
            if ~cellfun(@isempty,trend_response)
            else
                msgbox('Spinoff generation terminated by user')
                return
            end
            trend_i = sscanf(trend_response{1},'%i');
            if isempty(trend_i)
                msgbox('No Trend selected')
                return
            elseif ~any(trend_i == 1:sum(possible_trends_bool))
                msgbox('Selection out of range')
                return
            end
            trend_i = possible_trends_i(trend_i); % Referencing actual trend now.
        end
        
        temp_str = sprintf('%s\n','Select temperature to use');
        unique_temps = unique([handles.data.trends(trend_i).temperatures.name],'stable');
        for i = 1:length(unique_temps)
            temp_str = sprintf('%s\n',[temp_str,num2str(i),'. ',unique_temps{i}]);
        end
        temp_response = inputdlg(temp_str,'Select Temperature');
        if ~cellfun(@isempty,temp_response)
        else
            msgbox('Spinoff generation terminated by user')
            return
        end
        temp_i = sscanf(temp_response{1},'%i');
        if isempty(temp_i)
            msgbox('No Temperature selected')
            return
        elseif temp_i > length(unique_temps)
            msgbox('Selection out of range')
            return
        end
        
        selected_temp_name = unique_temps{temp_i};
        
        temp_names = reshape([handles.data.trends(trend_i).temperatures(:).name],...
            length(handles.data.trends(trend_i).temperatures(1).name),...
            length(handles.data.trends(trend_i).temperatures))';
        
        good_data_bool = strcmp(selected_temp_name,temp_names);
        
        temperature_data = reshape([handles.data.trends(trend_i).temperatures.value],...
            length(handles.data.temperature(1).name),length(handles.data.trends(trend_i).temperatures))';
        temperature_data = temperature_data(good_data_bool);
        
        contains_temp = any(good_data_bool,2);
        
        trend_dates = handles.data.trends(trend_i).dates(contains_temp);
        
        good_data_bool = ~isnan(temperature_data);
        
        temperature_data = temperature_data(good_data_bool);
        trend_dates = trend_dates(good_data_bool);
        
        type_response = questdlg('Select type of temperature slopes to use:','Temperature Type','By year','By n events','By year');
        
        if strcmp(type_response,'By year')
            type = 1;
        elseif strcmp(type_response,'By n events')
            type = 2;
            n_events_response = inputdlg('Select number of events per line of best fit','Temperature analysis by n events');
            n_events = sscanf(n_events_response{1},'%i');
            if isempty(n_events)
                msgbox('Spinoff generation terminated by user')
                return
            end
        else
            msgbox('Temperature analysis terminated by user')
            return
        end
        
        for i = 1:length(handles.figdata)
            param_i = handles.figdata(i).params(2);
            
            trend_data = handles.data.trends(trend_i).vals(param_i,:,param_i)';
            trend_data = trend_data(good_data_bool);
            
            temperature_title = {title_1;...
                [possible_trends_str{trend_i} ' vs ' selected_temp_name]};
            y_label_txt = handles.data.full_params{param_i};
            
            if type == 1
                Spinoff_Plot_Temperature_by_year(trend_data,temperature_data,trend_dates,temperature_title,y_label_txt)
            elseif type == 2
                Spinoff_Plot_Temperature_by_event(trend_data,temperature_data,trend_dates,n_events,temperature_title,y_label_txt)
            end
        end
    case 7 % Temperature Compensated Peaks
        handles = update_figure_handles(handles);
        
        bool_good_conditions = true;
        if isempty(handles.data.trends) || isempty(handles.figdata)
            bool_good_conditions = false;
        elseif sum([handles.data.trends.is_visible]) < 1
            bool_good_conditions = false;
        end
        
        if ~bool_good_conditions
            msgbox('Error: There must be at least 1 visible trend to generate this plot')
            return
        end
        
        possible_trends_str = get(handles.Trend_list,'string');
        possible_trends_bool = [handles.data.trends.is_visible];
        possible_trends_i = find(possible_trends_bool);
        trend_str_cell = [{'Select Trend on which to perform temperature analysis'};...
            possible_trends_str(possible_trends_bool)];
        trend_str = sprintf('%s\n',trend_str_cell{1});
        for i = 2:length(trend_str_cell)
            trend_str = sprintf('%s\n',[trend_str,num2str(i-1),'. ', trend_str_cell{i}]);
        end
        if length(possible_trends_i) == 1
            trend_i = possible_trends_i;
            uiwait(msgbox(sprintf('Trend selected:\n%s',trend_str_cell{2})))
        else
            trend_response = inputdlg(trend_str,'Select Trend');
            if ~cellfun(@isempty,trend_response)
            else
                msgbox('Spinoff generation terminated by user')
                return
            end
            trend_i = sscanf(trend_response{1},'%i');
            if isempty(trend_i)
                msgbox('No Trend selected')
                return
            elseif ~any(trend_i == 1:sum(possible_trends_bool))
                msgbox('Selection out of range')
                return
            end
            trend_i = possible_trends_i(trend_i); % Referencing actual trend now.
        end
        
        temp_str = sprintf('%s\n','Select temperature to use');
        unique_temps = unique([handles.data.trends(trend_i).temperatures.name],'stable');
        for i = 1:length(unique_temps)
            temp_str = sprintf('%s\n',[temp_str,num2str(i),'. ',unique_temps{i}]);
        end
        temp_response = inputdlg(temp_str,'Select Temperature');
        if ~cellfun(@isempty,temp_response)
        else
            msgbox('Spinoff generation terminated by user')
            return
        end
        temp_i = sscanf(temp_response{1},'%i');
        if isempty(temp_i)
            msgbox('No Temperature selected')
            return
        elseif temp_i > length(unique_temps)
            msgbox('Selection out of range')
            return
        end
        
        selected_temp_name = unique_temps{temp_i};
        
        temp_names = reshape([handles.data.trends(trend_i).temperatures(:).name],...
            length(handles.data.trends(trend_i).temperatures(1).name),...
            length(handles.data.trends(trend_i).temperatures))';
        
        good_data_bool = strcmp(selected_temp_name,temp_names);
        
        temperature_data = reshape([handles.data.trends(trend_i).temperatures.value],...
            length(handles.data.temperature(1).name),length(handles.data.trends(trend_i).temperatures))';
        temperature_data = temperature_data(good_data_bool);
        
        contains_temp = any(good_data_bool,2);
        
        trend_dates = handles.data.trends(trend_i).dates(contains_temp);
        
        for i = 1:length(handles.figdata)
            param_i = handles.figdata(i).params(2);
            
            dependance = handles.data.parameters(param_i).unitdep;
            if dependance
                unit_i = handles.data.filters(dependance).values{handles.data.active_filter_vals(param_i)};
                unit_z = handles.data.parameters(param_i).unit{unit_i};
            else
                unit_z = handles.data.parameters(param_i).unit;
            end
            
            correction_response = ...
                inputdlg(sprintf('Specify compensation to use (unit: %s / degree C) for:\n%s',...
                unit_z,possible_trends_str{trend_i}));
            
            if ~cellfun(@isempty,correction_response)
            else
                msgbox('Spinoff generation terminated by user')
                return
            end
            
            correction = sscanf(correction_response{1},'%f');
            if isempty(correction)
                msgbox('No Temperature correction specified')
                return
            end
            correction_unit = sprintf('%s / degree C',unit_z);
            
            correction_temp_response = ...
                inputdlg(sprintf('Specify Temperature to compensate to (in degrees C). Commonly done to 0 degrees C'));
            
            if ~cellfun(@isempty,correction_temp_response)
            else
                msgbox('Spinoff generation terminated by user')
                return
            end
            
            correction_temperature = sscanf(correction_temp_response{1},'%f');
            if isempty(correction_temperature)
                msgbox('No Temperature specified')
                return
            end
            
            trend_data = handles.data.trends(trend_i).vals(param_i,:,param_i)';
            trend_data = trend_data(any(good_data_bool,2));
            
            trend_type = handles.data.trends(trend_i).type(1:end-1);
            
            temperature_title = {title_1;...
                ['Temperature Compensated (to ' num2str(correction_temperature)...
                ' degrees C) ' trend_type ' ' handles.data.parameters(param_i).string 's vs Event']};
            y_label_txt = [trend_type ' ' handles.data.full_params{param_i}];
            
            legend_txt = [trend_type ' ' handles.data.parameters(param_i).string 's'];
            
            Spinoff_Plot_Temperature_Compensated_Manually(...
                trend_data,temperature_data,trend_dates,correction,correction_unit,...
                correction_temperature,temperature_title,y_label_txt,legend_txt)
        end
        
    case 8 % Smoothing Function
    handles = update_figure_handles(handles);
    %Update column Input
    time_col = find(strcmpi({handles.data.parameters.string},'Time'),1);
    pos_col = find(strcmpi({handles.data.parameters.string},'Position'),1);
    rate_col = find(strcmpi({handles.data.parameters.string},'Motor Rate'),1);
    cur_col = find(strcmpi({handles.data.parameters.string},'Current'),1);
    abs_cur_col = find(strcmpi({handles.data.parameters.string},'Absolute Current'),1);
    abs_rate_col = find(strcmpi({handles.data.parameters.string},'Absolute Motor Rate'),1);
    load_param_col = find(strcmpi({handles.data.parameters.string},'Rigidization Force'),1);

    %Gets what the user input on Z-axis
    column_int = get(handles.Z_axis_selection,'value');
    for i = 1:length(handles.figdata)
            param_i = handles.figdata(i).params(2);

               
            smt_title = {title_1;lim_str;...
                ['Smoothing Function of ' ...
                handles.data.parameters(param_i).string]};
            %GIVES DATES
            smt_data = handles.data.filtered_data(on_off_vals);
            y_label_txt = handles.data.full_params{param_i};
    end
 
    for i = 1:length(smt_data)
        is_num = all(~isnan(smt_data{i}),2);
        smt_data{i} = smt_data{i}(is_num,[time_col,pos_col,cur_col,abs_cur_col,rate_col,abs_rate_col]);
    end

    Spinoff_Plot_Smoothing_Function(smt_data,spin_dates,smt_title, y_label_txt, column_int)

    
%     'FFT',...
%     'Standard Deviation Trend',...
%     'Trend of the Slopes',...
%     'Difference of 2 Trend Lines',...
%     'Temperature vs Event',...
%     'Trend vs Temperature Scatter Plot',...
%     'Temperature Compensated Trend Line',...
%     'Smoothing Function',...
end


% --- Executes on button press in Overview_button.
%Changes view of 3D plots to X-Y-Z plane
function Overview_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Overview_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
handles = goto_view([-25,40],handles);
guidata(handles.output,handles)


% --- Executes on button press in Profile_button.
%Changes view of 3D plots to X-Z plane
function Profile_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Profile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
handles = goto_view([0,0],handles);
guidata(handles.output,handles)


% --- Executes on button press in Trending_view_button.
%Changes view of 3D plots to X-Y plane
function Trending_view_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Trending_view_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
handles = goto_view([90,0],handles);
guidata(handles.output,handles)


% --- Executes on button press in Birds_eye_button.
%Changes view of 3D plots to Y-Z plane
function Birds_eye_button_Callback(hObject, eventdata, handles)
%% hObject    handle to Birds_eye_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
if handles.settings.enhanced_view
    handles = goto_view([90,-90],handles);
else
    handles = goto_view([0,90],handles);
end
guidata(handles.output,handles)



% --- Executes on button press in Load_session_menu.
%Loads added session
function Load_session_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Load_session_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = disable_filtering(handles);
handles = disable_3d_plotting(handles);

handles.figdata = [];

set(handles.Import_status,'ForegroundColor','blue')
set(handles.Import_status,...
    'String','Status: Loading previous trending session')

if isempty(eventdata)
    [file,directory] = uigetfile([handles.loading.session_dir '\*.mat'],...
        'Select session to load');
    handles.loading.session_dir = directory;
else
    file = eventdata;
    directory = '';
end
if ~file
    update_import_status(handles);
    msgbox('User cancelled loading previous session')
    return
end
try
    load([directory,file])
    handles.data = [];
    handles.data = saved_session;
catch err
    update_import_status(handles);
    msgbox('Loading previous session failed')
    rethrow(err)
end

if exist('loading','var')
    handles.loading = loading;
    handles.settings = settings;
end
% Set default values to non-existent variables to ensure version
% compatibility
if ~isfield(handles.settings,'raw_data_visible')
    handles.settings.raw_data_visible = true;
end
if ~isfield(handles.settings,'enhanced_view')
    handles.settings.enhanced_view = false;
end
if ~isfield(handles.settings,'grids_visible')
    handles.settings.grids_visible = true;
end
if ~isfield(handles.settings,'baselines_visible')
    handles.settings.baselines_visible = false;
end
if ~isfield(handles.settings,'animations_enabled')
    handles.settings.animations_enabled = false;
end
if ~isfield(handles.settings,'current_view')
    handles.settings.current_view = [-25,40];
end
if ~isfield(handles.settings,'reversed_y_axis')
    handles.settings.reversed_y_axis = false;
end
if ~isfield(handles.settings,'colour_scheme')
    handles.settings.colour_scheme = 'Z gradient';
end
if ~isfield(handles.settings,'marker_style')
    handles.settings.marker_style = 'none';
end
if ~isfield(handles.settings,'marker_size')
    handles.settings.marker_size = 6;
end
if ~isfield(handles.settings,'line_style')
    handles.settings.line_style = '-';
end
if ~isfield(handles.settings,'line_width')
    handles.settings.line_width = 0.5;
end

handles = enable_filtering(handles);
handles = reset_filters(handles);
handles.data = saved_session;
if exist('on_off_vals','var')
    set(handles.Event_table,'data',on_off_vals)
else
    handles.data.filtered_dates = {};
end
handles = update_data_limit_list(handles);
handles = update_restrictions_list(handles);

if isempty(eventdata)
    msgbox(sprintf('A previous trending session has been loaded:\n%s',...
        [directory,file]))
end
handles = update_import_status(handles);

guidata(handles.output,handles)



% --- Executes on button press in Save_session_menu.
%Saves a session in G:Drive
function Save_session_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Save_session_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = msgbox('Saving session. Please wait');

on_off_vals = get(handles.Event_table,'data');
saved_session = handles.data;
loading = handles.loading;
for i = 1:length(loading.folders)
    [loading.folders(i).files.rawdata] = deal({});
end
settings = handles.settings;

if ~isempty(eventdata) % Auto saving using filename in "eventdata"
    save(eventdata,'saved_session','on_off_vals','loading','settings')
else % Save using UI
    uisave({'saved_session','on_off_vals','loading','settings'},...
        [handles.loading.session_dir '\Session_name.mat'])
end

close(h)



% --- Executes on button press in Reverse_Y_Axis_Menu.
%Reverses Y-axis order
function Reverse_Y_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Reverse_Y_Axis_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

if strcmpi(get(handles.Reverse_Y_Axis_Menu,'Checked'),'On')
    set(handles.Reverse_Y_Axis_Menu,'Checked','Off')
    handles.settings.reversed_y_axis = false;
else
    set(handles.Reverse_Y_Axis_Menu,'Checked','On')
    handles.settings.reversed_y_axis = true;
end

handles = update_plot_settings(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Z_Gradient_menu.
%Adds a gradient that defines the difference in z-axis reflected in the
%data
function Z_Gradient_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Z_Gradient_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
handles.settings.colour_scheme = 'Z gradient';
handles = update_plot_data(handles);
handles = update_figure_settings_menu(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Solid_Lines_menu.
%Changes plot to have solid lines with different colours based on event date.
function Solid_Lines_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Solid_Lines_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
handles.settings.colour_scheme = 'Solid line';
handles = update_plot_data(handles);
handles = update_figure_settings_menu(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Y_Gradient_menu.
function Y_Gradient_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Y_Gradient_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);
handles.settings.colour_scheme = 'Y gradient';
handles = update_plot_data(handles);
handles = update_figure_settings_menu(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Solid line style menu.
%allows the user to select the type of line; solid, dotted, dashed, dash-dot and no line
function Solid_line_line_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to selected menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

menu_handles = [...
    handles.Solid_line_menu_solid,...
    handles.Solid_line_menu_no_line,...
    handles.Solid_line_menu_dashed,...
    handles.Solid_line_menu_dotted,...
    handles.Solid_line_menu_dash_dot,...
    ];
line_styles = {'-','none','--',':','-. '};

Selected_menu = menu_handles == hObject;

handles.settings.line_style = line_styles{Selected_menu};

set(menu_handles(Selected_menu),'Checked','On')
set(menu_handles(~Selected_menu),'Checked','Off')

handles = update_plot_settings(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Set Line Width menu.
%Width controls the width of the line and varies from 0.5 – 10 units
function Line_width_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to selected menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

menu_handles = [...
    handles.Line_width_menu_05,...
    handles.Line_width_menu_1,...
    handles.Line_width_menu_2,...
    handles.Line_width_menu_3,...
    handles.Line_width_menu_4,...
    handles.Line_width_menu_5,...
    handles.Line_width_menu_6,...
    handles.Line_width_menu_7,...
    handles.Line_width_menu_8,...
    handles.Line_width_menu_9,...
    handles.Line_width_menu_10,...
    handles.Line_width_menu_11,...
    handles.Line_width_menu_12,...
    ];
widths = [0.5,1,2,3,4,5,6,7,8,9,10,11,12];

Selected_menu = menu_handles == hObject;

set(menu_handles(Selected_menu),'Checked','On')
set(menu_handles(~Selected_menu),'Checked','Off')

handles.settings.line_width = widths(Selected_menu);

handles = update_plot_settings(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Solid line marker style menu.
%offers different types of icons to be used for the points. 
%These include: the plus sign (+), circle (o), asterisk (*), point (.), 
%cross (x), square, diamond or none
function Solid_line_marker_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to selected menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

menu_handles = [...
    handles.Solid_line_menu_none,...
    handles.Solid_line_menu_plus,...
    handles.Solid_line_menu_circle,...
    handles.Solid_line_menu_asterisk,...
    handles.Solid_line_menu_point,...
    handles.Solid_line_menu_cross,...
    handles.Solid_line_menu_square,...
    handles.Solid_line_menu_diamond,...
    ];
marker_styles = {'none','+','o','*','.','x','square','diamond'};

Selected_menu = menu_handles == hObject;

set(menu_handles(Selected_menu),'Checked','On')
set(menu_handles(~Selected_menu),'Checked','Off')

handles.settings.marker_style = marker_styles{Selected_menu};

handles = update_plot_settings(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on button press in Set Marker Size menu.
%allows the user to change the size of the icon from 2 – 48 units
function Marker_width_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to selected menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = update_figure_handles(handles);

menu_handles = [...
    handles.Marker_width_menu_2,...
    handles.Marker_width_menu_4,...
    handles.Marker_width_menu_5,...
    handles.Marker_width_menu_6,...
    handles.Marker_width_menu_7,...
    handles.Marker_width_menu_8,...
    handles.Marker_width_menu_10,...
    handles.Marker_width_menu_12,...
    handles.Marker_width_menu_18,...
    handles.Marker_width_menu_24,...
    handles.Marker_width_menu_48,...
    ];
sizes = [2,4,5,6,7,8,10,12,18,24,48];

Selected_menu = menu_handles == hObject;

set(menu_handles(Selected_menu),'Checked','On')
set(menu_handles(~Selected_menu),'Checked','Off')

handles.settings.marker_size = sizes(Selected_menu);

handles = update_plot_settings(handles);

guidata(get(handles.Figure_settings_menu,'parent'),handles)



% --- Executes on Add_automatic_session_menu button press.
%Saves current session in \\vsgroups\Groups\SE\Trending 
%\Raw_data_and_Trending_GUI\Automatic

function Add_auto_session_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Add_automatic_session_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get name of session to save
sesh_name = inputdlg(['Specify the name of this Automatic session as '...
    'it should appear in the GUI''s dropdown list'],'Auto Session Name');
if isempty(sesh_name)
    msgbox('Adding Autosession aborted')
    return
end

response = questdlg([sprintf(['Confirm the creation of a new automatic '...
    'session named:\n%s\n\nExtracting data from:'],sesh_name{:}),...
    sprintf('\n\n%s',handles.loading.folders.name)],...
    'Confirm New Auto-Session','Confirm','Cancel','Confirm');

if ~strcmp(response,'Confirm')
    msgbox('Adding Autosession aborted')
    return
end

% Save session
filename = ['\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI'...
    '\Automatic_Updating_Scripts\Automatic_Sessions\',sesh_name{:},'.mat'];
Save_session_menu_Callback(hObject, filename, handles)



% --- Executes on Remove_automatic_session_menu button press.
%Removes previously saved automatic sessions from directory. 
function Remove_auto_session_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Remove_automatic_session_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get name of session to remove
session_names = get(handles.Automatic_Sessions_list,'string');
indexed_names = strcat(cellstr(num2str((1:length(...
    session_names))')),{': '},session_names);

dlg_response = inputdlg(['Specify the number of the Automatic session '...
    'to remove:',sprintf('\n%s',indexed_names{:})],'Auto Session Name');

if isempty(dlg_response)
    msgbox('Autosession removal terminated by user')
    return
elseif length(str2double(dlg_response{1})) ~= 1
    msgbox(['Exactly one number must be entered. Autosession removal '...
        'terminated.'])
    return
elseif str2double(dlg_response{1}) < 1 || ...
        str2double(dlg_response{1}) > length(session_names)
    msgbox('Invalid number. Autosession removal terminated.')
    return
end

session_i = str2double(dlg_response{1});

response = questdlg(sprintf(['Confirm the removal of the automatic '...
    'session named:\n%s'],session_names{session_i}),...
    'Confirm Auto-Session Removal','Confirm','Cancel','Confirm');

if ~strcmp(response,'Confirm')
    msgbox('Autosession Removal aborted')
    return
end

% Remove session
filename = ['\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI'...
    '\Automatic_Updating_Scripts\Automatic_Sessions\',...
    session_names{session_i},'.mat'];
delete(filename)
msgbox('Removal confirmed. The session will not appear in future GUI usages')



% --- Executes on Add_automatic_3D_plot_menu button press.
%Saves a 3D plot in \\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI\Automatic_Updating_Scripts\Automatic_Plots 
%and offers the option to email it (not fully integrated, still need MDA email account and script to send emails through it). 
function Add_auto_3d_plot_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Add_automatic_3D_plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get name of session to save
sesh_name = inputdlg(['Specify the name of the Automatic Plot''s '...
    'session (instructions)'],'Auto Plot Name');
if isempty(sesh_name)
    msgbox('Adding Auto Plot aborted')
    return
end

% Get email address to send results to
email_addresses = inputdlg({sprintf(['Specify the email addresses to '...
    'send the updated plots to. Hit cancel to skip this feature\n'...
    'Recipient 1:']),'Recipient 2:','Recipient 3:','Recipient 4:',...
    'Recipient 5:'},'Email Recipients');
email_addresses = email_addresses(~cellfun(@isempty,email_addresses));

figparams = {handles.figdata.params};

response = questdlg([sprintf(['Confirm the creation of a new automatic '...
    'plot named:\n%s\nThat updates the %i open figures, emails the '...
    'results to:\n'],sesh_name{:},length(figparams)),sprintf(...
    '\n%s',email_addresses{:}),sprintf('\n\nExtracting data from:'),...
    sprintf('\n%s',handles.loading.folders.name)],...
    'Confirm New Auto-Plotting Session','Confirm','Cancel','Confirm');

if ~strcmp(response,'Confirm')
    msgbox('Adding Auto Plot aborted')
    return
end

handles.data.auto_plot.emails = email_addresses;
handles.data.auto_plot.figparams = figparams;
% Save session
filename = ['\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI'...
    '\Automatic_Updating_Scripts\Automatic_Plots\',sesh_name{:},'.mat'];
Save_session_menu_Callback(hObject, filename, handles)



% --- Executes on Remove_automatic_3D_plot_menu button press.
%Removes previously saved automatic plot from directory.
function Remove_auto_3d_plot_menu_Callback(hObject, eventdata, handles)
%% hObject    handle to Remove_automatic_3D_plot_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get name of Auto plot to remove
auto_plots = dir(['\\vsgroups\Groups\SE\Trending\'...
    'Raw_data_and_Trending_GUI\Automatic_Updating_Scripts\'...
    'Automatic_Plots']);
auto_plots = auto_plots(~[auto_plots.isdir]);
indexed_names = strcat(cellstr(num2str((1:length(...
    {auto_plots.name}))')),{': '},{auto_plots.name}');

dlg_response = inputdlg(['Specify the number of the Automatic Plot '...
    'to remove:',sprintf('\n%s',indexed_names{:})],'Auto Plot Name');

if isempty(dlg_response)
    msgbox('Auto plot removal terminated by user')
    return
elseif length(str2double(dlg_response{1})) ~= 1
    msgbox(['Exactly one number must be entered. Auto plot removal '...
        'terminated.'])
    return
elseif str2double(dlg_response{1}) < 1 || ...
        str2double(dlg_response{1}) > length(auto_plots)
    msgbox('Invalid number. Auto plot removal terminated.')
    return
end

session_i = str2double(dlg_response{1});

response = questdlg(sprintf(['Confirm the removal of the automatic '...
    'plot named:\n%s'],auto_plots(session_i).name),...
    'Confirm Auto-Plot Removal','Confirm','Cancel','Confirm');

if ~strcmp(response,'Confirm')
    msgbox('Auto plot Removal aborted')
    return
end

% Remove session
filename = ['\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI'...
    '\Automatic_Updating_Scripts\Automatic_Plots\',...
    auto_plots(session_i).name];
delete(filename)
msgbox('Removal confirmed. The session will not appear in future GUI usages')




function Key_press(hObject,eventdata,handles)
%%
if strcmp(eventdata.Modifier,'control')
    switch eventdata.Key
        case 'a'
        case 'c'
        case 'd'
            assignin('base','handles',handles)
            if isfield(handles,'data')
                assignin('base','data',handles.data)
            end
    %         assignin('base','raw_data',handles.raw_data)
            if isfield(handles,'figdata')
                assignin('base','figdata',handles.figdata)
            end
            fprintf('Data retrieved from GUI at: %i/%02i/%02i %02i:%02i:%06.3f\n',clock)
    end
elseif strcmp(eventdata.Modifier,'alt')
    switch eventdata.Key
        case 'c'
            switch handles.settings.colour_scheme
                case 'Z gradient'
                    handles.settings.colour_scheme = 'Solid line';
                    handles = update_plot_settings(handles);
                case 'Solid line'
                    handles.settings.colour_scheme = 'Y gradient';
                    handles = update_plot_data(handles);
                case 'Y gradient'
                    handles.settings.colour_scheme = 'Z gradient';
                    handles = update_plot_data(handles);
                otherwise
                    msgbox('Changing colour scheme failed')
            end
            handles = update_figure_settings_menu(handles);
            guidata(hObject,handles)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SUPPORT FUNCTIONS                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = initialize_auto_sessions(handles)
%%
handles.loading.autosessions = dir(...
    ['\\vsgroups\Groups\SE\Trending\Raw_data_and_Trending_GUI'...
    '\Automatic_Updating_Scripts\Automatic_Sessions']);
handles.loading.autosessions = ...
    handles.loading.autosessions(~[handles.loading.autosessions.isdir]);
session_names = {handles.loading.autosessions.name};
for i = 1:length(session_names)
    if ~strcmp(session_names{i}(end-3:end),'.mat')
        msgbox('Error: Automatic Sessions contains an unexpected file')
        error('Error: Automatic Sessions contains an unexpected file')
    end
    session_names{i} = session_names{i}(1:end-4);
end
set(handles.Automatic_Sessions_list,'string',session_names)


function handles = update_import_status(handles)
%%
% Update the Folder list in the GUI
set(handles.Folder_list,'value',[]);
folder_text = {handles.loading.folders.name};
for i = 1:length(folder_text)
    if length(folder_text{i}) > 60
        folder_text{i} = ['...' folder_text{i}(end-55:end)];
    end
end
set(handles.Folder_list,'Max',length(folder_text)+1)
set(handles.Folder_list,'String',folder_text)

% Determine available files, and how many have already been
% imported/processed
num_imported = 0;
num_new = 0;
num_processed = 0;
num_prev_imported = 0;
num_failed = 0;
if ~isempty(handles.data)
    files_found = false(size(handles.data.eventdates));
else
    files_found = [];
end
for folder_i = 1:length(handles.loading.folders)
    folder = handles.loading.folders(folder_i).name;
    newfiles = dir(folder);
    newfiles = newfiles(~[newfiles.isdir])'; % Exclude folders as potential text files
    
    if isempty(handles.loading.folders(folder_i).files)
        handles.loading.folders(folder_i).files = struct(...
            'name',{},'date',{},'bytes',{},'isdir',{},'datenum',...
            {},'isimported',{},'isprocessed',{},...
            'failed_processing',{},'rawdata',{});
    end
    
    for i = 1:length(newfiles)
        if all(newfiles(i).name(end-3:end)=='.txt')
            full_name = [handles.loading.folders(folder_i).name,...
                '\',newfiles(i).name];
            file_i = find(strcmp(newfiles(i).name,...
                {handles.loading.folders(folder_i).files.name}));
            if length(file_i) > 1
                error(['Error occured during file detection. '...
                    'Multiple files with the same name found.'])
            elseif isempty(file_i)
                file_i = length(handles.loading.folders(folder_i).files)+1;
                newfiles(i).isimported = false;
                newfiles(i).isprocessed = false;
                newfiles(i).failed_processing = false;
                newfiles(i).rawdata = struct('data',{},'textdata',{});
                handles.loading.folders(folder_i).files(file_i) = ...
                    newfiles(i);
            elseif ~all(files_found)
                files_found(strcmp(full_name,handles.data.eventsources...
                    )) = true;
            end
            handles.loading.folders(folder_i).files(file_i).isimported =...
                ~isempty(...
                handles.loading.folders(folder_i).files(file_i).rawdata);
            if isempty(handles.data)
                handles.loading.folders(folder_i).files(...
                    file_i).isprocessed = false;
                handles.loading.folders(folder_i).files(...
                    file_i).failed_processing = false;
            else
                handles.loading.folders(folder_i).files(...
                    file_i).isprocessed = any(strcmp(full_name,...
                    handles.data.eventsources));
            end
        end
    end
    if any([handles.loading.folders(folder_i).files.isimported] & ...
            ~([handles.loading.folders(folder_i).files.isprocessed] |...
            [handles.loading.folders(folder_i).files.failed_processing]))
        error('Found a file that has been imported but not processed')
    end
    num_imported = num_imported + sum(...
        [handles.loading.folders(folder_i).files.isimported]);
    num_new = num_new + sum(...
        ~([handles.loading.folders(folder_i).files.isprocessed] | ...
        [handles.loading.folders(folder_i).files.failed_processing]));
    num_processed = num_processed + sum(...
        [handles.loading.folders(folder_i).files.isprocessed]);
    num_prev_imported = num_prev_imported + sum(...
        [handles.loading.folders(folder_i).files.isprocessed] & ...
        ~[handles.loading.folders(folder_i).files.isimported]) + sum(...
        [handles.loading.folders(folder_i).files.failed_processing] & ...
        ~[handles.loading.folders(folder_i).files.isimported]);
    num_failed = num_failed + sum(...
        [handles.loading.folders(folder_i).files.failed_processing]);
end
num_to_remove = sum(~files_found);
if isempty(num_to_remove)
    num_to_remove = 0;
end

line_1 = sprintf('%i data files imported, %i new, %i old'...
    ,num_imported,num_new,num_prev_imported);
line_2 = sprintf('%i files processed, %i failed, %i to be removed'...
    ,num_processed,num_failed,num_to_remove);
set(handles.Import_status,'string',{line_1;line_2})

if ~num_processed
    set(handles.Import_status,'foregroundcolor','red')
    if isempty(handles.data)
        handles = disable_filtering(handles);
        handles = disable_3d_plotting(handles);
    end
elseif num_to_remove || num_new
    set(handles.Import_status,'foregroundcolor',[1,0.5,0])
    handles = enable_filtering(handles);
else
    set(handles.Import_status,'foregroundcolor',[0,0.5,0])
    handles = enable_filtering(handles);
end



function handles = disable_filtering(handles)
%%
set(handles.Filter_selection,'Enable','off')
set(handles.Filter_status,'Enable','off')
set(handles.Filter_value_Selection,'Enable','off')
set(handles.Add_filter_button,'Enable','off')
set(handles.Apply_filters_button,'Enable','off')
set(handles.Contain_selection,'Enable','off')
set(handles.Restriction_data_selection,'Enable','off')
set(handles.Min_restriction,'Enable','off')
set(handles.Max_restriction,'Enable','off')
set(handles.Add_restriction_button,'Enable','off')
set(handles.Restriction_list,'Enable','off')
set(handles.Remove_restriction_button,'Enable','off')
set(handles.Data_limit_selection,'Enable','off')
set(handles.Min_limit,'Enable','off')
set(handles.Max_limit,'Enable','off')
set(handles.Add_data_limit_button,'Enable','off')
set(handles.Limit_list,'Enable','off')
set(handles.Remove_data_limit_button,'Enable','off')



function handles = reset_filters(handles)
%%
handles.data.full_params = {handles.data.parameters.string};
for i = 1:length(handles.data.full_params)
    dependance = handles.data.parameters(i).unitdep;
    if dependance
        unit_i = handles.data.filters(dependance).values{1};
        unit = handles.data.parameters(i).unit{unit_i};
    else
        unit = handles.data.parameters(i).unit;
    end
    handles.data.full_params{i} = ...
        [handles.data.full_params{i}, ' (' unit ')'];
end

handles.data.active_filter_vals = ones(1,length(handles.data.filters)); % Default filter setting is first value.
handles.data.selected_filter_vals = ones(1,length(handles.data.filters));

set(handles.Filter_selection,'Value',1);
set(handles.Filter_selection,'String',{handles.data.filters.name});

set(handles.Filter_value_Selection,'Value',1);
set(handles.Filter_value_Selection,'String',...
    handles.data.filters(1).options);

set(handles.Restriction_data_selection,'Value',1);
set(handles.Restriction_data_selection,'String',handles.data.full_params);
handles.data.active_restrictions = {};
handles.data.selected_restrictions = {};

set(handles.Data_limit_selection,'Value',1);
set(handles.Data_limit_selection,'String',handles.data.full_params);
handles.data.active_data_limits = {};
handles.data.selected_data_limits = {};

status_str = cell(length(handles.data.filters),1);
for i = 1:length(handles.data.filters)
    status_str{i} = [handles.data.filters(i).name,': ',...
        handles.data.filters(i).options{1}];
end
set(handles.Filter_status,'String',[{'Default values:'};status_str])

handles.figdata = [];
handles.data.trends = [];



function handles = enable_filtering(handles)
%%
set(handles.Filter_selection,'Enable','on')
set(handles.Filter_status,'Enable','on')
set(handles.Filter_value_Selection,'Enable','on')
set(handles.Add_filter_button,'Enable','on')
set(handles.Apply_filters_button,'Enable','on')
set(handles.Contain_selection,'Enable','on')
set(handles.Restriction_data_selection,'Enable','on')
set(handles.Min_restriction,'Enable','on')
set(handles.Max_restriction,'Enable','on')
set(handles.Add_restriction_button,'Enable','on')
set(handles.Restriction_list,'Enable','on')
set(handles.Remove_restriction_button,'Enable','on')
set(handles.Data_limit_selection,'Enable','on')
set(handles.Min_limit,'Enable','on')
set(handles.Max_limit,'Enable','on')
set(handles.Add_data_limit_button,'Enable','on')
set(handles.Limit_list,'Enable','on')
set(handles.Remove_data_limit_button,'Enable','on')



function handles = disable_3d_plotting(handles)
%%
set(handles.Event_table,'Enable','off')

set(handles.X_axis_selection,'Enable','off')
set(handles.Z_axis_selection,'Enable','off')
set(handles.Generate_3d_button,'Enable','off')

set(handles.Raw_data_check,'value',1)

handles = disable_full_plotting(handles);



function handles = enable_3D_plotting(handles)
%%
set(handles.Event_table,'Enable','on')

prev_x_options = get(handles.X_axis_selection,'String');
prev_x_val = get(handles.X_axis_selection,'Value');
set(handles.X_axis_selection,'Value',1);
set(handles.X_axis_selection,'String',handles.data.full_params);
if isequal(prev_x_options,get(handles.X_axis_selection,'String'))
    set(handles.X_axis_selection,'Value',prev_x_val);
end
set(handles.X_axis_selection,'Enable','on')

prev_z_options = get(handles.Z_axis_selection,'String');
prev_z_val = get(handles.Z_axis_selection,'Value');
set(handles.Z_axis_selection,'Value',1);
set(handles.Z_axis_selection,'String',handles.data.full_params);
if isequal(prev_z_options,get(handles.Z_axis_selection,'String'))
    set(handles.Z_axis_selection,'Value',prev_z_val);
end
set(handles.Z_axis_selection,'Enable','on')

set(handles.Generate_3d_button,'Enable','on')



function handles = disable_full_plotting(handles)
%%
set(handles.Max_radio,'Enable','off')
set(handles.Avg_radio,'Enable','off')
set(handles.Min_radio,'Enable','off')
set(handles.Generate_trend_button,'Enable','off')
set(handles.Trend_list,'Value',[])
set(handles.Trend_list,'Enable','off')
set(handles.Trend_list,'string','')
set(handles.Toggle_trend_button,'Enable','off')
set(handles.Clear_trend_button,'Enable','off')
set(handles.Spinoff_button,'Enable','off')
set(handles.Spinoff_selection,'Enable','off')
set(handles.Overview_button,'Enable','off')
set(handles.Profile_button,'Enable','off')
set(handles.Birds_eye_button,'Enable','off')
set(handles.Trending_view_button,'Enable','off')
set(handles.Raw_data_check,'Enable','off')
set(handles.Enhancedviewcheck,'Enable','off')
set(handles.Grid_check,'Enable','off')
set(handles.Baseline_check,'Enable','off')
set(handles.Enable_animations_check,'Enable','off')
set(handles.Figure_settings_menu,'Enable','off')

set(handles.Export_button,'Enable','off')
set(handles.Add_automatic_3D_plot_menu,'Enable','off')



function handles = enable_full_plotting(handles)
%%
set(handles.Max_radio,'Enable','on')
set(handles.Avg_radio,'Enable','on')
set(handles.Min_radio,'Enable','on')
set(handles.Generate_trend_button,'Enable','on')
set(handles.Trend_list,'Enable','on')
set(handles.Toggle_trend_button,'Enable','on')
set(handles.Clear_trend_button,'Enable','on')

set(handles.Spinoff_button,'Enable','on')
set(handles.Spinoff_selection,'Enable','on')

set(handles.Overview_button,'Enable','on')
set(handles.Profile_button,'Enable','on')
set(handles.Birds_eye_button,'Enable','on')
set(handles.Trending_view_button,'Enable','on')
set(handles.Raw_data_check,'Enable','on')
set(handles.Enhancedviewcheck,'Enable','on')
set(handles.Grid_check,'Enable','on')
set(handles.Baseline_check,'Enable','on')
set(handles.Enable_animations_check,'Enable','on')
set(handles.Figure_settings_menu,'Enable','on')
set(handles.Add_automatic_3D_plot_menu,'Enable','on')




function handles = update_restrictions_list(handles)
%%
musts = {' must have a value between',' must not have a value between'};
boxtxt = cell(size(handles.data.selected_restrictions,1),1);
for i = 1:size(handles.data.selected_restrictions,1)
    param_name = handles.data.parameters(handles.data.selected_restrictions{i,2}).string;
    musttxt = musts{2-handles.data.selected_restrictions{i,1}};
    if isempty(handles.data.selected_restrictions{i,3})
        mintxt = 'inf';
    else
        mintxt = num2str(handles.data.selected_restrictions{i,3});
    end
    if isempty(handles.data.selected_restrictions{i,4})
        maxtxt = 'inf';
    else
        maxtxt = num2str(handles.data.selected_restrictions{i,4});
    end
    boxtxt{i} = [param_name,musttxt,' (',mintxt,',',maxtxt,')'];
    bool_exist = false;
    for j = 1:size(handles.data.active_restrictions,1)
        if isequal(handles.data.active_restrictions(j,:),...
                handles.data.selected_restrictions(i,:))
            bool_exist = true;
            break
        end
    end
    if ~bool_exist
        boxtxt{i} = ['(New) ',boxtxt{i}];
    end
end

box_i = length(boxtxt);
for i = 1:size(handles.data.active_restrictions,1)
    bool_is_selected = false;
    for j = 1:size(handles.data.selected_restrictions,1)
        if isequal(handles.data.active_restrictions(i,:),...
                handles.data.selected_restrictions(j,:))
            bool_is_selected = true;
            break
        end
    end
    if ~bool_is_selected
        param_name = handles.data.parameters(handles.data.active_restrictions{i,2}).string;
        musttxt = musts{2-handles.data.active_restrictions{i,1}};
        if isempty(handles.data.active_restrictions{i,3})
            mintxt = 'inf';
        else
            mintxt = num2str(handles.data.active_restrictions{i,3});
        end
        if isempty(handles.data.active_restrictions{i,4})
            maxtxt = 'inf';
        else
            maxtxt = num2str(handles.data.active_restrictions{i,4});
        end
        box_i = box_i+1;
        boxtxt{box_i} = ['(Removed) ',param_name,musttxt,...
            ' (',mintxt,',',maxtxt,')'];
    end
end
if isempty(boxtxt)
    boxtxt = 'None';
    set(handles.Restriction_list,'value',1)
    set(handles.Restriction_list,'max',1)
else
    set(handles.Restriction_list,'max',length(boxtxt)+1)
    set(handles.Restriction_list,'value',[])
end
set(handles.Restriction_list,'String',boxtxt)
if isequal(handles.data.selected_restrictions,handles.data.active_restrictions)
    set(handles.Restriction_list,'ForegroundColor',[0,0.5,0])
else
    set(handles.Restriction_list,'ForegroundColor','red')
end



function handles = update_data_limit_list(handles)
%%
boxtxt = cell(size(handles.data.selected_data_limits,1),1);
for i = 1:size(handles.data.selected_data_limits,1)
    param_name = handles.data.parameters(handles.data.selected_data_limits{i,1}).string;
    if isempty(handles.data.selected_data_limits{i,2})
        mintxt = 'inf';
    else
        mintxt = num2str(handles.data.selected_data_limits{i,2});
    end
    if isempty(handles.data.selected_data_limits{i,3})
        maxtxt = 'inf';
    else
        maxtxt = num2str(handles.data.selected_data_limits{i,3});
    end
    boxtxt{i} = [param_name,': (',mintxt,',',maxtxt,')'];
    bool_exist = false;
    for j = 1:size(handles.data.active_data_limits,1)
        if isequal(handles.data.active_data_limits(j,:),...
                handles.data.selected_data_limits(i,:))
            bool_exist = true;
            break
        end
    end
    if ~bool_exist
        boxtxt{i} = ['(New) ',boxtxt{i}];
    end
end

box_i = length(boxtxt);
for i = 1:size(handles.data.active_data_limits,1)
    bool_is_selected = false;
    for j = 1:size(handles.data.selected_data_limits,1)
        if isequal(handles.data.active_data_limits(i,:),...
                handles.data.selected_data_limits(j,:))
            bool_is_selected = true;
            break
        end
    end
    if ~bool_is_selected
        param_name = handles.data.parameters(handles.data.active_data_limits{i,1}).string;
        if isempty(handles.data.active_data_limits{i,2})
            mintxt = 'inf';
        else
            mintxt = num2str(handles.data.active_data_limits{i,2});
        end
        if isempty(handles.data.active_data_limits{i,3})
            maxtxt = 'inf';
        else
            maxtxt = num2str(handles.data.active_data_limits{i,3});
        end
        box_i = box_i+1;
        boxtxt{box_i} = ['(Removed) ',param_name,...
            ': (',mintxt,',',maxtxt,')'];
    end
end
if isempty(boxtxt)
    set(handles.Limit_list,'value',1)
    boxtxt = 'None';
    set(handles.Limit_list,'max',1)
else
    set(handles.Limit_list,'max',length(boxtxt)+1)
    set(handles.Limit_list,'value',[])
end
set(handles.Limit_list,'String',boxtxt)
if isequal(handles.data.selected_data_limits,handles.data.active_data_limits)
    set(handles.Limit_list,'ForegroundColor',[0,0.5,0])
else
    set(handles.Limit_list,'ForegroundColor','red')
end



function handles = update_plot_settings(handles)
%%
for i = 1:length(handles.figdata)
    set(handles.figdata(i).axesID,'position','default')
    on_off_vals = get(handles.Event_table,'data');
    on_off_vals = [on_off_vals{2:end,2}];
    % Update raw data visibility (3d or 2d)
    if handles.settings.raw_data_visible
        if any(strcmp(handles.settings.colour_scheme,{'Z gradient','Y gradient'}));
            set(handles.figdata(i).h3dlines(on_off_vals),'visible','on')
            set(handles.figdata(i).h3dlines(~on_off_vals),'visible','off')
            set(handles.figdata(i).h2dlines,'visible','off')
            set(handles.figdata(i).axesID,'position','default')
        else
            set(handles.figdata(i).h2dlines(on_off_vals),'visible','on')
            set(handles.figdata(i).h2dlines(~on_off_vals),'visible','off')
            set(handles.figdata(i).h3dlines,'visible','off')
            if any(on_off_vals)
            elseif isempty(handles.data.trends)
                
            elseif ~any([handles.data.trends.is_visible])
                legend(handles.figdata(i).axesID,'off')
            end
        end
    else
        set(handles.figdata(i).h3dlines,'visible','off')
        set(handles.figdata(i).h2dlines,'visible','off')
    end
    % Update trendline and baseline visibility
    for j = 1:length(handles.figdata(i).htrends)
        if handles.settings.baselines_visible
            set(handles.figdata(i).hbaselines,'visible','on')
            set(handles.figdata(i).htrends,'visible','off')
        else
            set(handles.figdata(i).htrends,'visible','on')
            set(handles.figdata(i).hbaselines,'visible','off')
        end
    end
    % Update enhancement visibilities
    if handles.settings.enhanced_view
        set(handles.figdata(i).axesID,'position',[0.06,0.20,0.92,0.72])
        ylabel(handles.figdata(i).axesID,'')
        set(handles.figdata(i).axesID,'yticklabel','')
        set(handles.figdata(i).enhanced_h_lines,'visible','on')
        set(handles.figdata(i).enhanced_h_years,'visible','on')
        set(handles.figdata(i).enhanced_h_y_txt,'visible','on')
    else
        ylabel(handles.figdata(i).axesID,'Event Date')
        set(handles.figdata(i).axesID,'yticklabel',...
            handles.data.filtered_dates)
        set(handles.figdata(i).enhanced_h_lines,'visible','off')
        set(handles.figdata(i).enhanced_h_years,'visible','off')
        set(handles.figdata(i).enhanced_h_y_txt,'visible','off')
    end
    % Update grid visibility
    if handles.settings.grids_visible
        set(handles.figdata(i).axesID,'xgrid','on')
        set(handles.figdata(i).axesID,'zgrid','on')
    else
        set(handles.figdata(i).axesID,'xgrid','off')
        set(handles.figdata(i).axesID,'zgrid','off')
    end  
    % Update line styles
    set(handles.figdata(i).h2dlines,'linestyle',handles.settings.line_style)
    % Update line widths
    set(handles.figdata(i).h2dlines,'LineWidth',handles.settings.line_width)
    % Update marker types
    set(handles.figdata(i).h2dlines,'marker',handles.settings.marker_style)
    % Update Marker sizes
    set(handles.figdata(i).h2dlines,'MarkerSize',handles.settings.marker_size)
    % Change y direction
    if handles.settings.reversed_y_axis
        set(handles.figdata(i).axesID,'YDir','Reverse')
    else
        set(handles.figdata(i).axesID,'YDir','Normal')
    end
    % Update legend (Last to make sure position is good)
    if handles.settings.raw_data_visible && ...
            strcmpi(handles.settings.colour_scheme,'Solid line') && any(on_off_vals) % 2d lines get top priority for legends
        pos = get(handles.figdata(i).axesID,'position');
        legend(handles.figdata(i).h2dlines(on_off_vals),...
            handles.data.filtered_dates(on_off_vals),...
            'location','eastoutside')
        set(handles.figdata(i).axesID,'position',[pos(1:2),0.66,pos(4)])
    elseif isempty(handles.data.trends)
        legend(handles.figdata(i).axesID,'off')
    elseif any([handles.data.trends.is_visible]) && handles.settings.baselines_visible % Baselines + Trendlines get second priority
        hlines = handles.figdata(i).hbaselines;
        hlines = hlines';
        hlines = hlines(:);
        leg_txt = cell(5,length(handles.figdata(i).trend_txt));
        leg_txt(1,:) = handles.figdata(i).trend_txt;
        for j = 1:size(leg_txt,2)
            leg_txt{1,j} = handles.figdata(i).trend_txt{j};
            leg_txt{2,j} = sprintf('Overall Average: %.3f',handles.figdata(i).baseline_vals(j,2));
            leg_txt{3,j} = sprintf('Upper 3%s bound: %.3f','\sigma',handles.figdata(i).baseline_vals(j,1));
            leg_txt{4,j} = sprintf('Lower 3%s bound: %.3f','\sigma',handles.figdata(i).baseline_vals(j,3));
            leg_txt{5,j} = sprintf('Outliers (%i)',handles.figdata(i).baseline_vals(j,4));
        end
        leg_txt = leg_txt(:);
        legend(hlines,leg_txt,'location','northeast')
    elseif any([handles.data.trends.is_visible]) && ~handles.settings.baselines_visible % Trendlines get third priority
        legend([handles.figdata(i).htrends],handles.figdata(i).trend_txt,'location','northeast')
    else
        legend(handles.figdata(i).axesID,'off')
    end
end



function handles = update_figure_handles(handles)
%%
if isfield(handles.figdata,'figID')
    open_figs = ishandle([handles.figdata.axesID]);
    handles.figdata = handles.figdata(open_figs);
end
if isempty(handles.figdata)
    handles = disable_full_plotting(handles);
end
% Rename uncontrolled figures with 'GUI Controlled' in their title
open_figs = findall(0,'type','figure');
for i = 1:length(open_figs)
    figname = get(open_figs(i),'name');
    if strfind(figname,'GUI Controlled')
        if isempty(handles.figdata)
            set(open_figs(i),'name',['Uncontrolled',figname(15:end)])
        elseif ~any(open_figs(i) == [handles.figdata.figID])
            set(open_figs(i),'name',['Uncontrolled',figname(15:end)])
        end
    end
end
        



function handles = update_plot_data(handles)
%% Regenerates all plot data
%% Get dates to display on the y axis, depending on visible trendlines
if ~isempty(handles.data.trends)
    dates = [handles.data.filtered_dates,...
        handles.data.trends([handles.data.trends.is_visible]).dates];
    timevals = zeros(size(dates));
    
    for i = length(dates):-1:1
        if any(strcmp(dates{i},dates(1:i-1)))
            dates(i) = [];
            timevals(i) = [];
        else
            str = dates{i};
            year = str2double(str(1:4));
            if ~isstrprop(str(7),'digit')
                date = str2double(str(6));
                hour = str2double(str(8:9));
                minute = str2double(str(11:12));
            elseif ~isstrprop(str(8),'digit')
                date = str2double(str(6:7));
                hour = str2double(str(9:10));
                minute = str2double(str(12:13));
            else
                date = str2double(str(6:8));
                hour = str2double(str(10:11));
                minute = str2double(str(13:14));
            end
            timevals(i) = year*1000+date+hour/100+minute/10000;
        end
    end
    [~,sort_I] = sort(timevals);
    dates = dates(sort_I);
else
    dates = handles.data.filtered_dates;
end

trends = handles.data.trends;

for i = 1:length(handles.figdata)

    selected_params = handles.figdata(i).params;
    
    txt = handles.data.processed_fcn;
    for j = 1:length(handles.data.filters)
        txt = [txt(1:end),handles.data.filters(j).titletxt{...
            handles.data.active_filter_vals(j)}];
    end
    handles.figdata(i).titletxt = {txt;...
        [handles.data.parameters(selected_params(2)).string ' vs ' ...
        handles.data.parameters(selected_params(1)).string]};
    title(handles.figdata(i).axesID,handles.figdata(i).titletxt)
    set(handles.figdata(i).figID,'name',['GUI Controlled 3D Figure: ' handles.figdata(i).titletxt{1}])

    handles.figdata(i).xaxistxt = ...
        handles.data.full_params{selected_params(1)};
    xlabel(handles.figdata(i).axesID,handles.figdata(i).xaxistxt)

    handles.figdata(i).zaxistxt = ...
        handles.data.full_params{selected_params(2)};
    zlabel(handles.figdata(i).axesID,handles.figdata(i).zaxistxt)
    
    delete(handles.figdata(i).h3dlines)
    delete(handles.figdata(i).h2dlines)
    delete(handles.figdata(i).enhanced_h_lines)
    delete([handles.figdata(i).enhanced_h_years',...
        handles.figdata(i).enhanced_h_y_txt'])
    handles.figdata(i).h3dlines = [];
    handles.figdata(i).h2dlines = [];
    handles.figdata(i).enhanced_h_lines = [];
    handles.figdata(i).enhanced_h_years = [];
    handles.figdata(i).enhanced_h_y_txt = [];
    hold(handles.figdata(i).axesID,'all')
    % Plot filtered 2d and 3d data
    colours = get(handles.figdata(i).axesID,'ColorOrder');
    for j = 1:length(handles.data.filtered_data)
        event_i = find(strcmp(handles.data.eventdates,...
            handles.data.filtered_dates(j)),1,'first');
        colour_i = floor(rem(event_i-1,size(colours,1)))+1;
        xvals = handles.data.filtered_data{j}(:,selected_params(1))';
        yvals = find(strcmp(handles.data.filtered_dates{j},dates))*...
            ones(size(xvals));
        zvals = handles.data.filtered_data{j}(:,selected_params(2))';
        if strcmp(handles.settings.colour_scheme,'Y gradient')
            grad = yvals;
        else
            grad = zvals;
        end
        handles.figdata(i).h3dlines(j) = surface(...
            [xvals;xvals],[yvals;yvals],[zvals;zvals],[grad;grad],...
            'parent',handles.figdata(i).axesID,...
            'facecol','no','edgecol','interp','linew',2);
        handles.figdata(i).h2dlines(j) = plot3(...
            handles.figdata(i).axesID,xvals,yvals,zvals,...
            'color',colours(colour_i,:));
    end
    set(handles.figdata(i).h2dlines,'visible','off')
    
    set(handles.figdata(i).axesID,'ytick',1:length(dates))
    set(handles.figdata(i).axesID,'ylim',[0.5,length(dates)+0.5])
    set(handles.figdata(i).axesID,'yticklabel',dates)
    
    % Generate Trend lines
    delete(handles.figdata(i).htrends)
    handles.figdata(i).htrends = [];
    delete(handles.figdata(i).hbaselines)
    handles.figdata(i).hbaselines = [];
    handles.figdata(i).baseline_vals = [];
    colours = get(handles.figdata(i).axesID,'ColorOrder');
    colours = colours([1,2,4:end],:); % Exclude red, the colour of the outliers

    % Plot trendlines and their baselines
    for j = 1:length(trends)
        if trends(j).is_visible
            colour_i = floor(rem(j-1,size(colours,1)))+1;
            params = handles.figdata(i).params;
            xvals = trends(j).vals(params(2),:,params(1));
            yvals = zeros(size(xvals));
            for k = 1:length(yvals)
                yvals(k) = find(strcmp(trends(j).dates{k},dates));
            end
            zvals = trends(j).vals(params(2),:,params(2));
            handles.figdata(i).htrends(end+1) = plot3(handles.figdata(i).axesID,...
                xvals,yvals,zvals,'-*','color',colours(colour_i,:));
            
            avg_x = mean(xvals); % plot baselines in the middle of the trend on the x axis
            outliers = true;
            inliers = true(1,length(zvals));
            while outliers
                avg_z = mean(zvals(inliers));
                upper_bnd = avg_z + 3*std(zvals(inliers));
                lower_bnd = avg_z - 3*std(zvals(inliers));
                new_inliers = zvals < upper_bnd & zvals > lower_bnd;
                if all(new_inliers == inliers)
                    outliers = false;
                else
                    inliers = new_inliers;
                end
            end
            handles.figdata(i).baseline_vals(end+1,:) = [upper_bnd,avg_z,lower_bnd,sum(~inliers)];
            
            handles.figdata(i).hbaselines(end+1,1) = plot3(handles.figdata(i).axesID,... % Inliers
                xvals(inliers),yvals(inliers),zvals(inliers),'-*','color',colours(colour_i,:));
            handles.figdata(i).hbaselines(end,2) = plot3(handles.figdata(i).axesID,... % Average value
                [avg_x,avg_x],[yvals(1),yvals(end)],[avg_z,avg_z],'--','color',[0,0.5,0]);
            handles.figdata(i).hbaselines(end,3) = plot3(handles.figdata(i).axesID,... % Upper bound
                [avg_x,avg_x],[yvals(1),yvals(end)],[upper_bnd,upper_bnd],'--m');
            handles.figdata(i).hbaselines(end,4) = plot3(handles.figdata(i).axesID,... % Lower bound
                [avg_x,avg_x],[yvals(1),yvals(end)],[lower_bnd,lower_bnd],'--c');
            if any(~inliers)
                handles.figdata(i).hbaselines(end,5) = plot3(handles.figdata(i).axesID,...
                    xvals(~inliers),yvals(~inliers),zvals(~inliers),'*r');
            else
                handles.figdata(i).hbaselines(end,5) = plot3(handles.figdata(i).axesID,...
                    nan,nan,nan,'*r');
            end
        end
    end
    
    % Set axis limits to active data limits if they exist
    x_param = handles.figdata(i).params(1);
    z_param = handles.figdata(i).params(2);
    prev_view = get(handles.figdata(i).axesID,'View');
    set(handles.figdata(i).axesID,'View',[0,0]);
    xlim(handles.figdata(i).axesID,'auto')
    zlim(handles.figdata(i).axesID,'auto')
    xlims = get(handles.figdata(i).axesID,'xlim');
    zlims = get(handles.figdata(i).axesID,'zlim');
    set(handles.figdata(i).axesID,'View',prev_view);
    for j = 1:size(handles.data.active_data_limits,1)
        if handles.data.active_data_limits{j,1} == x_param
            if ~isempty(handles.data.active_data_limits{j,2})
                xlims(1) = handles.data.active_data_limits{j,2};
            end
            if ~isempty(handles.data.active_data_limits{j,3})
                xlims(2) = handles.data.active_data_limits{j,3};
            end
        elseif handles.data.active_data_limits{j,1} == z_param
            if ~isempty(handles.data.active_data_limits{j,2})
                zlims(1) = handles.data.active_data_limits{j,2};
            end
            if ~isempty(handles.data.active_data_limits{j,3})
                zlims(2) = handles.data.active_data_limits{j,3};
            end
        end
    end
    zlims(2) = zlims(2) + eps;
    % Add Enhancements
    yticks = get(handles.figdata(i).axesID,'YTick');
    xpos = xlims(1) - 0.01*(xlims(2)-xlims(1));
    zpos = zlims(1) - 0.01*(zlims(2)-zlims(1));
    
    handles.figdata(i).enhanced_h_y_txt = ... % Rotated y axis labels
        text(xpos*ones(size(yticks)),yticks,zpos*ones(size(yticks)),...
        dates,'Parent',handles.figdata(i).axesID,...
        'HorizontalAlignment','right','rotation',90,'FontSize',8);
    
    years = zeros(size(dates));
    yvals = zeros(length(dates),2);
    for j = 1:length(dates)
        years(j) = str2double(dates{j}(1:4));
        if ~any(years(j)==yvals(:,2))
            yvals(j,:) = [j-0.5,years(j)];
        end
    end
    yvals = yvals(yvals(:,1)~=0,:);
    
    for j = 2:size(yvals,1)
        handles.figdata(i).enhanced_h_lines(j-1) = ... % Add year separation grid
            plot3([xlims(1),xlims],...
            ones(1,3)*yvals(j,1),[zlims,zlims(2)],':','color','black',...
            'Parent',handles.figdata(i).axesID);
    end

    xpos = (xlims(1)+0.05*(xlims(2)-xlims(1)))*ones(size(yvals,1),1);
    ypos = (yvals(:,1)+[yvals(2:end,1);length(dates)+0.5])/2;
    zpos = (zlims(1)+0.05*(zlims(2)-zlims(1)))*ones(size(yvals,1),1);
    handles.figdata(i).enhanced_h_years = ... % Add year labels
        text(xpos,ypos,zpos,num2str(yvals(:,2)),...
        'rotation',90,'Parent',handles.figdata(i).axesID);
    
    set(handles.figdata(i).axesID,'xlim',xlims)
    set(handles.figdata(i).axesID,'zlim',zlims)
end

% Update Active Trend Lines box (and legend text for trends)
set(handles.Trend_list,'Max',length(handles.data.trends)+1)
trend_txt = cell(size(trends));
event_groups = cell(size(trends));
for i = 1:length(trends)
    group_found = false;
    for j = 1:length(event_groups)
        if isequal(trends(i).dates,event_groups{j})
            group_found = true;
        end
    end
    if ~group_found
        event_groups{i} = trends(i).dates;
    end
end
event_groups = event_groups(~cellfun(@isempty,event_groups));

for i = 1:length(trends)
    % Generate string identifier based on trend settings
    trend_str = [trends(i).type,' - '];
    % Display unique filters
    for j = 1:length(trends(i).filters)
        disp_filter = false;
        for k = 1:length(trends)
            if trends(k).filters(j) ~= trends(i).filters(j)
                disp_filter = true;
            end
        end
        if handles.data.active_filter_vals(j) ~= trends(i).filters(j)
            disp_filter = true;
        end
        if disp_filter
            trend_str = [trend_str,handles.data.filters(j).name,':',handles.data.filters(j).options{trends(i).filters(j)},', '];
        end
    end
    % Display unique limits
    for j = 1:size(trends(i).limits,1)
        disp_lim = false;
        for k = find(1:length(trends) ~= i)
            lim_found = false;
            for l = 1:size(trends(k).limits,1)
                if isequal(trends(i).limits(j,:),trends(k).limits(l,:))
                    lim_found = true;
                end
            end
            if ~lim_found
                disp_lim = true;
            end
        end
        if disp_lim
            min_txt = num2str(trends(i).limits{j,2});
            if isempty(min_txt)
                min_txt = 'inf';
            end
            max_txt = num2str(trends(i).limits{j,3});
            if isempty(max_txt)
                max_txt = 'inf';
            end
            trend_str = [trend_str,handles.data.parameters(trends(i).limits{j,1}).string,...
                ' limits:(',min_txt,',',max_txt,'), '];
        end
    end
    % Display unique restrictions
    for j = 1:size(trends(i).restrictions,1)
        disp_rest = false;
        for k = find(1:length(trends) ~= i)
            rest_found = false;
            for l = 1:size(trends(k).restrictions,1)
                if isequal(trends(i).restrictions(j,:),trends(k).restrictions(l,:))
                    rest_found = true;
                end
            end
            if ~rest_found
                disp_rest = true;
            end
        end
        if disp_rest
            musts = {' must not contain (',' must contain ('};
            min_txt = num2str(trends(i).restrictions{j,3});
            if isempty(min_txt)
                min_txt = 'inf';
            end
            max_txt = num2str(trends(i).restrictions{j,4});
            if isempty(max_txt)
                max_txt = 'inf';
            end
            trend_str = [trend_str,handles.data.parameters(trends(i).restrictions{j,2}).string,...
                musts{1+trends(i).restrictions{j,1}},min_txt,',',max_txt,'), '];
        end
    end
    trend_str = trend_str(1:end-2);
    if length(event_groups) > 1
        for j = 1:length(event_groups)
            if isequal(event_groups{j},trends(i).dates)
                group_i = j;
            end
        end
        trend_str = [trend_str, sprintf(', Group %i (%i events)',group_i,length(trends(i).dates))];
    end
    if ~trends(i).is_visible
        trend_str = ['(Hidden) - ',trend_str];
    end
    trend_txt(i) = {trend_str};
end

any_visible_trends = false;
for j = 1:length(handles.figdata)
    if ~isempty(trends)
        handles.figdata(j).trend_txt = trend_txt([trends.is_visible]);
        any_visible_trends = any([trends.is_visible]);
    end
end

if any_visible_trends
    set(handles.Export_button,'Enable','on')
else
    set(handles.Export_button,'Enable','off')
end

set(handles.Trend_list,'string',trend_txt)

handles = update_plot_settings(handles);



function handles = goto_view(az_el,handles)
%%
handles.anim_frames = 50;
handles.fram_dur = 0.02;

if get(handles.Enable_animations_check,'Value')
    az = cell(size(handles.figdata));
    el = az;
    for i = 1:length(handles.figdata)
        [az{i},el{i}] = view(handles.figdata(i).axesID);
        az{i} = linspace(az{i},az_el(1),handles.anim_frames);
        el{i} = linspace(el{i},az_el(2),handles.anim_frames);
    end
    for i = 1:size(az{i},2)
        for j = 1:length(handles.figdata)
            view(handles.figdata(j).axesID,az{j}(i),el{j}(i))
        end
        pause(handles.fram_dur)
    end
else
    for i = 1:length(handles.figdata)
        view(handles.figdata(i).axesID,az_el(1),az_el(2))
    end
end
handles.settings.current_view = az_el;



function handles = update_figure_settings_menu(handles)
%%
menus = [handles.Z_Gradient_menu,handles.Solid_Lines_menu,handles.Y_Gradient_menu];
switch handles.settings.colour_scheme
    case 'Z gradient'
        set(menus(1),'Checked','on')
        set(menus(2:3),'Checked','off')
        set(handles.Solid_line_settings_menu,'enable','off')
    case 'Solid line'
        set(menus(2),'Checked','on')
        set(menus([1,3]),'Checked','off')
        set(handles.Solid_line_settings_menu,'enable','on')
    case 'Y gradient'
        set(menus(3),'Checked','on')
        set(menus(1:2),'Checked','off')
        set(handles.Solid_line_settings_menu,'enable','off')
    otherwise
        msgbox('Something went wrong')
        error('Failed to update colour menu')
end

