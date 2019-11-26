% Thomas Pearson April 21, 2015
% This function returns processed SPDM LEE data for use with the
% Trending_GUI V2
%
%
% Inputs: 
%       raw_data            Data that has been imported using the
%                           Trending_GUI.
%
%       bool_query_version  Flag to return only the version info of this
%                           function
%
% Outputs:
%       processed_data      A structure recognized by the Trending_GUI with
%                           the following fields:
%           .process_fcn            'SPDM LEE'
%           .process_fcn_version    The current version of this function
%
%           .filters(j).name        Name of filter (string)
%           .filters(j).options{m}  Text options for that filter (cell)
%           .filters(j).values{m}   Corresponding integer values (cell of
%                                   vectors e.g.{[1,2],[1],[2]}). Each
%                                   vector represents a list of possible
%                                   values in filterdata{i}
%           .filters(j).titletxt{m} Text to be added to the title of plots
%                                   (string)
%
%           .filterdata{i}          nxj array containing values that
%                                   correspond to filters.values. "n" is
%                                   the number of data points in the i'th
%                                   file
%
%           .parameters(k).string   Name of the parameter to use in the
%                                   axis labels (include units).
%           .parameters(k).unit     Cell containing possible units for the
%                                   paramater, or string if unitdep is 0.
%           .parameters(k).unitdep  Index of filter the unit is dependant
%                                   on. 0 means there is no dependance
%
%           .parameterdata{i}       nxk array containing the processed data
%                                   of the i'th file
%
%           .eventdates{i}          The event date (cell of strings)
%
%           .temperature(i).name{t} A 1xt cell containing the names of the
%                                   possible temperatures.
%           .temperature(i).value(t) A vector containing the corresponding
%                                   temperature values. Each event will
%                                   have a single temperature per mechanism
%
% The expected (mandatory) PUIs in raw_data are as follows:
%{
Active_Effector_Active_Mechanism	C_ActEffActiveMech
Active_Effector_Mechanism_Position	C_ActEffMechPos
Active_Effector_Measured_Motor_Current	C_ActEffMeasMotorCur
Active_Effector_Derived_Motor_Rate	C_ActEffDeriveMotRat
Computed_Active_Eff_Rigidization_Force_Derived_Torque	C_ActEffRigFrcTorq
SPDM_LEE_SACU_1_Capture_With_Latch	CMRC13SW0PTQJ
SPDM_LEE_SACU_2_Capture_With_Latch	CMRC13SW0Q5UJ
SPDM_LEE_SACU_1_LCS_State	CMRC13SW0PSTJ
SPDM_LEE_SACU_2_LCS_State	CMRC13SW0Q4XJ
Computed_Command_Status_ID_10Hz	C_CmdStatusID
%}
%
% The optional PUIs that are also useful are:
%{
SPDM_Arm_1_Shell_Temp	CMRC13SW0Q0PT
SPDM_Arm_1_Latches_Temp	CMRC13SW0Q0LT
SPDM_Arm_1_Carriage_Temp	CMRC13SW0Q0KT
SPDM_Arm_1_Snare_MM_Motor_Temp	CMRC13SW0Q0OT
SPDM_Arm_1_Rigidize_MM_Motor_Temp	CMRC13SW0Q0NT
SPDM_Arm_1_Latch_MM_Motor_Temp	CMRC13SW0Q0MT
SPDM_Arm_2_Shell_Temp	CMRC13SW0QCTT
SPDM_Arm_2_Latches_Temp	CMRC13SW0QCPT
SPDM_Arm_2_Carriage_Temp	CMRC13SW0QCOT
SPDM_Arm_2_Snare_MM_Motor_Temp	CMRC13SW0QCST
SPDM_Arm_2_Rigidize_MM_Motor_Temp	CMRC13SW0QCRT
SPDM_Arm_2_Latch_MM_Motor_Temp	CMRC13SW0QCQT
%}
%
% PUIs that are not currently used, but may be useful for future
% implementations are:
%{
Computed_Command_Response_ID_10Hz	C_CmdResponseID
%}

function processed_data = process_SPDM_LEE_data(raw_data,bool_query_version)

processed_data.processed_fcn = 'SPDM LEE'; % Start of the title for all plots
processed_data.processed_fcn_version = 2.1;

if bool_query_version
    return
end

processed_data.filters = struct(...
    'name',{'Mechanism',...
            'Active SACU',...
            'Grapple Fixture Type',...
            'Nominal/Off-Nominal'},...
    'options',{ {'Snare','Carriage','Latch','Umbilical'},...
                {'All','SACU 1','SACU 2'},...
                {'All','PDGF','FRGF'},...
                {'Both','Nominal','Off-Nominal'}},...
    'values',{  {1,2,3,3},... % These also act as indeces to select units for position. For now, Latching is indistinguishable from Mating
                {[0,1,2],1,2},...% Zero indicates unknown state
                {[0,1,2],1,2},...% Zero indicates unknown state
                {[0,1],1,0}},...% 1 = nominal, 0 = off-nominal
    'titletxt',{{' Snare',' Rigidize',' Latch',' Mate'},...
                {'',' (SACU 1)',' (SACU 2)'},...
                {'',' on PDGFs',' on FRGFs'},...
                {'',' - Nominal Captures',' - Off-Nominal Captures'}});

mech_filter_col = 1;
str_filter_col = 2;
GF_filter_col = 3;
nom_filter_col = 4;

processed_data.parameters = struct(...
    'string',{...
        'Time',...
        'Position',...
        'Current',...
        'Absolute Current',...
        'Motor Rate',...
        'Absolute Motor Rate',...
        'Rigidization Force'},...
    'unit',{...
        's',...
        {'Degrees','Inches','Inches','Inches'},...
        'A',...
        'A',...
        'Rad/s',...
        'Rad/s',...
        'lbf'},...
    'unitdep',{...
        0,...
        mech_filter_col,... % Only position units are dependant on mechanism.
        0,...
        0,...
        0,...
        0,...
        0});

time_param_col = 1;
pos_param_col = 2;
cur_param_col = 3;
abs_cur_param_col = 4;
rate_param_col = 5;
abs_rate_param_col = 6;
load_param_col = 7;


num_files = length(raw_data);
bool_good_data = true(1,num_files);

temperature_names = {'Shell Temperature','Carriage Temperature','Latch Temperature',...
    'SMM Motor Temperature','RMM Motor Temperature','LMM Motor Temperature'};

temperature_PUIs = {'CMRC13SW0Q0PT','CMRC13SW0Q0KT','CMRC13SW0Q0LT',...
    'CMRC13SW0Q0OT','CMRC13SW0Q0NT','CMRC13SW0Q0MT';...
    'CMRC13SW0QCTT','CMRC13SW0QCOT','CMRC13SW0QCPT',...
    'CMRC13SW0QCST','CMRC13SW0QCRT','CMRC13SW0QCQT'};

processed_data.parameterdata = {};
processed_data.filterdata = {};
processed_data.eventdates = {};
processed_data.temperature = {};

for i = 1:num_files
    n_pts = size(raw_data(i).data,1);
    n_filters = length(processed_data.filters);
    
    processed_data.filterdata{i} = zeros(n_pts,n_filters);
    
    % Determine operational string
    SACU1_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0PSTJ'),1,'first');
    SACU2_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0Q4XJ'),1,'first');
    if length([SACU1_col,SACU2_col]) ~= 2
        error('Error: SPDM_LEE_SACU_1_LCS_State: CMRC13SW0PSTJ and/or SPDM_LEE_SACU_2_LCS_State: CMRC13SW0Q4XJ were not found')
    end
    SACU1_active = any(raw_data(i).data(:,SACU1_col) == 2);
    SACU2_active = any(raw_data(i).data(:,SACU2_col) == 2);
    if SACU1_active && SACU2_active
        error('Error: Both strings (prime and redun) are active')
    elseif ~(SACU1_active || SACU2_active)
        processed_data.filterdata{i}(:,str_filter_col) = false;
        active_SACU = [];
    else
        active_SACU = find([SACU1_active,SACU2_active],1,'first');
        processed_data.filterdata{i}(:,str_filter_col) = active_SACU;
    end
    
    % Extract temperatures
    for j = 1:length(temperature_PUIs)
        if any(strcmpi(temperature_PUIs{active_SACU,j},raw_data(i).textdata(1,2:end)))
            temp_col = find(strcmp(raw_data(i).textdata(1,2:end),temperature_PUIs{active_SACU,j}),1,'first');
            processed_data.temperature(i).name{j} = temperature_names{j};
            processed_data.temperature(i).value(j) = ... % replace 999... type data with the mean temp value
                mean(raw_data(i).data(raw_data(i).data(:,temp_col)<10000,temp_col));
        else
            processed_data.temperature(i).name{j} = '';
            processed_data.temperature(i).value(j) = NaN;
        end
    end
    
    % Determine GF type (FRGF/PDGF)
    if active_SACU == 1
        gf_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0PTQJ'),1,'first');
    elseif active_SACU == 2
        gf_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0Q5UJ'),1,'first');
    else
        error('Something went wrong with active SACUs')
    end
    
    if isempty(gf_col)
        error('Error: SPDM_LEE_SACU_1/2_Capture_With_Latch: CMRC13SW0PTQJ/CMRC13SW0Q5UJ was not found')
    end
    pdgf_active = sum(raw_data(i).data(:,gf_col) == 1);
    frgf_active = sum(raw_data(i).data(:,gf_col) == 0);
    if pdgf_active > frgf_active
        processed_data.filterdata{i}(:,GF_filter_col) = 1;
    elseif frgf_active
        processed_data.filterdata{i}(:,GF_filter_col) = 2;
    else
        processed_data.filterdata{i}(:,GF_filter_col) = false;
    end
    
    % Fix data due to changing puis over different dates
    str = char(raw_data(i).textdata(3,1));
    year = str2double(str(1:4));
    if ~isstrprop(str(7),'digit')
        date = str2double(str(6));
    elseif ~isstrprop(str(8),'digit')
        date = str2double(str(6:7));
    else
        date = str2double(str(6:8));
    end
    time = year*1000+date;
    time_col = ...
        [find(strcmpi(raw_data(i).textdata(1,2:end),'Unix_Time'),1,'first'),...
        find(strcmpi(raw_data(i).textdata(1,2:end),'Duration'),1,'first')];
    if length(time_col) ~= 1
        error('Error with Time PUI (Unix_Time/Duration)')
    end
    
    if time < 2014250 % From 2008 day 50 until 2014 day 250
        
        cur_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMeasMotorCur'),1,'first');
        rate_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffDeriveMotRat'),1,'first');
        load_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffRigFrcTorq'),1,'first');
        
        pos_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMechPos'),1,'first');
        mech_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffActiveMech'),1,'first');
        
        if length([pos_col,mech_col]) ~= 2
            error('Error with position or mechanism PUI')
        end
        
        position_data = raw_data(i).data(:,pos_col);
        mech_data = raw_data(i).data(:,mech_col) + 1;
        
        raw_data(i).data(:,cur_col) = raw_data(i).data(:,cur_col)*16;
        
    else % 2014 day 250 to present
        
        cur_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMeasMotorCur'),1,'first');
        rate_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffDeriveMotRat'),1,'first');
        load_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffRigFrcTorq'),1,'first');
        
        pos_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMechPos'),1,'first');
        mech_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffActiveMech'),1,'first');
        
        if length([pos_col,mech_col]) ~= 2
            error('Error with position or mechanism PUI')
        end
        
        position_data = raw_data(i).data(:,pos_col);
        mech_data = raw_data(i).data(:,mech_col) + 1;
        
    end
    
    if length([cur_col,rate_col,load_col]) ~= 3
        error('Error with current, motor rate, or rig force PUIs')
    end
    
    processed_data.filterdata{i}(:,mech_filter_col) = mech_data;
    
    % Extract good_rows
    bool_good_rows = all([position_data,mech_data],2) & ... % position and mechanism must be non-zero
        all([raw_data(i).data(:,[cur_col,rate_col,load_col]),position_data,mech_data] < 100000,2); % No 999999999 type data
    
    processed_data.parameterdata{i}(:,...
        [pos_param_col,time_param_col,cur_param_col,rate_param_col,load_param_col]) = ...
        [position_data,raw_data(i).data(:,[time_col,cur_col,rate_col,load_col])];
    processed_data.parameterdata{i} = processed_data.parameterdata{i}(bool_good_rows,:);
    processed_data.filterdata{i} = processed_data.filterdata{i}(bool_good_rows,:);
    
    % Remove faulty data where the active mechanism pui is 0.
    if ~isempty(processed_data.parameterdata{i})
        mech_data = processed_data.filterdata{i}(:,mech_filter_col);
        mech_starts = [1;find(diff(mech_data))+1]; % indices after change
        mech_ends = [mech_starts(2:end,:)-1;length(mech_data)];
        mech_lengths = diff([mech_starts;length(mech_data)+1]);
        
        problem_mechs = mech_lengths < 10 & (mech_data(mech_starts)==1); % Segments with active mech == 0 that are shorter than 10 data points
        
        problem_starts = mech_starts(problem_mechs,:);
        problem_ends = mech_ends(problem_mechs,:);
        
        bool_good_rows2 = true(size(mech_data));
        for j = 1:length(problem_starts)
            if problem_starts(j) == 1 || problem_ends(j) == length(mech_data) % Remove if at the start or end of file
                bool_good_rows2(problem_starts(j):problem_ends(j)) = false;
            elseif mech_data(problem_starts(j)-1) == mech_data(problem_ends(j)+1) % Fix if its in the middle of good data
                mech_data(problem_starts(j):problem_ends(j)) = mech_data(problem_starts(j)-1);
            else % Remove otherwise
                bool_good_rows2(problem_starts(j):problem_ends(j)) = false;
            end
        end
        processed_data.filterdata{i}(:,mech_filter_col) = mech_data;
        processed_data.parameterdata{i} = processed_data.parameterdata{i}(bool_good_rows2,:);
        processed_data.filterdata{i} = processed_data.filterdata{i}(bool_good_rows2,:);
    end
    
    processed_data.parameterdata{i}(:,[abs_cur_param_col,abs_rate_param_col]) = ...
        abs(processed_data.parameterdata{i}(:,[cur_param_col,rate_param_col]));
    
    if ~isempty(processed_data.parameterdata{i})
        % Filter for nominal captures only. Nominal defined as:
        % 1. Completed, not aborted or rejected
        % 2. The rigidization load does not drop by more than 100 lbf (during
        %       latching or rigidization)
        % 3. For rigidization, the peak load must be between 1210 and 1272

        stat_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_CmdStatusID'),1,'first');

        rejected_i = find(raw_data(i).data(:,stat_col)==4,1,'last');
        aborted_i = find(raw_data(i).data(:,stat_col)==7,1,'last');
        completed_i = find(raw_data(i).data(:,stat_col)==11,1,'last');

        status_i_s = [0,0,0];
        if ~isempty(rejected_i)
            status_i_s(1) = rejected_i;
        end
        if ~isempty(aborted_i)
            status_i_s(2) = aborted_i;
        end
        if ~isempty(completed_i)
            status_i_s(3) = completed_i;
        end
        final_stat_i = max(status_i_s);
        processed_data.filterdata{i}(:,nom_filter_col) = true;

        if final_stat_i == rejected_i
            processed_data.filterdata{i}(:,nom_filter_col) = false;
        elseif final_stat_i == aborted_i
            processed_data.filterdata{i}(:,nom_filter_col) = false;
        elseif final_stat_i == completed_i
        else
            error('Completion status not discovered')
        end

        forces = processed_data.parameterdata{i}(:,load_param_col);

        [peak_force,peak_i] = max(forces);
        if any(peak_force-forces(peak_i:end) > 100)
            fprintf('File %i of %i sees a drop of %f lbf and will be considered off-nominal.\n'...
                ,i,num_files,max(peak_force-forces(peak_i:end)))
            processed_data.filterdata{i}(:,nom_filter_col) = false;
        end

        % Remove events with gap greater than 3 seconds?
        time_vals = processed_data.parameterdata{i}(:,time_param_col);
        if any(diff(time_vals) > 3)
            fprintf('File %i of %i has a data gap larger than 3 seconds and will be considered off-nominal.\n'...
                ,i,num_files)
            processed_data.filterdata{i}(:,nom_filter_col) = false;
        end

        % Peak force must be between 1210 and 1272 lbf (For rigidization
        % only)
        peak_force = max(processed_data.parameterdata{i}(processed_data.filterdata{i}(:,mech_filter_col)==2,load_param_col));
        if ~isempty(peak_force)
            if peak_force < 1210 || peak_force > 1272
                fprintf('File %i of %i has a peak force of %f lbf and will be considered off-nominal.\n'...
                    ,i,num_files,peak_force)
                processed_data.filterdata{i}(processed_data.filterdata{i}(:,mech_filter_col)==2,nom_filter_col) = false;
            end
        end
    
        % Update status
        processed_data.eventdates{i} = raw_data(i).textdata{find(bool_good_rows,1,'first')+1,1};
        processed_data.parameterdata{i}(:,time_param_col) = ...
            processed_data.parameterdata{i}(:,time_param_col) - ...
            processed_data.parameterdata{i}(1,time_param_col); % Start time = 0
    else
        bool_good_data(i) = false;
        fprintf('File %i of %i removed. No good data found\n',...
            i,num_files)
    end
end

processed_data.parameterdata = processed_data.parameterdata(bool_good_data);
processed_data.filterdata = processed_data.filterdata(bool_good_data);
processed_data.eventdates = processed_data.eventdates(bool_good_data);
processed_data.temperature = processed_data.temperature(bool_good_data);

% Sort raw_data
time = zeros(length(processed_data.eventdates),1);

for i = 1:length(time)
    str = processed_data.eventdates{i};
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
    time(i) = year*1000+date+hour/100+minute/10000;
end

[~,sort_I] = sort(time);

processed_data.parameterdata = processed_data.parameterdata(sort_I);
processed_data.filterdata = processed_data.filterdata(sort_I);
processed_data.eventdates = processed_data.eventdates(sort_I);
processed_data.temperature = processed_data.temperature(sort_I);

end
