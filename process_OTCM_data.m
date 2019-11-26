% This function returns processed OTCM data for use with the 
% Trending_GUI V2
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
%           .process_fcn            'SPDM OTCM'
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
Computed_Meq_ID_(10Hz)	C_MeqID
%}
%
% The optional PUIs that are also useful are:
%{
SPDM_Arm_2_OTCM_Housing_Temp_2	CMRC13SW0QG5T
SPDM_Arm_2_OTCM_Housing_Temp_1	CMRC13SW0QG4T
SPDM_Arm_1_OTCM_Housing_Temp_2	CMRC13SW0Q41T
SPDM_Arm_1_OTCM_Housing_Temp_1	CMRC13SW0Q40T
SPDM_Arm_2_OTCM_Torquer_MM_Rate	CMRC13SW0QC8R
SPDM_Arm_2_OTCM_Umbilical_MM_Rate	CMRC13SW0QC6R
SPDM_Arm_2_OTCM_Advance_MM_Rate	CMRC13SW0QC4R
SPDM_Arm_2_OTCM_OEU_Gripper_MM_Rate	CMRC13SW0QC2R
SPDM_Arm_1_OTCM_Torquer_MM_Rate	CMRC13SW0Q04R
SPDM_Arm_1_OTCM_Umbilical_MM_Rate	CMRC13SW0Q02R
SPDM_Arm_1_OTCM_Advance_MM_Rate	CMRC13SW0Q00R
SPDM_Arm_1_OTCM_OEU_Gripper_MM_Rate	CMRC13SW0PZYR
%}
%
% PUIs that are not currently used, but may be useful for future
% implementations are:
%{
MSS_Payload_Status_OCS_SPDM_Arm_2_OTCM	CMRC13SW0FFYJ	10565626	UNKNOWN	UNKNOWN	0	0
MSS_Payload_Status_OCS_SPDM_Arm_1_OTCM	CMRC13SW0FFWJ	10565628	UNKNOWN	UNKNOWN	0	0
Computed_Command_Response_ID_(10Hz)	C_CmdResponseID	1281759	N/A	N/A	0	0
Computed_Command_Status_ID_(10Hz)	C_CmdStatusID	1281760	N/A	N/A	0	0
SPDM_Arm_1_OTCM_Umbilical_Mated	CMRC13SW0PZPJ	10563810	UNKNOWN			
SPDM_Arm_1_OTCM_Torquer_MM_In_Motion	CMRC13SW0PZIJ	10563817	UNKNOWN			
SPDM_Arm_1_OTCM_Advance_MM_In_Motion	CMRC13SW0PZHJ	10563818	UNKNOWN			
SPDM_Arm_1_OTCM_OEU_Gripper_MM_In_Motion	CMRC13SW0PZGJ	10563819	UNKNOWN			
SPDM_Arm_1_OTCM_Umbilical_MM_In_Motion	CMRC13SW0PZFJ	10563820	UNKNOWN			
SPDM_1_OTCM_Umbilical_Mechanism_Absolute_Position	CMRC13SW14QGH	10562614	CM	CM	0	0
SPDM_1_OTCM_Advance_Mechanism_Absolute_Position	CMRC13SW14QFH	10562615	CM	CM	0	0
SPDM_1_OTCM_Torquer_Mechanism_Absolute_Position	CMRC13SW14QEF	10562616	REV	REV	0	0
SPDM_1_OTCM_Gripper_Mechanism_Absolute_Position	CMRC13SW14QDH	10562617	CM	CM	0	0
SPDM_Arm_1_OTCM_TMM_Motor_Winding_Temp	CMRC13SW0Q45T	10563657	CNT			
SPDM_Arm_1_OTCM_UMM_Motor_Winding_Temp	CMRC13SW0Q44T	10563658	CNT			
SPDM_Arm_1_OTCM_AMM_Motor_Winding_Temp	CMRC13SW0Q43T	10563659	CNT			
SPDM_Arm_1_OTCM_Gripper_MM_Motor_Winding_Temp	CMRC13SW0Q42T	10563660	CNT			
SPDM_Arm_1_OTCM_OEU_Derived_Torque	CMRC13SW0Q07G	10563792	N-M			
SPDM_Arm_2_OTCM_Umbilical_Mated	CMRC13SW0QBTJ	10563400	UNKNOWN			
SPDM_Arm_2_OTCM_Torquer_MM_In_Motion	CMRC13SW0QBMJ	10563407	UNKNOWN			
SPDM_Arm_2_OTCM_Advance_MM_In_Motion	CMRC13SW0QBLJ	10563408	UNKNOWN			
SPDM_Arm_2_OTCM_OEU_Gripper_MM_In_Motion	CMRC13SW0QBKJ	10563409	UNKNOWN			
SPDM_Arm_2_OTCM_Umbilical_MM_In_Motion	CMRC13SW0QBJJ	10563410	UNKNOWN			
SPDM_2_OTCM_Umbilical_Mechanism_Absolute_Position	CMRC13SW14QKH	10562610	CM	CM	0	0
SPDM_2_OTCM_Advance_Mechanism_Absolute_Position	CMRC13SW14QJH	10562611	CM	CM	0	0
SPDM_2_OTCM_Torquer_Mechanism_Absolute_Position	CMRC13SW14QIF	10562612	REV	REV	0	0
SPDM_2_OTCM_Gripper_Mechanism_Absolute_Position	CMRC13SW14QHH	10562613	CM	CM	0	0
SPDM_Arm_2_OTCM_TMM_Motor_Winding_Temp	CMRC13SW0QG9T	10563247	CNT			
SPDM_Arm_2_OTCM_UMM_Motor_Winding_Temp	CMRC13SW0QG8T	10563248	CNT			
SPDM_Arm_2_OTCM_AMM_Motor_Winding_Temp	CMRC13SW0QG7T	10563249	CNT			
SPDM_Arm_2_OTCM_Gripper_MM_Motor_Winding_Temp	CMRC13SW0QG6T	10563250	CNT			
SPDM_Arm_2_OTCM_OEU_Derived_Torque	CMRC13SW0QCBG	10563382	N-M			
%}
%

function processed_data = process_OTCM_data(raw_data,bool_query_version)

processed_data.processed_fcn = 'SPDM OTCM'; % Start of the title for all plots
processed_data.processed_fcn_version = 2.1;

if bool_query_version
    return
end

processed_data.filters = struct(...
    'name',{'OTCM 1/2',...
            'Mechanism',...
            'Advancer Slopes',...
            'Gripper Secondary Peaks'},...
    'options',{ {'Both','OTCM 1','OTCM 2'},...
                {'Gripper','Advancer','Torquer','Umbilical'},...
                {'Off','On'},...
                {'Off','On'}},...
    'values',{  {[1,2],1,2},...
                {1,2,3,4},... % These also act as indeces to select units for position
                {[0,1],1},...
                {[0,1],1}},...
    'titletxt',{{' 1 and 2',' 1',' 2'},...
                {' Gripper',' Advancer',' Torquer',' Umbilical'}...
                {'',' Slopes'}...
                {'',' Secondary Peaks'}});

OTCM_filter_col = 1;
mech_filter_col = 2;
advance_slopes_filter_col = 3;
grip_peaks_filter_col = 4;

processed_data.parameters = struct(...
    'string',{...
        'Time',...
        'Position',...
        'Current',...
        'Absolute Current',...
        'Motor Rate',...
        'Absolute Motor Rate'},...
        'unit',{...
        's',...
        {'cm','cm','Rev','cm'},...
        'A',...
        'A',...
        'Rad/s',...
        'Rad/s'},...
        'unitdep',{...
        0,...
        mech_filter_col,... % Only position units are dependant on mechanism.
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


num_files = length(raw_data);
bool_good_data = true(1,num_files);

temperature_names = {'OTCM 1 Housing Temp 1',...
    'OTCM 1 Housing Temp 2','OTCM 2 Housing Temp 1',...
    'OTCM 2 Housing Temp 2'};

temperature_PUIs = {'CMRC13SW0Q40T','CMRC13SW0Q41T',...
    'CMRC13SW0QG4T','CMRC13SW0QG5T'};

processed_data.parameterdata = {};
processed_data.filterdata = {};
processed_data.eventdates = {};
processed_data.temperature = {};

for i = 1:num_files
    n_pts = size(raw_data(i).data,1);
    n_filters = length(processed_data.filters);
    n_params = length(processed_data.parameters);
    
    processed_data.filterdata{i} = zeros(n_pts,n_filters);
    
    % Determine operational SACU
    OTCM_state_col = find(strcmp(raw_data(i).textdata(1,2:end),'C_MeqID'),1,'first');
    if isempty(OTCM_state_col)
        error('Error: C_MeqID not found')
    end
    OTCM_1_state = any(raw_data(i).data(:,OTCM_state_col) == 6);
    OTCM_2_state = any(raw_data(i).data(:,OTCM_state_col) == 7);
    if OTCM_1_state && ~OTCM_2_state
        Operational_OTCM = 1;
    elseif OTCM_2_state && ~OTCM_1_state
        Operational_OTCM = 2;
    elseif OTCM_2_state && OTCM_1_state
        Operational_OTCM = 0;
        rate1_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0PZYR'),1,'first'); % Gripper 1 motor rate
        rate2_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0QC2R'),1,'first'); % Gripper 2 motor rate
        if length([rate1_col,rate2_col]) == 2
            if any(0.1 < abs(raw_data(i).data(:,rate1_col)) & raw_data(i).data(:,rate1_col) < 10000)
                Operational_OTCM = 1;
                if any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                    Operational_OTCM = 3; % Both have non-zero gripper motor rates
                end
            elseif any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                Operational_OTCM = 2;
            end
        end
        if ~Operational_OTCM % If Gripper rate did not discover operational OTCM, try Advancer
            rate1_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0Q00R'),1,'first'); % Advancer 1 motor rate
            rate2_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0QC4R'),1,'first'); % Advancer 2 motor rate
            if length([rate1_col,rate2_col]) == 2
                if any(0.1 < abs(raw_data(i).data(:,rate1_col)) & raw_data(i).data(:,rate1_col) < 10000)
                    Operational_OTCM = 1;
                    if any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                        Operational_OTCM = 3; % Both advancer rates are non-zero
                    end
                elseif any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                    Operational_OTCM = 2;
                end
            end
        end
        if ~Operational_OTCM % If Advancer rate did not discover operational OTCM, try Torquer
            rate1_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0Q04R'),1,'first'); % Torquer 1 rate
            rate2_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0QC8R'),1,'first'); % Torquer 2 rate
            if length([rate1_col,rate2_col]) == 2
                if any(0.1 < abs(raw_data(i).data(:,rate1_col)) & raw_data(i).data(:,rate1_col) < 10000)
                    Operational_OTCM = 1;
                    if any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                        Operational_OTCM = 3;
                    end
                elseif any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                    Operational_OTCM = 2;
                end
            end
        end
        if ~Operational_OTCM % If Torquer rate did not discover operational OTCM, try Umbilical
            rate1_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0Q02R'),1,'first'); % Umbilical 1 rate
            rate2_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0QC6R'),1,'first'); % Umbilical 2 rate
            if length([rate1_col,rate2_col]) == 2
                if any(0.1 < abs(raw_data(i).data(:,rate1_col)) & raw_data(i).data(:,rate1_col) < 10000)
                    Operational_OTCM = 1;
                    if any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                        Operational_OTCM = 3;
                    end
                elseif any(0.1 < abs(raw_data(i).data(:,rate2_col)) & raw_data(i).data(:,rate2_col) < 10000)
                    Operational_OTCM = 2;
                end
            end
        end
        if Operational_OTCM == 3
            msgbox(sprintf('Error: Both OTCM still determined to be active based on motor rates'))
        end
    else
        fprintf('WARNING: No OTCMs Operational in this file\n')
        bool_good_data(i) = false;
    end
    processed_data.filterdata{i}(:,OTCM_filter_col) = Operational_OTCM;
    mech_col = find(strcmp(raw_data(i).textdata(1,2:end),'C_ActEffActiveMech'),1,'first');
    if isempty(mech_col)
        error('Error with PUI: C_ActEffActiveMech')
    end
    processed_data.filterdata{i}(:,mech_filter_col) = raw_data(i).data(:,mech_col)+1;
    
    bad_mechs = processed_data.filterdata{i}(:,mech_filter_col) > 100;
    for j = find(bad_mechs)
        prev_i = find(~bad_mechs(1:j),1,'last');
        next_i = find(~bad_mechs(j:end),1,'first');
        if processed_data.filterdata{i}(prev_i,mech_filter_col) == ...
                processed_data.filterdata{i}(next_i,mech_filter_col)
            processed_data.filterdata{i}(j,mech_filter_col) = ...
                processed_data.filterdata{i}(next_i,mech_filter_col);
        end
    end
    
    % Extract temperatures
    for j = 1:length(temperature_PUIs)
        if any(strcmpi(temperature_PUIs{j},raw_data(i).textdata(1,2:end)))
            temp_col = find(strcmp(raw_data(i).textdata(1,2:end),temperature_PUIs{j}),1,'first');
            processed_data.temperature(i).name{j} = temperature_names{j};
            processed_data.temperature(i).value(j) = ... % replace 999... type data with the mean temp value
                mean(raw_data(i).data(raw_data(i).data(:,temp_col)<10000,temp_col));
        else
            processed_data.temperature(i).name{j} = '';
            processed_data.temperature(i).value(j) = NaN;
        end
    end
    
    bool_good_rows = true(n_pts,1);
    % Extract good_rows
    time_col = ...
        [find(strcmpi(raw_data(i).textdata(1,2:end),'Unix_Time'),1,'first'),...
        find(strcmpi(raw_data(i).textdata(1,2:end),'Duration'),1,'first')];
    if length(time_col) ~= 1
        error('Error with Time PUI (Unix_Time/Duration)')
    end
    pos_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMechPos'),1,'first');
    cur_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMeasMotorCur'),1,'first');
    rate_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffDeriveMotRat'),1,'first');
    if length([pos_col,cur_col,rate_col]) ~= 3
        error('Error with PUIs (C_ActEffMechPos/C_ActEffMeasMotorCur/C_ActEffDeriveMotRat/C_ActEffRigFrcTorq)')
    end
    if bool_good_data(i)
        for j = 1:n_pts
            if any(~raw_data(i).data(j,pos_col)) || ... % Position must be non zero
                    any(raw_data(i).data(j,[pos_col,cur_col,rate_col]) == 999999999) % position, current, and rate must be good data
                bool_good_rows(j) = false;
            end
        end
%         % Remove isolated good data points (bad data before and after)
%         bool_good_rows([false;diff(diff(bool_good_rows))==-2;false]) = false;
        
        if any(bool_good_rows)
            processed_data.parameterdata{i} = zeros(sum(bool_good_rows),n_params);
            processed_data.parameterdata{i}(:,...
                [time_param_col,pos_param_col,cur_param_col,rate_param_col])...
                = raw_data(i).data(bool_good_rows,[...
                time_col,pos_col,cur_col,rate_col]);
            processed_data.parameterdata{i}(:,time_param_col) = ...
                processed_data.parameterdata{i}(:,time_param_col) - ...
                processed_data.parameterdata{i}(1,time_param_col); % Set t=0 at first point
            processed_data.filterdata{i} = processed_data.filterdata{i}(bool_good_rows,:);
            % Remove faulty data where the active mechanism pui is 0.
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
%{
            
            bool_bad_rows = [false; ... % Active mechanism changes twice in 2 data points
                (diff(processed_data.filterdata{i}(:,mech_filter_col)) & ...
                [~diff(abs(diff(processed_data.filterdata{i}(:,mech_filter_col))));false])];
            processed_data.parameterdata{i} = processed_data.parameterdata{i}(~bool_bad_rows,:);
            processed_data.filterdata{i} = processed_data.filterdata{i}(~bool_bad_rows,:);
            if processed_data.filterdata{i}(end-1,mech_filter_col) ~= processed_data.filterdata{i}(end,mech_filter_col) % Active mechanism changes at the end
                processed_data.parameterdata{i} = processed_data.parameterdata{i}(1:end-1,:);
                processed_data.filterdata{i} = processed_data.filterdata{i}(1:end-1,:);
            end
%}
            processed_data.eventdates{i} = raw_data(i).textdata{find(bool_good_rows,1,'first')+1,1};

            processed_data.parameterdata{i}(:,[abs_cur_param_col,abs_rate_param_col]) = ...
                abs(processed_data.parameterdata{i}(:,[cur_param_col,rate_param_col]));
        else
            bool_good_data(i) = false;
        end
    end
    
    if bool_good_data(i)
        % Advancer slopes data
        processed_data.filterdata{i}(:,advance_slopes_filter_col) = false;
        ad_bool = processed_data.filterdata{i}(:,mech_filter_col)==2;
        cur_data = processed_data.parameterdata{i}(ad_bool,cur_param_col);
        pos_data = processed_data.parameterdata{i}(ad_bool,pos_param_col);
        [min_val,min_i] = min(cur_data);
        last_i = [];
        if min_val
            last_i = find(cur_data < min_val + 0.1 & (1:length(cur_data) > min_i)',1,'first')-1; % Filter 1: Last current is after peak and less than 0.1A greater
        end
        if last_i
            first_i = find(cur_data > -0.20 & (1:length(cur_data) < last_i)',1,'last'); % Filter 2: First current is last one above -0.2A
            if isempty(first_i)
                first_i = 1;
            end
            first_i = find(pos_data > pos_data(first_i) + 0.5 & (1:length(cur_data) > first_i)',1,'first'); % Filter 3: Shave 0.5 cm off the front
            last_i = find(pos_data < pos_data(last_i) - 0.1 & (1:length(cur_data) < last_i)',1,'last'); % Filter 4: Shave 0.1 cm off the end

            if abs(pos_data(first_i)-pos_data(last_i)) > 0.7 % Final filter: Ensure there is 0.7 cm worth of data
                offset = find(ad_bool,1,'first')-1;
                processed_data.filterdata{i}(first_i+offset:last_i+offset,advance_slopes_filter_col) = true;
            end
        end
        
        % Gripper secondary peaks data
        processed_data.filterdata{i}(:,grip_peaks_filter_col) = false;
        gr_bool = processed_data.filterdata{i}(:,mech_filter_col)==1;
        cur_data = processed_data.parameterdata{i}(gr_bool,cur_param_col);
        pos_data = processed_data.parameterdata{i}(gr_bool,pos_param_col);
        rate_data = processed_data.parameterdata{i}(gr_bool,rate_param_col);
        [min_val,first_i] = min(cur_data);
        if ~isempty(cur_data)
            if min_val < -0.8 && pos_data(first_i) > 1.5 && ~any(rate_data<-200)
                first_i = find([diff(cur_data);0]<0 & (1:length(cur_data) > first_i)',1,'first');
    %             [~,max_i] = max(cur_data(first_i:end));
    %             first_i = first_i + max_i -1;
                offset = size(processed_data.parameterdata{i},1) - sum(gr_bool);
                if any(cur_data(first_i:end) > -0.8) && any(cur_data(first_i:end) < -0.4)
                    processed_data.filterdata{i}(first_i+offset:length(cur_data)+offset,grip_peaks_filter_col) = true;
                end
            end
        end
    else
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
