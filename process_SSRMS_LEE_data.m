% Thomas Pearson March 19, 2015
% This function returns processed SSRMS LEE data for use with the
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
%           .process_fcn            'SSRMS LEE'
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
From GMT 2008/050 to present:

Active_Effector_Active_Mechanism	C_ActEffActiveMech
Active_Effector_Mechanism_Position	C_ActEffMechPos
Active_Effector_Measured_Motor_Current	C_ActEffMeasMotorCur
Active_Effector_Derived_Motor_Rate	C_ActEffDeriveMotRat
Computed_Active_Eff_Rigidization_Force_Derived_Torque	C_ActEffRigFrcTorq
Operating_Base_OCS_SSRMS	CMRC13SW018DU
Capture_with_Latch_SSRMS	CMRC13SW00P0U
Fst_State_Ssrms_Primary_OCS	CMRC13SW018ZU
Fst_State_Ssrms_Redundant_OCS	CMRC13SW0190U
Computed_Command_Status_ID_10Hz	C_CmdStatusID
Lee_Rigidize_Control_Type_SSRMS CMRC13SW15DGU

From GMT 2002/250 to GMT 2008/050:

Active_Effector_Active_Mechanism	C_ActEffActiveMech
Active_Effector_Mechanism_Position	C_ActEffMechPos
Active_Effector_Measured_Motor_Current	C_ActEffMeasMotorCur
Active_Effector_Derived_Motor_Rate	C_ActEffDeriveMotRat
Active_Effector_Rigid_Force	C_ActEffRigidForce
Operating_Base_OCS_SSRMS	CMRC13SW018DU
Capture_with_Latch_SSRMS	CMRC13SW00P0U
Fst_State_Ssrms_Primary_OCS	CMRC13SW018ZU
Fst_State_Ssrms_Redundant_OCS	CMRC13SW0190U
Computed_Command_Status_ID_10Hz	C_CmdStatusID

From start of life (2001) to GMT 2002/250

Measured_Motor_Current_SSRMS_LEE	CMRC13SW00P4C
Derived_Motor_Rate_SSRMS_LEE	CMRC13SW00P3R
SSRMS_LEE_Rigidization_Force	CMRC13SW00OUG
Latch_Position_SSRMS_LEE	CMRC13SW00P7H
Carriage_Position_SSRMS_LEE	CMRC13SW00P6H
Snare_Position_SSRMS_LEE	CMRC13SW00P8H
Snare_in_Motion_SSRMS_LEE	CMRC13SW00OHU
Carriage_in_Motion_SSRMS_LEE	CMRC13SW00OKU
Latches_in_Motion_SSRMS_LEE	CMRC13SW00OPU
Umbilical_In_Motion_SSRMS_LEE	CMRC13SW00ORU
Operating_Base_OCS_SSRMS	CMRC13SW018DU
Capture_with_Latch_SSRMS	CMRC13SW00P0U
Fst_State_Ssrms_Primary_OCS	CMRC13SW018ZU
Fst_State_Ssrms_Redundant_OCS	CMRC13SW0190U
Computed_Command_Status_ID_10Hz	C_CmdStatusID

%}
%
% The optional PUIs that are also useful are:
%{
Shell_Temperature_SSRMS	CMRC13SW00U3T
Latches_Temperature_SSRMS	CMRC13SW00TZT
Carriage_Temperature_SSRMS	CMRC13SW00TYT
SMM_Motor_Temperature_SSRMS	CMRC13SW00U2T
RMM_Motor_Temperature_SSRMS	CMRC13SW00U1T
LMM_Motor_Temperature_SSRMS	CMRC13SW00U0T
Calibration_Status_Ssrms_Tip_Lee_OCS	CMRC13SW01BYU
Execution_Run_Speed_SSRMS	CMRC13SW00OYU
Ssrms_Gf_Id_SSRMS_Evnt	CMRC13SW01VTU
Latches_Disengaged_SSRMS_LEE	CMRC13SW00OMU
SSRMS_Base_Location_OCS_SSRMS	CMRC13SW018CU
%}
%
% PUIs that are not currently used, but may be useful for future
% implementations are:
%{
Computed_Command_Response_ID_10Hz	C_CmdResponseID
%}

function processed_data = process_SSRMS_LEE_data(raw_data,bool_query_version)

processed_data.processed_fcn = 'SSRMS LEE'; % Start of the title for all plots
processed_data.processed_fcn_version = 2.2;

if bool_query_version
    return
end

processed_data.filters = struct(...
    'name',{'LEE A/B',...
            'Mechanism',...
            'Prime/Redundant',...
            'Grapple Fixture Type',...
            'Nominal/Off-Nominal',...
            'Rigidize Control Type'},...
    'options',{ {'Both','LEE A','LEE B'},...
                {'Snare','Carriage','Latch','Umbilical'},...
                {'All','Prime','Redundant'},...
                {'All','PDGF','FRGF'},...
                {'Both','Nominal','Off-Nominal'},...
                {'Force','Position'}},...
    'values',{  {[1,2],1,2},...
                {1,2,3,3},... % These also act as indeces to select units for position. For now, Latching is indistinguishable from Mating
                {[0,1,2],1,2},...% Zero indicates unknown state
                {[0,1,2],1,2},...% Zero indicates unknown state
                {[0,1],1,0},...% 1 = nominal, 0 = off-nominal
                {0,1}},...%0 =force, 1 = position
    'titletxt',{{' A and B',' A',' B'},...
                {' Snare',' Rigidize',' Latch',' Mate'},...
                {'',' (Prime)',' (Redundant)'},...
                {'',' on PDGFs',' on FRGFs'},...
                {'',' - Nominal Captures',' - Off-Nominal Captures'},...
                {'Force','Position'}});

LEE_filter_col = 1;
mech_filter_col = 2;
str_filter_col = 3;
GF_filter_col = 4;
nom_filter_col = 5;
rig_filter_col = 6;

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

temperature_PUIs = {'CMRC13SW00U3T','CMRC13SW00TYT','CMRC13SW00TZT',...
    'CMRC13SW00U2T','CMRC13SW00U1T','CMRC13SW00U0T'};

processed_data.parameterdata = {};
processed_data.filterdata = {};
processed_data.eventdates = {};
processed_data.temperature = {};

for i = 1:num_files
    n_pts = size(raw_data(i).data,1);
    n_filters = length(processed_data.filters);
    processed_data.filterdata{i} = zeros(n_pts,n_filters); %creates array of zeroes 
    
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
    
%     %Rigidize Control Type - Force & Position
%     force_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW15DGU'),1,'first');
%     position_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW15DGU'),1,'first');
%     if length([force_col,position_col]) ~= 2
%         error('Error: Lee_Rigidize_Control_Type_SSRMS: CMRC13SW15DGU was not found')
%     end
%     force_active = any(raw_data(i).data(:,force_active) == 0);
%     position_active = any(raw_data(i).data(:,position_col) == 1);
%     if force_active && position_col
%         error('Error: Both force and position are active')
%     elseif ~(force_active || position_col)
%         processed_data.filterdata{i}(:,rig_filter_col) = false;
%     else
%         processed_data.filterdata{i}(:,rig_filter_col) = find([force_active,position_active],1,'first');
%     end
        
    % Determine operational LEE
    lee_check_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW018DU'),1,'first');
    if isempty(lee_check_col)
        error('Error: Operating_Base_OCS_SSRMS: CMRC13SW018DU was not found')
    end
    lee_a_active = any(raw_data(i).data(:,lee_check_col) == 5);
    lee_b_active = any(raw_data(i).data(:,lee_check_col) == 0);
    if lee_a_active && lee_b_active
        warning('Warning: Both LEEs are active')
        lee_a_active = find(raw_data(i).data(:,lee_check_col) == 5);
        lee_b_active = find(raw_data(i).data(:,lee_check_col) == 0);
        fprintf(['LEE A has %i data points indicating its active and '...
            'LEE B has %i data points indicating its active\n'],...
            length(lee_a_active),length(lee_b_active))
        if length(lee_a_active) > length(lee_b_active)
            lee_a_active = true;
            lee_b_active = false;
            fprintf('Assuming LEE A is active')
        else
            lee_a_active = false;
            lee_b_active = true;
            fprintf('Assuming LEE B is active')
        end
        processed_data.filterdata{i}(:,LEE_filter_col) = find([lee_a_active,lee_b_active],1,'first');
        pause(1)
    elseif ~(lee_a_active || lee_b_active)
        processed_data.filterdata{i}(:,str_filter_col) = false;
    else
        processed_data.filterdata{i}(:,LEE_filter_col) = find([lee_a_active,lee_b_active],1,'first');
    end
    
    % Determine operational string
    prime_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW018ZU'),1,'first');
    redun_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW0190U'),1,'first');
    if length([prime_col,redun_col]) ~= 2
        error('Error: Fst_State_Ssrms_Primary_OCS: CMRC13SW018ZU and/or Fst_State_Ssrms_Redundant_OCS: CMRC13SW0190U were not found')
    end
    prime_active = any(raw_data(i).data(:,prime_col) == 5);
    redun_active = any(raw_data(i).data(:,redun_col) == 5);
    if prime_active && redun_active
        error('Error: Both strings (prime and redun) are active')
    elseif ~(prime_active || redun_active)
        processed_data.filterdata{i}(:,str_filter_col) = false;
    else
        processed_data.filterdata{i}(:,str_filter_col) = find([prime_active,redun_active],1,'first');
    end
    
    % Determine GF type (FRGF/PDGF)
    gf_col = find(strcmp(raw_data(i).textdata(1,2:end),'CMRC13SW00P0U'),1,'first');
    if isempty(gf_col)
        error('Error: Capture_with_Latch_SSRMS: CMRC13SW00P0U was not found')
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
    
    if time < 2002250 % before 2002 day 250
        
        cur_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00P4C'),1,'first');
        rate_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00P3R'),1,'first');
        load_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00OUG'),1,'first');
        raw_data(i).data(:,load_col) = raw_data(i).data(:,load_col) * 0.224808943; % N to lbf conversion
        
        % Snare_Position_SSRMS_LEE
        snare_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00P8H'),1,'first');
        % Snare_in_Motion_SSRMS_LEE
        snare_motion_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00OHU'),1,'first');
        snare_motion = raw_data(i).data(:,snare_motion_col);
        snare_motion(snare_motion == 999999999) = false;
        snare_motion = logical(snare_motion);
        % Carriage_Position_SSRMS_LEE
        carri_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00P6H'),1,'first');
        % Carriage_in_Motion_SSRMS_LEE
        carri_motion_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00OKU'),1,'first');
        carri_motion = raw_data(i).data(:,carri_motion_col);
        carri_motion(carri_motion == 999999999) = false;
        carri_motion = logical(carri_motion);
        % Latch_Position_SSRMS_LEE
        latch_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00P7H'),1,'first');
        % Latches_in_Motion_SSRMS_LEE
        latch_motion_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00OPU'),1,'first');
        latch_motion = raw_data(i).data(:,latch_motion_col);
        latch_motion(latch_motion == 999999999) = false;
        latch_motion = logical(latch_motion);
        % Umbilical_In_Motion_SSRMS_LEE
        umbil_motion_col = find(strcmpi(raw_data(i).textdata(1,2:end),'CMRC13SW00ORU'),1,'first');
        umbil_motion = raw_data(i).data(:,umbil_motion_col);
        umbil_motion(umbil_motion == 999999999) = false;
        umbil_motion = logical(umbil_motion);
        
        if length([snare_col,snare_motion_col,carri_col,carri_motion_col,latch_col,latch_motion_col,umbil_motion_col]) ~= 7
            error('Error with the 1Hz PUIs used for events before 2002/250')
        end
        
        position_data = zeros(size(raw_data(i).data,1),1);
        position_data(snare_motion) = raw_data(i).data(snare_motion,snare_col);
        position_data(carri_motion) = raw_data(i).data(carri_motion,carri_col);
        position_data(latch_motion) = raw_data(i).data(latch_motion,latch_col);
        position_data(umbil_motion) = raw_data(i).data(umbil_motion,latch_col);
        position_data(snare_motion+carri_motion+latch_motion+umbil_motion > 1) = 0; % Remove data if multiple mechanisms in motion.
        
        mech_data = snare_motion+2*carri_motion+3*latch_motion+3*umbil_motion;
        mech_data(snare_motion+carri_motion+latch_motion+umbil_motion > 1) = 0;
        
    elseif time<2008050 % From 2002 day 250 until 2008 day 50
        
        cur_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMeasMotorCur'),1,'first');
        rate_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffDeriveMotRat'),1,'first');
        load_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffRigidForce'),1,'first');
        
        raw_data(i).data(raw_data(i).data(:,load_col)~=0,load_col) = ...
            raw_data(i).data(raw_data(i).data(:,load_col)~=0,load_col)  * 9.763653125 + 752.2539; % Convert to lbf
        
        pos_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffMechPos'),1,'first');
        mech_col = find(strcmpi(raw_data(i).textdata(1,2:end),'C_ActEffActiveMech'),1,'first');
        
        % Remove data points that have mechanism position > 30 but indicate
        % carriage is active.
        bool_rm_pts = (raw_data(i).data(:,mech_col)==1 & raw_data(i).data(:,pos_col)>30);
        
        raw_data(i).data = raw_data(i).data(~bool_rm_pts,:);
        raw_data(i).textdata = raw_data(i).textdata([true;~bool_rm_pts],:);
        processed_data.filterdata{i} = processed_data.filterdata{i}(~bool_rm_pts,:);
        
        if length([pos_col,mech_col]) ~= 2
            error('Error with position or mechanism PUI')
        end
        
        position_data = raw_data(i).data(:,pos_col);
        mech_data = raw_data(i).data(:,mech_col) + 1;
        
    elseif time < 2014250 % From 2008 day 50 until 2014 day 250
        
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
