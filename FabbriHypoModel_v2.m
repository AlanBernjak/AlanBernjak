% function runs the Fabbri model using randomly sampled input parameters Ko,gKr and IVshift
% Latin hypercube sampling with 30 randomly sampled input parameters is used

% Outputs include AP waveform time series, all input parameters and
% features of the AP waveform and are saved in the 'testFabbri' structure

% 'plot_Fabbri2' plots the CL values in the Ko,gKr,IVshift parameters space

% uses 'findfiducial4.m' script to calculate the features of the AP
% waveform

% sympathetic and parasympathetic activities can be included separately by 
% using the 'vagal' and 'symp' variables

clear 
close all

tic

nParams = 3;
nSampls = 30;

[X, parameters, FabbriOutput, biomarkers] = generateDesignData(nParams, nSampls);
% [parameters, FabbriOutput, biomarkers] = generateDesignData2('Data3d_2000_60_2');

t=toc

save testFabbri2 X parameters FabbriOutput biomarkers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, parameters, FabbriOutput, biomarkers] = generateDesignData(nParams, nSampls)

    startval.Ko = 4.0;
    startval.x_g_Kr = 1;
    startval.IHerg_shift_pa = 0;
    startval.IHerg_shift_pi = 0; % set equal to above later in the initconst function
        
    startval.ACh = 0;
    startval.shift_Ify = 0;
    startval.P_CaL = 0.4578;    
    startval.Pup_basal = 5;
    
    startval.Iso = 0;   
    startval.gKs = 0.00065;
    startval.Isoshift_Ify = 0;
    startval.Iso_i_NaK = 1;
    startval.Isoshift_dL = 0;
    startval.Isoshift_slopedL = 0;
    startval.Isoshift_ninf = 0;
    
%     startval.tau_dL = ;
%     startval.tau_n = ;
    
%     vagal = 'TRUE1'; % switching vagal stimulus ON
%     vagal = 'TRUE2';
    vagal = 'FALSE';   % switching vagal stimulus OFF
    
%     symp = 'TRUE1';  % switching sympathetic stimulus ON
%     symp = 'TRUE2';
    symp = 'FALSE';    % switching sympathetic stimulus OFF
    
    %%%%%
    
    pFields = {'Ko','x_g_Kr','IHerg_shift_pa'};
    
    X = lhsdesign(nSampls,nParams,'iterations',1000); % Latin hypercube sampling
%     X = [X X];
    for sample = 1 : nSampls
        
        disp(sprintf('design data sample %d of %d',sample, nSampls));
                
        parameters(sample).ACh = startval.ACh;
        parameters(sample).shift_Ify = startval.shift_Ify;
        parameters(sample).Iso = startval.Iso;
        parameters(sample).Isoshift_Ify =startval.Isoshift_Ify;
        parameters(sample).Iso_i_NaK = startval.Iso_i_NaK;
        parameters(sample).Isoshift_dL = startval.Isoshift_dL;
        parameters(sample).Isoshift_slopedL = startval.Isoshift_slopedL;
        parameters(sample).Isoshift_ninf = startval.Isoshift_ninf;
        parameters(sample).P_CaL = startval.P_CaL;
        parameters(sample).Pup_basal = startval.Pup_basal;   
        
        parameters(sample).IHerg_shift_pa = startval.IHerg_shift_pa;
        parameters(sample).IHerg_shift_pi = startval.IHerg_shift_pi;
        parameters(sample).gKs = startval.gKs;
    
        Xrow = X(sample,:);

%       % predefined (non-random) input parameters 
%         glu = 1:-0.2:0;%[1 0.8 0.6 0.4 0.2 0];
%         xgKr = 0.7+glu*0.3;
%         hergshift = 0:-0.5:-3; %0:-0.2:-5.4;% -5.4 %0:-0.2:-4.2;
%         xgKr = 0.0011*V+0.8372;
%         Kovec=0.6:0.1:6;
%         Kovec=[4 3.5 3];
        
        % changing Ko between 3 and 4
        for k = 1 %:nParams      
            
            % fixed Ko
%           % parameters(sample).(pFields{k}) = 4;  

%           % randomly sampled Ko between 3 and 4 
            parameters(sample).(pFields{k}) = 3.0 + Xrow(k) * (4.0-3.0);
%             
%           % default Ko
%           parameters(sample).(pFields{k}) = startval.Ko;

            
        end
        
        % changing x_g_Kr between 0.7 and 1
        for k = 2 %:nParams            
            
            % fixed gKr
            % parameters(sample).(pFields{k}) = 1; 
            
%           % randomly sampled gKr between 0.7 and 1
            parameters(sample).(pFields{k}) = 0.7 + Xrow(k)*(1-0.7); 
%             
%           % default gKr
%           parameters(sample).(pFields{k}) = startval.x_g_Kr;
           
        end
        
        % changing IVshift between 0 and -3 mV
        for k = 3%:nParams            
            
            % fixed IVshift
            % parameters(sample).(pFields{k}) = 0;
            
%           % randomly sampled IVshift between -3 and 0
            parameters(sample).(pFields{k}) = -3 + Xrow(k)*3.0; 
%
%           %default IVshift
%           parameters(sample).(pFields{k}) = startval.IHerg_shift_pa;

        end
               
        
        if strcmp(vagal,'TRUE1')
            parameters(sample).ACh = 1e-5; % in mM!
                        
        elseif strcmp(vagal,'TRUE2')
            parameters(sample).ACh = 2.5e-5; % in mM!
%             parameters(sample).ACh = 1.5e-5; % in mM! % brings HR down to
%             ~40bpm at Ko=4.0mmol/L
            
        elseif strcmp(symp,'TRUE1')
            
            parameters(sample).Iso = 1;
            parameters(sample).P_CaL = startval.P_CaL*(1+0.23);
            parameters(sample).Pup_basal = startval.Pup_basal*(1+0.25); 
            parameters(sample).gKs = startval.gKs*(1+0.20);
            parameters(sample).Isoshift_Ify = 7.5;
            parameters(sample).Iso_i_NaK = 1+0.2;
            parameters(sample).Isoshift_dL = -8;
            parameters(sample).Isoshift_slopedL = -27; % correct
            parameters(sample).Isoshift_ninf = -14;
            
        elseif strcmp(symp,'TRUE2')
            
            parameters(sample).Iso = 1;
            parameters(sample).P_CaL = startval.P_CaL*(1+0.968);
            parameters(sample).Pup_basal = startval.Pup_basal*(1+0.25); 
            parameters(sample).gKs = startval.gKs*(1+0.92);
            parameters(sample).Isoshift_Ify = 12;
            parameters(sample).Iso_i_NaK = 1+0.92;
            parameters(sample).Isoshift_dL = -12.8;
            parameters(sample).Isoshift_slopedL = -29.7; % correct
            parameters(sample).Isoshift_ninf = -22;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        parameters(sample).tspan_end = 60; % duration of simulation
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        output = RunFabbriA2017(parameters(sample));
        
        biomarkers(sample,:) = findfiducial4(output.time,output.volt,0,0);
        
%         FabbriOutput = 0;
        FabbriOutput(sample).time = output.time; % only saving AP waveform
        FabbriOutput(sample).volt = output.volt;
        
% saving currents to the output structre        
%         C=57; %C=57pF;
                
%         FabbriOutput(sample).If = output.ALGEBRAIC(:,60)*1000/C;
%         FabbriOutput(sample).INaCa = output.ALGEBRAIC(:,76)*1000/C;
%         FabbriOutput(sample).ICaL = output.ALGEBRAIC(:,84)*1000/C;
%         FabbriOutput(sample).ICaT = output.ALGEBRAIC(:,85)*1000/C;
%         FabbriOutput(sample).IKr = output.ALGEBRAIC(:,90)*1000/C;
%         FabbriOutput(sample).IKs = output.ALGEBRAIC(:,96)*1000/C;
%         FabbriOutput(sample).INaK = output.ALGEBRAIC(:,62)*1000/C;
%         FabbriOutput(sample).INa = output.ALGEBRAIC(:,80)*1000/C;
%         FabbriOutput(sample).Ito = output.ALGEBRAIC(:,88)*1000/C;
%         FabbriOutput(sample).IKACh = output.ALGEBRAIC(:,99)*1000/C;
%         FabbriOutput(sample).IKur = output.ALGEBRAIC(:,86)*1000/C;
        
%         FabbriOutput(sample).Cai = output.STATES(:,18);
%         FabbriOutput(sample).Casub = output.STATES(:,2);
%         FabbriOutput(sample).CajSR = output.STATES(:,16);
%         FabbriOutput(sample).CanSR = output.STATES(:,17);
    
        
        clear output
              
        
    end
end

% only for use with predefined input parameters
function [parameters, FabbriOutput, biomarkers] = generateDesignData2(filename_data)

    startval.Ko = 4.0;
    startval.x_g_Kr = 1;
    startval.IHerg_shift_pa = 0;
    startval.IHerg_shift_pi = 0;
        
    startval.ACh = 0;
    startval.shift_Ify = 0;
    startval.P_CaL = 0.4578;    
    startval.Pup_basal = 5;
    
    startval.Iso = 0;   
    startval.gKs = 0.00065;
    startval.Isoshift_Ify = 0;
    startval.Iso_i_NaK = 1;
    startval.Isoshift_dL = 0;
    startval.Isoshift_slopedL = 0;
    startval.Isoshift_ninf = 0;
    
%     startval.tau_dL = ;
%     startval.tau_n = ;
    
%     vagal = 'TRUE1';
%     vagal = 'TRUE2';
    vagal = 'FALSE';
    
%     symp = 'TRUE1';
%     symp = 'TRUE2';
    symp = 'FALSE';
    
    %%%%%
    
    pFields = {'Ko','x_g_Kr','IHerg_shift_pa'};
    
    data = load(filename_data);
    
    
    for sample = 1:length(data.parameters)
        
        disp(sprintf('design data sample %d',sample));
                
        parameters(sample).ACh = startval.ACh;
        parameters(sample).shift_Ify = startval.shift_Ify;
        parameters(sample).Iso = startval.Iso;
        parameters(sample).Isoshift_Ify =startval.Isoshift_Ify;
        parameters(sample).Iso_i_NaK = startval.Iso_i_NaK;
        parameters(sample).Isoshift_dL = startval.Isoshift_dL;
        parameters(sample).Isoshift_slopedL = startval.Isoshift_slopedL;
        parameters(sample).Isoshift_ninf = startval.Isoshift_ninf;
        parameters(sample).P_CaL = startval.P_CaL;
        parameters(sample).Pup_basal = startval.Pup_basal;   
        
%         parameters(sample).IHerg_shift_pa = parameters(sample).IHerg_shift_pa;
%         parameters(sample).IHerg_shift_pi = parameters(sample).IHerg_shift_pa;
        parameters(sample).gKs = startval.gKs;
    
        parameters(sample).Ko = data.parameters(sample).Ko;
        parameters(sample).x_g_Kr = data.parameters(sample).x_g_Kr; 
        parameters(sample).IHerg_shift_pa = data.parameters(sample).IHerg_shift_pa;
        

        if strcmp(vagal,'TRUE1')
            parameters(sample).ACh = 1e-5; % in mM!
                        
        elseif strcmp(vagal,'TRUE2')
            parameters(sample).ACh = 2.5e-5; % in mM!
%             parameters(sample).ACh = 1.5e-5; % in mM! % brings HR down to
%             ~40bpm at Ko=4.0mmol/L
            
        elseif strcmp(symp,'TRUE1')
            
            parameters(sample).Iso = 1;
            parameters(sample).P_CaL = startval.P_CaL*(1+0.23);
            parameters(sample).Pup_basal = startval.Pup_basal*(1+0.25); 
            parameters(sample).gKs = startval.gKs*(1+0.20);
            parameters(sample).Isoshift_Ify = 7.5;
            parameters(sample).Iso_i_NaK = 1+0.2;
            parameters(sample).Isoshift_dL = -8;
            parameters(sample).Isoshift_slopedL = -27; % correct
            parameters(sample).Isoshift_ninf = -14;
            
        elseif strcmp(symp,'TRUE2')
            
            parameters(sample).Iso = 1;
            parameters(sample).P_CaL = startval.P_CaL*(1+0.968);
            parameters(sample).Pup_basal = startval.Pup_basal*(1+0.25); 
            parameters(sample).gKs = startval.gKs*(1+0.92);
            parameters(sample).Isoshift_Ify = 12;
            parameters(sample).Iso_i_NaK = 1+0.92;
            parameters(sample).Isoshift_dL = -12.8;
            parameters(sample).Isoshift_slopedL = -29.7; % correct
            parameters(sample).Isoshift_ninf = -22;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        parameters(sample).tspan_end = 60;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        output = RunFabbriA2017(parameters(sample));
        
        biomarkers(sample,:) = findfiducial4(output.time,output.volt,0,0);
        
%         FabbriOutput(sample).time = output.time;
%         FabbriOutput(sample).volt = output.volt;
        
        FabbriOutput = 0;
        
%         plot_concs = 1;
%         if plot_concs == 1
%              plotconcentrations
%         end         
        
        clear output
              
        
    end
end

function output = RunFabbriA2017(parameters)
    
    [VOI, STATES, ALGEBRAIC, CONSTANTS] = FabbriA2017CellML(parameters);
    
    output.time = VOI;
    output.volt = STATES(:,1);
    output.STATES = STATES;
    output.CONSTANTS = CONSTANTS;
    output.ALGEBRAIC = ALGEBRAIC;
    
end

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = FabbriA2017CellML(parameters)
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.   
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel(parameters);
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =101;
end
% There are a total of 33 entries in each of the rate and state variable arrays.
% There are a total of 116 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel(parameters)
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts(parameters);
    
    tspan_end = parameters.tspan_end;

    % Set timespan to solve over 
    tspan = [0, tspan_end];

    % Set numerical accuracy options for ODE solver
% %     options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);
    options = odeset('RelTol', 1e-07, 'AbsTol', 1e-06, 'MaxStep', 0.001);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
%     figure();    
%     plot(VOI, STATES(:,1));
    
%     plot(VOI, STATES);
%     xlabel(LEGEND_VOI);
%     l = legend(LEGEND_STATES);
%     set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (second)');
    LEGEND_CONSTANTS(:,1) = strpad('R in component Membrane (joule_per_kilomole_kelvin)');
    LEGEND_CONSTANTS(:,2) = strpad('T in component Membrane (kelvin)');
    LEGEND_CONSTANTS(:,3) = strpad('F in component Membrane (coulomb_per_mole)');
    LEGEND_CONSTANTS(:,4) = strpad('C in component Membrane (microF)');
    LEGEND_CONSTANTS(:,92) = strpad('RTONF in component Membrane (millivolt)');
    LEGEND_ALGEBRAIC(:,60) = strpad('i_f in component i_f (nanoA)');
    LEGEND_ALGEBRAIC(:,62) = strpad('i_NaK in component i_NaK (nanoA)');
    LEGEND_ALGEBRAIC(:,76) = strpad('i_NaCa in component i_NaCa (nanoA)');
    LEGEND_ALGEBRAIC(:,80) = strpad('i_Na in component i_Na (nanoA)');
    LEGEND_ALGEBRAIC(:,90) = strpad('i_Kr in component i_Kr (nanoA)');
    LEGEND_ALGEBRAIC(:,96) = strpad('i_Ks in component i_Ks (nanoA)');
    LEGEND_ALGEBRAIC(:,88) = strpad('i_to in component i_to (nanoA)');
    LEGEND_ALGEBRAIC(:,84) = strpad('i_CaL in component i_CaL (nanoA)');
    LEGEND_ALGEBRAIC(:,85) = strpad('i_CaT in component i_CaT (nanoA)');
    LEGEND_ALGEBRAIC(:,99) = strpad('i_KACh in component i_KACh (nanoA)');
    LEGEND_ALGEBRAIC(:,86) = strpad('i_Kur in component i_Kur (nanoA)');
    LEGEND_ALGEBRAIC(:,10) = strpad('V in component Membrane (millivolt)');
    LEGEND_CONSTANTS(:,5) = strpad('clamp_mode in component Membrane (dimensionless)');
    LEGEND_ALGEBRAIC(:,6) = strpad('V_clamp in component Voltage_clamp (millivolt)');
    LEGEND_STATES(:,1) = strpad('V_ode in component Membrane (millivolt)');
    LEGEND_ALGEBRAIC(:,101) = strpad('i_tot in component Membrane (nanoA)');
    LEGEND_CONSTANTS(:,6) = strpad('t_holding in component Voltage_clamp (second)');
    LEGEND_CONSTANTS(:,7) = strpad('t_test in component Voltage_clamp (second)');
    LEGEND_CONSTANTS(:,8) = strpad('V_test in component Voltage_clamp (millivolt)');
    LEGEND_CONSTANTS(:,9) = strpad('V_holding in component Voltage_clamp (millivolt)');
    LEGEND_CONSTANTS(:,10) = strpad('ACh in component Rate_modulation_experiments (millimolar)');
    LEGEND_CONSTANTS(:,11) = strpad('Iso_1_uM in component Rate_modulation_experiments (dimensionless)');
    LEGEND_ALGEBRAIC(:,19) = strpad('Nai in component Nai_concentration (millimolar)');
    LEGEND_CONSTANTS(:,12) = strpad('Nao in component Ionic_values (millimolar)');
    LEGEND_CONSTANTS(:,13) = strpad('Ki in component Ionic_values (millimolar)');
    LEGEND_CONSTANTS(:,14) = strpad('Ko in component Ionic_values (millimolar)');
    LEGEND_STATES(:,2) = strpad('Ca_sub in component Ca_dynamics (millimolar)');
    LEGEND_CONSTANTS(:,15) = strpad('Cao in component Ionic_values (millimolar)');
    LEGEND_ALGEBRAIC(:,37) = strpad('E_Na in component Ionic_values (millivolt)');
    LEGEND_CONSTANTS(:,97) = strpad('E_K in component Ionic_values (millivolt)');
    LEGEND_ALGEBRAIC(:,1) = strpad('E_Ca in component Ionic_values (millivolt)');
    LEGEND_CONSTANTS(:,111) = strpad('V_sub in component Cell_parameters (millimetre3)');
    LEGEND_CONSTANTS(:,113) = strpad('V_i in component Cell_parameters (millimetre3)');
    LEGEND_ALGEBRAIC(:,50) = strpad('i_fNa in component i_f (nanoA)');
    LEGEND_ALGEBRAIC(:,83) = strpad('i_siNa in component i_CaL (nanoA)');
    LEGEND_STATES(:,3) = strpad('Nai_ in component Nai_concentration (millimolar)');
    LEGEND_CONSTANTS(:,16) = strpad('Nai_clamp in component Nai_concentration (dimensionless)');
    LEGEND_ALGEBRAIC(:,56) = strpad('i_fK in component i_f (nanoA)');
    LEGEND_CONSTANTS(:,17) = strpad('g_f in component i_f (microS)');
    LEGEND_CONSTANTS(:,93) = strpad('G_f in component i_f (microS)');
    LEGEND_CONSTANTS(:,104) = strpad('g_f_Na in component i_f (microS)');
    LEGEND_CONSTANTS(:,101) = strpad('G_f_Na in component i_f (microS)');
    LEGEND_CONSTANTS(:,102) = strpad('g_f_K in component i_f (microS)');
    LEGEND_CONSTANTS(:,98) = strpad('G_f_K in component i_f (microS)');
    LEGEND_CONSTANTS(:,18) = strpad('Km_f in component i_f (millimolar)');
    LEGEND_CONSTANTS(:,19) = strpad('alpha in component i_f (dimensionless)');
    LEGEND_STATES(:,4) = strpad('y in component i_f_y_gate (dimensionless)');
    LEGEND_CONSTANTS(:,20) = strpad('blockade in component i_f (dimensionless)');
    LEGEND_ALGEBRAIC(:,11) = strpad('tau_y in component i_f_y_gate (second)');
    LEGEND_ALGEBRAIC(:,30) = strpad('y_infinity in component i_f_y_gate (dimensionless)');
    LEGEND_CONSTANTS(:,96) = strpad('ACh_shift in component i_f_y_gate (millivolt)');
    LEGEND_CONSTANTS(:,100) = strpad('Iso_shift in component i_f_y_gate (millivolt)');
    LEGEND_CONSTANTS(:,21) = strpad('y_shift in component i_f_y_gate (millivolt)');
    LEGEND_CONSTANTS(:,22) = strpad('Km_Kp in component i_NaK (millimolar)');
    LEGEND_CONSTANTS(:,23) = strpad('Km_Nap in component i_NaK (millimolar)');
    LEGEND_CONSTANTS(:,24) = strpad('i_NaK_max in component i_NaK (nanoA)');
    LEGEND_CONSTANTS(:,103) = strpad('Iso_increase in component i_NaK (dimensionless)');
    LEGEND_CONSTANTS(:,25) = strpad('K_NaCa in component i_NaCa (nanoA)');
    LEGEND_ALGEBRAIC(:,73) = strpad('x1 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,69) = strpad('x2 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,74) = strpad('x3 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,75) = strpad('x4 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,64) = strpad('k41 in component i_NaCa (dimensionless)');
    LEGEND_CONSTANTS(:,105) = strpad('k34 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,72) = strpad('k23 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,71) = strpad('k21 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,68) = strpad('k32 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,63) = strpad('k43 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,66) = strpad('k12 in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,67) = strpad('k14 in component i_NaCa (dimensionless)');
    LEGEND_CONSTANTS(:,26) = strpad('Qci in component i_NaCa (dimensionless)');
    LEGEND_CONSTANTS(:,27) = strpad('Qn in component i_NaCa (dimensionless)');
    LEGEND_CONSTANTS(:,28) = strpad('Qco in component i_NaCa (dimensionless)');
    LEGEND_CONSTANTS(:,29) = strpad('K3ni in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,30) = strpad('Kci in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,31) = strpad('K1ni in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,32) = strpad('K2ni in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,33) = strpad('Kcni in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,34) = strpad('K3no in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,35) = strpad('K1no in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,36) = strpad('K2no in component i_NaCa (millimolar)');
    LEGEND_CONSTANTS(:,37) = strpad('Kco in component i_NaCa (millimolar)');
    LEGEND_ALGEBRAIC(:,70) = strpad('do in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,65) = strpad('di in component i_NaCa (dimensionless)');
    LEGEND_CONSTANTS(:,38) = strpad('blockade_NaCa in component i_NaCa (dimensionless)');
    LEGEND_ALGEBRAIC(:,78) = strpad('i_Na_ in component i_Na (nanoA)');
    LEGEND_ALGEBRAIC(:,79) = strpad('i_Na_L in component i_Na (nanoA)');
    LEGEND_CONSTANTS(:,39) = strpad('g_Na in component i_Na (microS)');
    LEGEND_CONSTANTS(:,40) = strpad('g_Na_L in component i_Na (microS)');
    LEGEND_ALGEBRAIC(:,77) = strpad('E_mh in component i_Na (millivolt)');
    LEGEND_STATES(:,5) = strpad('m in component i_Na_m_gate (dimensionless)');
    LEGEND_STATES(:,6) = strpad('h in component i_Na_h_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,47) = strpad('alpha_m in component i_Na_m_gate (per_second)');
    LEGEND_ALGEBRAIC(:,53) = strpad('beta_m in component i_Na_m_gate (per_second)');
    LEGEND_ALGEBRAIC(:,12) = strpad('m_infinity in component i_Na_m_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,58) = strpad('tau_m in component i_Na_m_gate (second)');
    LEGEND_CONSTANTS(:,41) = strpad('delta_m in component i_Na_m_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,31) = strpad('E0_m in component i_Na_m_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,32) = strpad('alpha_h in component i_Na_h_gate (per_second)');
    LEGEND_ALGEBRAIC(:,48) = strpad('beta_h in component i_Na_h_gate (per_second)');
    LEGEND_ALGEBRAIC(:,13) = strpad('h_infinity in component i_Na_h_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,54) = strpad('tau_h in component i_Na_h_gate (second)');
    LEGEND_ALGEBRAIC(:,81) = strpad('i_siCa in component i_CaL (nanoA)');
    LEGEND_ALGEBRAIC(:,82) = strpad('i_siK in component i_CaL (nanoA)');
    LEGEND_CONSTANTS(:,107) = strpad('ACh_block in component i_CaL (dimensionless)');
    LEGEND_CONSTANTS(:,42) = strpad('P_CaL in component i_CaL (nanoA_per_millimolar)');
    LEGEND_STATES(:,7) = strpad('dL in component i_CaL_dL_gate (dimensionless)');
    LEGEND_STATES(:,8) = strpad('fL in component i_CaL_fL_gate (dimensionless)');
    LEGEND_STATES(:,9) = strpad('fCa in component i_CaL_fCa_gate (dimensionless)');
    LEGEND_CONSTANTS(:,106) = strpad('Iso_increase in component i_CaL (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('dL_infinity in component i_CaL_dL_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,61) = strpad('tau_dL in component i_CaL_dL_gate (second)');
    LEGEND_ALGEBRAIC(:,49) = strpad('alpha_dL in component i_CaL_dL_gate (per_second)');
    LEGEND_ALGEBRAIC(:,59) = strpad('beta_dL in component i_CaL_dL_gate (per_second)');
    LEGEND_ALGEBRAIC(:,33) = strpad('adVm in component i_CaL_dL_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,55) = strpad('bdVm in component i_CaL_dL_gate (millivolt)');
    LEGEND_CONSTANTS(:,43) = strpad('k_dL in component i_CaL_dL_gate (millivolt)');
    LEGEND_CONSTANTS(:,44) = strpad('V_dL in component i_CaL_dL_gate (millivolt)');
    LEGEND_CONSTANTS(:,108) = strpad('Iso_shift_dL in component i_CaL_dL_gate (millivolt)');
    LEGEND_CONSTANTS(:,109) = strpad('Iso_slope_dL in component i_CaL_dL_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,15) = strpad('fL_infinity in component i_CaL_fL_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,34) = strpad('tau_fL in component i_CaL_fL_gate (second)');
    LEGEND_CONSTANTS(:,45) = strpad('shift_fL in component i_CaL_fL_gate (millivolt)');
    LEGEND_CONSTANTS(:,46) = strpad('k_fL in component i_CaL_fL_gate (millivolt)');
    LEGEND_CONSTANTS(:,47) = strpad('alpha_fCa in component i_CaL_fCa_gate (per_second)');
    LEGEND_ALGEBRAIC(:,2) = strpad('fCa_infinity in component i_CaL_fCa_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,8) = strpad('tau_fCa in component i_CaL_fCa_gate (second)');
    LEGEND_CONSTANTS(:,48) = strpad('Km_fCa in component i_CaL_fCa_gate (millimolar)');
    LEGEND_CONSTANTS(:,49) = strpad('P_CaT in component i_CaT (nanoA_per_millimolar)');
    LEGEND_STATES(:,10) = strpad('dT in component i_CaT_dT_gate (dimensionless)');
    LEGEND_STATES(:,11) = strpad('fT in component i_CaT_fT_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,16) = strpad('dT_infinity in component i_CaT_dT_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,35) = strpad('tau_dT in component i_CaT_dT_gate (second)');
    LEGEND_ALGEBRAIC(:,17) = strpad('fT_infinity in component i_CaT_fT_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,36) = strpad('tau_fT in component i_CaT_fT_gate (second)');
    LEGEND_CONSTANTS(:,50) = strpad('offset_fT in component i_CaT_fT_gate (second)');
    LEGEND_ALGEBRAIC(:,87) = strpad('j_SRCarel in component Ca_SR_release (millimolar_per_second)');
    LEGEND_STATES(:,12) = strpad('R in component Ca_SR_release (dimensionless)');
    LEGEND_STATES(:,13) = strpad('O in component Ca_SR_release (dimensionless)');
    LEGEND_STATES(:,14) = strpad('I in component Ca_SR_release (dimensionless)');
    LEGEND_STATES(:,15) = strpad('RI in component Ca_SR_release (dimensionless)');
    LEGEND_CONSTANTS(:,51) = strpad('ks in component Ca_SR_release (per_second)');
    LEGEND_CONSTANTS(:,52) = strpad('MaxSR in component Ca_SR_release (dimensionless)');
    LEGEND_CONSTANTS(:,53) = strpad('MinSR in component Ca_SR_release (dimensionless)');
    LEGEND_CONSTANTS(:,54) = strpad('EC50_SR in component Ca_SR_release (millimolar)');
    LEGEND_CONSTANTS(:,55) = strpad('HSR in component Ca_SR_release (dimensionless)');
    LEGEND_ALGEBRAIC(:,9) = strpad('koSRCa in component Ca_SR_release (per_millimolar2_second)');
    LEGEND_ALGEBRAIC(:,18) = strpad('kiSRCa in component Ca_SR_release (per_millimolar_second)');
    LEGEND_CONSTANTS(:,56) = strpad('koCa in component Ca_SR_release (per_millimolar2_second)');
    LEGEND_CONSTANTS(:,57) = strpad('kiCa in component Ca_SR_release (per_millimolar_second)');
    LEGEND_ALGEBRAIC(:,3) = strpad('kCaSR in component Ca_SR_release (dimensionless)');
    LEGEND_CONSTANTS(:,58) = strpad('kim in component Ca_SR_release (per_second)');
    LEGEND_CONSTANTS(:,59) = strpad('kom in component Ca_SR_release (per_second)');
    LEGEND_STATES(:,16) = strpad('Ca_jsr in component Ca_dynamics (millimolar)');
    LEGEND_ALGEBRAIC(:,4) = strpad('diff in component Ca_SR_release (millimolar)');
    LEGEND_ALGEBRAIC(:,5) = strpad('P_tot in component Ca_SR_release (dimensionless)');
    LEGEND_ALGEBRAIC(:,89) = strpad('j_Ca_dif in component Ca_intracellular_fluxes (millimolar_per_second)');
    LEGEND_ALGEBRAIC(:,92) = strpad('j_up in component Ca_intracellular_fluxes (millimolar_per_second)');
    LEGEND_ALGEBRAIC(:,95) = strpad('j_tr in component Ca_intracellular_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,60) = strpad('tau_dif_Ca in component Ca_intracellular_fluxes (second)');
    LEGEND_CONSTANTS(:,61) = strpad('tau_tr in component Ca_intracellular_fluxes (second)');
    LEGEND_CONSTANTS(:,99) = strpad('P_up in component Ca_intracellular_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,62) = strpad('P_up_basal in component Ca_intracellular_fluxes (millimolar_per_second)');
    LEGEND_CONSTANTS(:,94) = strpad('b_up in component Ca_intracellular_fluxes (dimensionless)');
    LEGEND_CONSTANTS(:,63) = strpad('K_up in component Ca_intracellular_fluxes (millimolar)');
    LEGEND_STATES(:,17) = strpad('Ca_nsr in component Ca_dynamics (millimolar)');
    LEGEND_STATES(:,18) = strpad('Cai in component Ca_dynamics (millimolar)');
    LEGEND_CONSTANTS(:,64) = strpad('slope_up in component Ca_intracellular_fluxes (millimolar)');
    LEGEND_CONSTANTS(:,65) = strpad('TC_tot in component Ca_buffering (millimolar)');
    LEGEND_CONSTANTS(:,66) = strpad('TMC_tot in component Ca_buffering (millimolar)');
    LEGEND_CONSTANTS(:,67) = strpad('CM_tot in component Ca_buffering (millimolar)');
    LEGEND_CONSTANTS(:,68) = strpad('CQ_tot in component Ca_buffering (millimolar)');
    LEGEND_ALGEBRAIC(:,94) = strpad('delta_fTC in component Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,97) = strpad('delta_fTMC in component Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,91) = strpad('delta_fCMs in component Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,100) = strpad('delta_fCMi in component Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,98) = strpad('delta_fCQ in component Ca_buffering (per_second)');
    LEGEND_ALGEBRAIC(:,7) = strpad('delta_fTMM in component Ca_buffering (per_second)');
    LEGEND_STATES(:,19) = strpad('fTMM in component Ca_buffering (dimensionless)');
    LEGEND_STATES(:,20) = strpad('fCMi in component Ca_buffering (dimensionless)');
    LEGEND_STATES(:,21) = strpad('fCMs in component Ca_buffering (dimensionless)');
    LEGEND_STATES(:,22) = strpad('fTC in component Ca_buffering (dimensionless)');
    LEGEND_STATES(:,23) = strpad('fTMC in component Ca_buffering (dimensionless)');
    LEGEND_STATES(:,24) = strpad('fCQ in component Ca_buffering (dimensionless)');
    LEGEND_CONSTANTS(:,69) = strpad('kf_TC in component Ca_buffering (per_millimolar_second)');
    LEGEND_CONSTANTS(:,70) = strpad('kf_TMM in component Ca_buffering (per_millimolar_second)');
    LEGEND_CONSTANTS(:,71) = strpad('kf_TMC in component Ca_buffering (per_millimolar_second)');
    LEGEND_CONSTANTS(:,72) = strpad('kf_CM in component Ca_buffering (per_millimolar_second)');
    LEGEND_CONSTANTS(:,73) = strpad('kf_CQ in component Ca_buffering (per_millimolar_second)');
    LEGEND_CONSTANTS(:,74) = strpad('kb_TC in component Ca_buffering (per_second)');
    LEGEND_CONSTANTS(:,75) = strpad('kb_TMC in component Ca_buffering (per_second)');
    LEGEND_CONSTANTS(:,76) = strpad('kb_TMM in component Ca_buffering (per_second)');
    LEGEND_CONSTANTS(:,77) = strpad('kb_CM in component Ca_buffering (per_second)');
    LEGEND_CONSTANTS(:,78) = strpad('kb_CQ in component Ca_buffering (per_second)');
    LEGEND_CONSTANTS(:,79) = strpad('Mgi in component Ca_buffering (millimolar)');
    LEGEND_CONSTANTS(:,112) = strpad('V_jsr in component Cell_parameters (millimetre3)');
    LEGEND_CONSTANTS(:,114) = strpad('V_nsr in component Cell_parameters (millimetre3)');
    LEGEND_CONSTANTS(:,110) = strpad('V_cell in component Cell_parameters (millimetre3)');
    LEGEND_CONSTANTS(:,80) = strpad('V_jsr_part in component Cell_parameters (dimensionless)');
    LEGEND_CONSTANTS(:,81) = strpad('V_i_part in component Cell_parameters (dimensionless)');
    LEGEND_CONSTANTS(:,82) = strpad('V_nsr_part in component Cell_parameters (dimensionless)');
    LEGEND_CONSTANTS(:,83) = strpad('R_cell in component Cell_parameters (micrometre)');
    LEGEND_CONSTANTS(:,84) = strpad('L_cell in component Cell_parameters (micrometre)');
    LEGEND_CONSTANTS(:,85) = strpad('L_sub in component Cell_parameters (micrometre)');
    LEGEND_CONSTANTS(:,86) = strpad('g_Kur in component i_Kur (microS)');
    LEGEND_STATES(:,25) = strpad('r_Kur in component i_Kur_rKur_gate (dimensionless)');
    LEGEND_STATES(:,26) = strpad('s_Kur in component i_Kur_sKur_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,38) = strpad('tau_r_Kur in component i_Kur_rKur_gate (second)');
    LEGEND_ALGEBRAIC(:,20) = strpad('r_Kur_infinity in component i_Kur_rKur_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,39) = strpad('tau_s_Kur in component i_Kur_sKur_gate (second)');
    LEGEND_ALGEBRAIC(:,21) = strpad('s_Kur_infinity in component i_Kur_sKur_gate (dimensionless)');
    LEGEND_CONSTANTS(:,87) = strpad('g_to in component i_to (microS)');
    LEGEND_STATES(:,27) = strpad('q in component i_to_q_gate (dimensionless)');
    LEGEND_STATES(:,28) = strpad('r in component i_to_r_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,22) = strpad('q_infinity in component i_to_q_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,40) = strpad('tau_q in component i_to_q_gate (second)');
    LEGEND_ALGEBRAIC(:,23) = strpad('r_infinity in component i_to_r_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,41) = strpad('tau_r in component i_to_r_gate (second)');
    LEGEND_CONSTANTS(:,88) = strpad('g_Kr in component i_Kr (microS)');
    %%% added mulitplication factor for gKr
    LEGEND_CONSTANTS(:,117) = strpad('x_g_Kr: multiplication factor for g_Kr in component i_Kr');
    LEGEND_CONSTANTS(:,118) = strpad('Voltage shift in component i_Kr pa gate due to lowered [glucose]o');
    LEGEND_CONSTANTS(:,133) = strpad('Voltage shift in component i_Kr pi gate due to lowered [glucose]o');
    LEGEND_STATES(:,29) = strpad('paS in component i_Kr_pa_gate (dimensionless)');
    LEGEND_STATES(:,30) = strpad('paF in component i_Kr_pa_gate (dimensionless)');
    LEGEND_STATES(:,31) = strpad('piy in component i_Kr_pi_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,24) = strpad('pa_infinity in component i_Kr_pa_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('alfapaF in component i_Kr_pa_gate (per_second)');
    LEGEND_ALGEBRAIC(:,26) = strpad('betapaF in component i_Kr_pa_gate (per_second)');
    LEGEND_ALGEBRAIC(:,42) = strpad('tau_paS in component i_Kr_pa_gate (second)');
    LEGEND_ALGEBRAIC(:,43) = strpad('tau_paF in component i_Kr_pa_gate (second)');
    LEGEND_ALGEBRAIC(:,44) = strpad('pi_infinity in component i_Kr_pi_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,27) = strpad('tau_pi in component i_Kr_pi_gate (second)');
    LEGEND_CONSTANTS(:,95) = strpad('g_Ks in component i_Ks (microS)');
    LEGEND_CONSTANTS(:,89) = strpad('g_Ks_ in component i_Ks (microS)');
    LEGEND_ALGEBRAIC(:,93) = strpad('E_Ks in component i_Ks (millivolt)');
    LEGEND_STATES(:,32) = strpad('n in component i_Ks_n_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,28) = strpad('n_infinity in component i_Ks_n_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,57) = strpad('tau_n in component i_Ks_n_gate (second)');
    LEGEND_CONSTANTS(:,115) = strpad('Iso_shift in component i_Ks_n_gate (millivolt)');
    LEGEND_ALGEBRAIC(:,45) = strpad('alpha_n in component i_Ks_n_gate (per_second)');
    LEGEND_ALGEBRAIC(:,51) = strpad('beta_n in component i_Ks_n_gate (per_second)');
    LEGEND_CONSTANTS(:,90) = strpad('g_KACh in component i_KACh (microS)');
    LEGEND_STATES(:,33) = strpad('a in component i_KACh_a_gate (dimensionless)');
    LEGEND_CONSTANTS(:,91) = strpad('ACh_on in component i_KACh (dimensionless)');
    LEGEND_CONSTANTS(:,116) = strpad('alpha_a in component i_KACh_a_gate (per_second)');
    LEGEND_ALGEBRAIC(:,29) = strpad('beta_a in component i_KACh_a_gate (per_second)');
    LEGEND_ALGEBRAIC(:,46) = strpad('a_infinity in component i_KACh_a_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,52) = strpad('tau_a in component i_KACh_a_gate (second)');
    LEGEND_RATES(:,1) = strpad('d/dt V_ode in component Membrane (millivolt)');
    LEGEND_RATES(:,3) = strpad('d/dt Nai_ in component Nai_concentration (millimolar)');
    LEGEND_RATES(:,4) = strpad('d/dt y in component i_f_y_gate (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt m in component i_Na_m_gate (dimensionless)');
    LEGEND_RATES(:,6) = strpad('d/dt h in component i_Na_h_gate (dimensionless)');
    LEGEND_RATES(:,7) = strpad('d/dt dL in component i_CaL_dL_gate (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt fL in component i_CaL_fL_gate (dimensionless)');
    LEGEND_RATES(:,9) = strpad('d/dt fCa in component i_CaL_fCa_gate (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt dT in component i_CaT_dT_gate (dimensionless)');
    LEGEND_RATES(:,11) = strpad('d/dt fT in component i_CaT_fT_gate (dimensionless)');
    LEGEND_RATES(:,12) = strpad('d/dt R in component Ca_SR_release (dimensionless)');
    LEGEND_RATES(:,13) = strpad('d/dt O in component Ca_SR_release (dimensionless)');
    LEGEND_RATES(:,14) = strpad('d/dt I in component Ca_SR_release (dimensionless)');
    LEGEND_RATES(:,15) = strpad('d/dt RI in component Ca_SR_release (dimensionless)');
    LEGEND_RATES(:,22) = strpad('d/dt fTC in component Ca_buffering (dimensionless)');
    LEGEND_RATES(:,23) = strpad('d/dt fTMC in component Ca_buffering (dimensionless)');
    LEGEND_RATES(:,19) = strpad('d/dt fTMM in component Ca_buffering (dimensionless)');
    LEGEND_RATES(:,20) = strpad('d/dt fCMi in component Ca_buffering (dimensionless)');
    LEGEND_RATES(:,21) = strpad('d/dt fCMs in component Ca_buffering (dimensionless)');
    LEGEND_RATES(:,24) = strpad('d/dt fCQ in component Ca_buffering (dimensionless)');
    LEGEND_RATES(:,18) = strpad('d/dt Cai in component Ca_dynamics (millimolar)');
    LEGEND_RATES(:,2) = strpad('d/dt Ca_sub in component Ca_dynamics (millimolar)');
    LEGEND_RATES(:,17) = strpad('d/dt Ca_nsr in component Ca_dynamics (millimolar)');
    LEGEND_RATES(:,16) = strpad('d/dt Ca_jsr in component Ca_dynamics (millimolar)');
    LEGEND_RATES(:,25) = strpad('d/dt r_Kur in component i_Kur_rKur_gate (dimensionless)');
    LEGEND_RATES(:,26) = strpad('d/dt s_Kur in component i_Kur_sKur_gate (dimensionless)');
    LEGEND_RATES(:,27) = strpad('d/dt q in component i_to_q_gate (dimensionless)');
    LEGEND_RATES(:,28) = strpad('d/dt r in component i_to_r_gate (dimensionless)');
    LEGEND_RATES(:,29) = strpad('d/dt paS in component i_Kr_pa_gate (dimensionless)');
    LEGEND_RATES(:,30) = strpad('d/dt paF in component i_Kr_pa_gate (dimensionless)');
    LEGEND_RATES(:,31) = strpad('d/dt piy in component i_Kr_pi_gate (dimensionless)');
    LEGEND_RATES(:,32) = strpad('d/dt n in component i_Ks_n_gate (dimensionless)');
    LEGEND_RATES(:,33) = strpad('d/dt a in component i_KACh_a_gate (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts(parameters)

    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];

    Ko = parameters.Ko;
    x_g_Kr = parameters.x_g_Kr;
    IHerg_shift_pa = parameters.IHerg_shift_pa;
%     IHerg_shift_pi = parameters.IHerg_shift_pi;
    IHerg_shift_pi = IHerg_shift_pa;
    
    ACh = parameters.ACh;
    shift_Ify = parameters.shift_Ify;
    Pup_basal = parameters.Pup_basal;
    P_CaL = parameters.P_CaL;
    Iso = parameters.Iso;
    Isoshift_Ify = parameters.Isoshift_Ify;
    gKs = parameters.gKs;
    Iso_i_NaK = parameters.Iso_i_NaK;
    Isoshift_dL = parameters.Isoshift_dL;
    Isoshift_slopedL = parameters.Isoshift_slopedL;
    Isoshift_ninf = parameters.Isoshift_ninf;
    
    CONSTANTS(:,1) = 8314.472;
    CONSTANTS(:,2) = 310;
    CONSTANTS(:,3) = 96485.3415;
    CONSTANTS(:,4) = 5.7e-5;
    CONSTANTS(:,5) = 0;
    STATES(:,1) = -47.787168;
    CONSTANTS(:,6) = 0.5;
    CONSTANTS(:,7) = 0.5;
    CONSTANTS(:,8) = -35;
    CONSTANTS(:,9) = -45;
    CONSTANTS(:,10) = ACh;
    CONSTANTS(:,11) = Iso;
    CONSTANTS(:,12) = 140;
    CONSTANTS(:,13) = 140;
    CONSTANTS(:,14) = Ko;
    STATES(:,2) = 6.226104e-5;
    CONSTANTS(:,15) = 1.8;
    STATES(:,3) = 5;
    CONSTANTS(:,16) = 1;
    CONSTANTS(:,17) = 0.00427;
    CONSTANTS(:,18) = 45;
    CONSTANTS(:,19) = 0.5927;
    STATES(:,4) = 0.009508;
    CONSTANTS(:,20) = 0;
    CONSTANTS(:,21) = shift_Ify;
    CONSTANTS(:,22) = 1.4;
    CONSTANTS(:,23) = 14;
    CONSTANTS(:,24) = 0.08105;
    CONSTANTS(:,25) = 3.343;
    CONSTANTS(:,26) = 0.1369;
    CONSTANTS(:,27) = 0.4315;
    CONSTANTS(:,28) = 0;
    CONSTANTS(:,29) = 26.44;
    CONSTANTS(:,30) = 0.0207;
    CONSTANTS(:,31) = 395.3;
    CONSTANTS(:,32) = 2.289;
    CONSTANTS(:,33) = 26.44;
    CONSTANTS(:,34) = 4.663;
    CONSTANTS(:,35) = 1628;
    CONSTANTS(:,36) = 561.4;
    CONSTANTS(:,37) = 3.663;
    CONSTANTS(:,38) = 0;
    CONSTANTS(:,39) = 0.0223;
    CONSTANTS(:,40) = 0;
    STATES(:,5) = 0.447724;
    STATES(:,6) = 0.003058;
    CONSTANTS(:,41) = 1e-5;
    CONSTANTS(:,42) = P_CaL;
    STATES(:,7) = 0.001921;
    STATES(:,8) = 0.846702;
    STATES(:,9) = 0.844449;
    CONSTANTS(:,43) = 4.3371;
    CONSTANTS(:,44) = -16.4508;
    CONSTANTS(:,45) = 0;
    CONSTANTS(:,46) = 0;
    CONSTANTS(:,47) = 0.0075;
    CONSTANTS(:,48) = 0.000338;
    CONSTANTS(:,49) = 0.04132;
    STATES(:,10) = 0.268909;
    STATES(:,11) = 0.020484;
    CONSTANTS(:,50) = 0;
    STATES(:,12) = 0.9308;
    STATES(:,13) = 6.181512e-9;
    STATES(:,14) = 4.595622e-10;
    STATES(:,15) = 0.069199;
    CONSTANTS(:,51) = 148041085.1;
    CONSTANTS(:,52) = 15;
    CONSTANTS(:,53) = 1;
    CONSTANTS(:,54) = 0.45;
    CONSTANTS(:,55) = 2.5;
    CONSTANTS(:,56) = 10000;
    CONSTANTS(:,57) = 500;
    CONSTANTS(:,58) = 5;
    CONSTANTS(:,59) = 660;
    STATES(:,16) = 0.409551;
    CONSTANTS(:,60) = 5.469e-5;
    CONSTANTS(:,61) = 0.04;
    CONSTANTS(:,62) = Pup_basal;
    CONSTANTS(:,63) = 0.000286113;
    STATES(:,17) = 0.435148;
    STATES(:,18) = 9.15641e-6;
    CONSTANTS(:,64) = 5e-5;
    CONSTANTS(:,65) = 0.031;
    CONSTANTS(:,66) = 0.062;
    CONSTANTS(:,67) = 0.045;
    CONSTANTS(:,68) = 10;
    STATES(:,19) = 0.653777;
    STATES(:,20) = 0.217311;
    STATES(:,21) = 0.158521;
    STATES(:,22) = 0.017929;
    STATES(:,23) = 0.259947;
    STATES(:,24) = 0.138975;
    CONSTANTS(:,69) = 88800;
    CONSTANTS(:,70) = 2277;
    CONSTANTS(:,71) = 227700;
    CONSTANTS(:,72) = 1.642e6;
    CONSTANTS(:,73) = 175.4;
    CONSTANTS(:,74) = 446;
    CONSTANTS(:,75) = 7.51;
    CONSTANTS(:,76) = 751;
    CONSTANTS(:,77) = 542;
    CONSTANTS(:,78) = 445;
    CONSTANTS(:,79) = 2.5;
    CONSTANTS(:,80) = 0.0012;
    CONSTANTS(:,81) = 0.46;
    CONSTANTS(:,82) = 0.0116;
    CONSTANTS(:,83) = 3.9;
    CONSTANTS(:,84) = 67;
    CONSTANTS(:,85) = 0.02;
    CONSTANTS(:,86) = 0.1539e-3;
    STATES(:,25) = 0.011845;
    STATES(:,26) = 0.845304;
    CONSTANTS(:,87) = 3.5e-3;
    STATES(:,27) = 0.430836;
    STATES(:,28) = 0.014523;
    STATES(:,29) = 0.283185;
    STATES(:,30) = 0.011068;
    STATES(:,31) = 0.709051;
    CONSTANTS(:,89) = gKs;
    STATES(:,32) = 0.1162;
    CONSTANTS(:,90) = 0.00345;
    STATES(:,33) = 0.00277;
    CONSTANTS(:,91) = 1;
    CONSTANTS(:,92) = ( CONSTANTS(:,1).*CONSTANTS(:,2))./CONSTANTS(:,3);
    CONSTANTS(:,93) = CONSTANTS(:,17)./(CONSTANTS(:,14)./(CONSTANTS(:,14)+CONSTANTS(:,18)));
%     deleted -0.25 from the response to iso - coded in separately through
%     P_up_basal (CO62)
%     CONSTANTS(:,94) = piecewise({CONSTANTS(:,11)>0.00000,  - 0.250000 , CONSTANTS(:,10)>0.00000, ( 0.700000.*CONSTANTS(:,10))./(9.00000e-05+CONSTANTS(:,10)) }, 0.00000);
    CONSTANTS(:,94) = piecewise({CONSTANTS(:,11)>0.00000,  0 , CONSTANTS(:,10)>0.00000, ( 0.700000.*CONSTANTS(:,10))./(9.00000e-05+CONSTANTS(:,10)) }, 0.00000);
    CONSTANTS(:,95) = piecewise({CONSTANTS(:,11)>0.00000,  1.20000.*CONSTANTS(:,89) }, CONSTANTS(:,89));
    CONSTANTS(:,96) = piecewise({CONSTANTS(:,10)>0.00000,  - 1.00000 - ( 9.89800.*power( 1.00000.*CONSTANTS(:,10), 0.618000))./(power( 1.00000.*CONSTANTS(:,10), 0.618000)+0.00122423) }, 0.00000);
    CONSTANTS(:,97) =  CONSTANTS(:,92).*log(CONSTANTS(:,14)./CONSTANTS(:,13));
    CONSTANTS(:,98) = CONSTANTS(:,93)./(CONSTANTS(:,19)+1.00000);
    CONSTANTS(:,99) =  CONSTANTS(:,62).*(1.00000 - CONSTANTS(:,94));
%     CONSTANTS(:,100) = piecewise({CONSTANTS(:,11)>0.00000, 7.50000 }, 0.00000);
    CONSTANTS(:,100) = Isoshift_Ify; %piecewise({CONSTANTS(:,11)>0.00000, 7.50000 }, 0.00000);
    CONSTANTS(:,101) =  CONSTANTS(:,19).*CONSTANTS(:,98);
    CONSTANTS(:,102) = ( CONSTANTS(:,98).*CONSTANTS(:,14))./(CONSTANTS(:,14)+CONSTANTS(:,18));
    CONSTANTS(:,103) = Iso_i_NaK;%piecewise({CONSTANTS(:,11)>0.00000, 1.20000 }, 1.00000);
    CONSTANTS(:,104) = ( CONSTANTS(:,101).*CONSTANTS(:,14))./(CONSTANTS(:,14)+CONSTANTS(:,18));
    CONSTANTS(:,105) = CONSTANTS(:,12)./(CONSTANTS(:,34)+CONSTANTS(:,12));
%     CONSTANTS(:,106) = piecewise({CONSTANTS(:,11)>0.00000, 1.23000 }, 1.00000);
    CONSTANTS(:,106) = 1; % changed and changes are coded in through the maximum conductance P_CaL (CO42)
    CONSTANTS(:,107) = ( 0.310000.*CONSTANTS(:,10))./(CONSTANTS(:,10)+9.00000e-05);
%     CONSTANTS(:,108) = piecewise({CONSTANTS(:,11)>0.00000,  - 8.00000 }, 0.00000);
    CONSTANTS(:,108) = Isoshift_dL;% piecewise({CONSTANTS(:,11)>0.00000,  - 8.00000 }, 0.00000);
%     CONSTANTS(:,109) = piecewise({CONSTANTS(:,11)>0.00000,  - 27.0000 }, 0.00000);
    CONSTANTS(:,109) = Isoshift_slopedL; %piecewise({CONSTANTS(:,11)>0.00000,  - 27.0000 }, 0.00000);
    CONSTANTS(:,110) =  1.00000e-09.* pi.*power(CONSTANTS(:,83), 2.00000).*CONSTANTS(:,84);
    CONSTANTS(:,111) =  1.00000e-09.*2.00000.* pi.*CONSTANTS(:,85).*(CONSTANTS(:,83) - CONSTANTS(:,85)./2.00000).*CONSTANTS(:,84);
    CONSTANTS(:,112) =  CONSTANTS(:,80).*CONSTANTS(:,110);
    CONSTANTS(:,113) =  CONSTANTS(:,81).*CONSTANTS(:,110) - CONSTANTS(:,111);
    CONSTANTS(:,114) =  CONSTANTS(:,82).*CONSTANTS(:,110);
%     CONSTANTS(:,115) = piecewise({CONSTANTS(:,11)>0.00000,  - 14.0000 }, 0.00000);
    CONSTANTS(:,115) = Isoshift_ninf; %piecewise({CONSTANTS(:,11)>0.00000,  - 14.0000 }, 0.00000); IKs n gate
    CONSTANTS(:,116) = (3.59880 - 0.0256410)./(1.00000+1.21550e-06./power( 1.00000.*CONSTANTS(:,10), 1.69510))+0.0256410;
    CONSTANTS(:,88) = 0.00424;
    % IKr - effect of glucose
    CONSTANTS(:,117) = x_g_Kr; % multiplication factor for gKr (default: 1)    
    CONSTANTS(:,118) = IHerg_shift_pa; % shift in activation curve IKr, gate pa (default: 0)    
    % If - effect of mutations
    CONSTANTS(:,119) = 0; % shift in activation curve If (-43.9)(default: 0)    
    CONSTANTS(:,120) = 0; % shift in time constant If (-43.9)(default: 0)    
    CONSTANTS(:,121) = 1; % factor max conductance in If due to mutation (0.34)(default: 1)    
    % Na activation gate m
    CONSTANTS(:,122) = 0; % shift in activation INa m gate (6 to 12.5) (default: 0)    
    CONSTANTS(:,123) = 0; % shift in time constant INa m gate (6 to 12.5) (default: 0)    
    CONSTANTS(:,124) = 1; % change in slope in activation INa m gate (18 to 78%: 1.18 to 1.78) (default: 1)    
    % Na inactivation gate h
    CONSTANTS(:,125) = 0; % shift in inactivation INa h gate (-9.7 to -15) (default: 0)    
%     CONSTANTS(:,126) = 0; % shift in time constant INa h gate (x2) (default: )    
    CONSTANTS(:,127) = 1; % change in slope in inactivation INa h gate (-6.8%: 1 to 0.93) (default: 1)    
    CONSTANTS(:,128) = 1; % change in current density INa (-40% to -100%) (default: 1)    
    % IKs - effect of mutations
    CONSTANTS(:,129) = 0; % shift in activation curve IKs, gate n; %(-43 to 25) (default: 0)    
%     CONSTANTS(:,130) = 0; % shift in time constant IKS, gate n (default: )    
    CONSTANTS(:,131) = 1; % shift in activation slope %( 1 to 1.53) (default: 1)        
    CONSTANTS(:,132) = 1; % (0.3) % current density IKs (-70%)(default: 1)
    
    %%%%
    % check also pa gate shift -- CONSTANTS(:,118);
    CONSTANTS(:,133) = IHerg_shift_pi; % shift in activation curve IKr, gate pi (default: 0)    
    
    
    if (isempty(STATES)), warning('Initial values for states not set'); end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)

    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,7) =  CONSTANTS(:,70).*CONSTANTS(:,79).*(1.00000 - (STATES(:,23)+STATES(:,19))) -  CONSTANTS(:,76).*STATES(:,19);
    RATES(:,19) = ALGEBRAIC(:,7);
    ALGEBRAIC(:,2) = CONSTANTS(:,48)./(CONSTANTS(:,48)+STATES(:,2));
    ALGEBRAIC(:,8) = ( 0.00100000.*ALGEBRAIC(:,2))./CONSTANTS(:,47);
    RATES(:,9) = (ALGEBRAIC(:,2) - STATES(:,9))./ALGEBRAIC(:,8);
    ALGEBRAIC(:,3) = CONSTANTS(:,52) - (CONSTANTS(:,52) - CONSTANTS(:,53))./(1.00000+power(CONSTANTS(:,54)./STATES(:,16), CONSTANTS(:,55)));
    ALGEBRAIC(:,9) = CONSTANTS(:,56)./ALGEBRAIC(:,3);
    ALGEBRAIC(:,18) =  CONSTANTS(:,57).*ALGEBRAIC(:,3);
    RATES(:,12) = ( CONSTANTS(:,58).*STATES(:,15) -  ALGEBRAIC(:,18).*STATES(:,2).*STATES(:,12)) - ( ALGEBRAIC(:,9).*power(STATES(:,2), 2.00000).*STATES(:,12) -  CONSTANTS(:,59).*STATES(:,13));
    RATES(:,13) = ( ALGEBRAIC(:,9).*power(STATES(:,2), 2.00000).*STATES(:,12) -  CONSTANTS(:,59).*STATES(:,13)) - ( ALGEBRAIC(:,18).*STATES(:,2).*STATES(:,13) -  CONSTANTS(:,58).*STATES(:,14));
    RATES(:,14) = ( ALGEBRAIC(:,18).*STATES(:,2).*STATES(:,13) -  CONSTANTS(:,58).*STATES(:,14)) - ( CONSTANTS(:,59).*STATES(:,14) -  ALGEBRAIC(:,9).*power(STATES(:,2), 2.00000).*STATES(:,15));
    RATES(:,15) = ( CONSTANTS(:,59).*STATES(:,14) -  ALGEBRAIC(:,9).*power(STATES(:,2), 2.00000).*STATES(:,15)) - ( CONSTANTS(:,58).*STATES(:,15) -  ALGEBRAIC(:,18).*STATES(:,2).*STATES(:,12));
    ALGEBRAIC(:,6) = piecewise({VOI>CONSTANTS(:,6)&VOI<CONSTANTS(:,6)+CONSTANTS(:,7), CONSTANTS(:,8) }, CONSTANTS(:,9));
    ALGEBRAIC(:,10) = piecewise({CONSTANTS(:,5)>=1.00000, ALGEBRAIC(:,6) }, STATES(:,1));
%     does not include CONSTANTS(:,21)
%     ALGEBRAIC(:,11) = 1.00000./(( 0.360000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100)))./(exp( 0.0660000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100))) - 1.00000)+( 0.100000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100)))./(1.00000 - exp(  - 0.200000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100))))) - 0.0540000;
    ALGEBRAIC(:,11) = 1.00000./(( 0.360000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120)))./(exp( 0.0660000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100) - CONSTANTS(:,21)- CONSTANTS(:,120))) - 1.00000)+( 0.100000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120)))./(1.00000 - exp(  - 0.200000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120))))) - 0.0540000;
%     ALGEBRAIC(:,30) = piecewise({ALGEBRAIC(:,10)< - (((80.0000 - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21)), 0.0132900+0.999210./(1.00000+exp(((((ALGEBRAIC(:,10)+97.1340) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21))./8.17520)) },  0.000250100.*exp( - (((ALGEBRAIC(:,10) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21))./12.8610));
    ALGEBRAIC(:,30) = piecewise({ALGEBRAIC(:,10)< - (((80.0000 - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21) - CONSTANTS(:,119)), 0.0132900+0.999210./(1.00000+exp(((((ALGEBRAIC(:,10)+97.1340) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21) - CONSTANTS(:,119))./8.17520)) },  0.000250100.*exp( - (((ALGEBRAIC(:,10) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21)- CONSTANTS(:,119))./12.8610));
    RATES(:,4) = (ALGEBRAIC(:,30) - STATES(:,4))./ALGEBRAIC(:,11);
    ALGEBRAIC(:,15) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+37.4000+CONSTANTS(:,45))./(5.30000+CONSTANTS(:,46))));
    ALGEBRAIC(:,34) =  0.00100000.*(44.3000+ 230.000.*exp( - power((ALGEBRAIC(:,10)+36.0000)./10.0000, 2.00000)));
    RATES(:,8) = (ALGEBRAIC(:,15) - STATES(:,8))./ALGEBRAIC(:,34);
    ALGEBRAIC(:,16) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+38.3000)./5.50000));
    ALGEBRAIC(:,35) = 0.00100000./( 1.06800.*exp((ALGEBRAIC(:,10)+38.3000)./30.0000)+ 1.06800.*exp( - (ALGEBRAIC(:,10)+38.3000)./30.0000));
    RATES(:,10) = (ALGEBRAIC(:,16) - STATES(:,10))./ALGEBRAIC(:,35);
    ALGEBRAIC(:,17) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+58.7000)./3.80000));
    ALGEBRAIC(:,36) = 1.00000./( 16.6700.*exp( - (ALGEBRAIC(:,10)+75.0000)./83.3000)+ 16.6700.*exp((ALGEBRAIC(:,10)+75.0000)./15.3800))+CONSTANTS(:,50);
    RATES(:,11) = (ALGEBRAIC(:,17) - STATES(:,11))./ALGEBRAIC(:,36);
    ALGEBRAIC(:,38) = 0.00900000./(1.00000+exp((ALGEBRAIC(:,10)+5.00000)./12.0000))+0.000500000;
    ALGEBRAIC(:,20) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+6.00000)./ - 8.60000));
    RATES(:,25) = (ALGEBRAIC(:,20) - STATES(:,25))./ALGEBRAIC(:,38);
    ALGEBRAIC(:,39) = 0.590000./(1.00000+exp((ALGEBRAIC(:,10)+60.0000)./10.0000))+3.05000;
    ALGEBRAIC(:,21) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+7.50000)./10.0000));
    RATES(:,26) = (ALGEBRAIC(:,21) - STATES(:,26))./ALGEBRAIC(:,39);
    ALGEBRAIC(:,22) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+49.0000)./13.0000));
    ALGEBRAIC(:,40) =  0.00100000.*0.600000.*(65.1700./( 0.570000.*exp(  - 0.0800000.*(ALGEBRAIC(:,10)+44.0000))+ 0.0650000.*exp( 0.100000.*(ALGEBRAIC(:,10)+45.9300)))+10.1000);
    RATES(:,27) = (ALGEBRAIC(:,22) - STATES(:,27))./ALGEBRAIC(:,40);
    ALGEBRAIC(:,23) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10) - 19.3000)./15.0000));
    ALGEBRAIC(:,41) =  0.00100000.*0.660000.*1.40000.*(15.5900./( 1.03700.*exp( 0.0900000.*(ALGEBRAIC(:,10)+30.6100))+ 0.369000.*exp(  - 0.120000.*(ALGEBRAIC(:,10)+23.8400)))+2.98000);
    RATES(:,28) = (ALGEBRAIC(:,23) - STATES(:,28))./ALGEBRAIC(:,41);
%     ALGEBRAIC(:,24) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+10.0144)./7.66070));    
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,24) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+10.0144-CONSTANTS(:,118))./7.66070));
%     ALGEBRAIC(:,42) = 0.846554./( 4.20000.*exp(ALGEBRAIC(:,10)./17.0000)+ 0.150000.*exp( - ALGEBRAIC(:,10)./21.6000));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,42) = 0.846554./( 4.20000.*exp((ALGEBRAIC(:,10)-CONSTANTS(:,118))./17.0000)+ 0.150000.*exp( - (ALGEBRAIC(:,10)-CONSTANTS(:,118))./21.6000));
    RATES(:,29) = (ALGEBRAIC(:,24) - STATES(:,29))./ALGEBRAIC(:,42);
%     ALGEBRAIC(:,43) = 1.00000./( 30.0000.*exp(ALGEBRAIC(:,10)./10.0000)+exp( - ALGEBRAIC(:,10)./12.0000));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,43) = 1.00000./( 30.0000.*exp((ALGEBRAIC(:,10)-CONSTANTS(:,118))./10.0000)+exp( - (ALGEBRAIC(:,10)-CONSTANTS(:,118))./12.0000));
    RATES(:,30) = (ALGEBRAIC(:,24) - STATES(:,30))./ALGEBRAIC(:,43);
%     ALGEBRAIC(:,44) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+28.6000)./17.1000));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,44) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+28.6000-CONSTANTS(:,133))./17.1000));
%     ALGEBRAIC(:,27) = 1.00000./( 100.000.*exp( - ALGEBRAIC(:,10)./54.6450)+ 656.000.*exp(ALGEBRAIC(:,10)./106.157));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,27) = 1.00000./( 100.000.*exp( - (ALGEBRAIC(:,10)-CONSTANTS(:,133))./54.6450)+ 656.000.*exp((ALGEBRAIC(:,10)-CONSTANTS(:,133))./106.157));
    RATES(:,31) = (ALGEBRAIC(:,44) - STATES(:,31))./ALGEBRAIC(:,27);
    ALGEBRAIC(:,29) =  10.0000.*exp( 0.0133000.*(ALGEBRAIC(:,10)+40.0000));
    ALGEBRAIC(:,46) = CONSTANTS(:,116)./(CONSTANTS(:,116)+ALGEBRAIC(:,29));
    ALGEBRAIC(:,52) = 1.00000./(CONSTANTS(:,116)+ALGEBRAIC(:,29));
    RATES(:,33) = (ALGEBRAIC(:,46) - STATES(:,33))./ALGEBRAIC(:,52);
%     ALGEBRAIC(:,13) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+69.8040)./4.45650));
    ALGEBRAIC(:,13) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+69.8040-CONSTANTS(:,125))./(4.45650*CONSTANTS(:,127))));
%     ALGEBRAIC(:,32) =  20.0000.*exp(  - 0.125000.*(ALGEBRAIC(:,10)+75.0000));
    ALGEBRAIC(:,32) =  20.0000.*exp(  - 0.125000.*(ALGEBRAIC(:,10)+75.0000-CONSTANTS(:,126)));
%     ALGEBRAIC(:,48) = 2000.00./( 320.000.*exp(  - 0.100000.*(ALGEBRAIC(:,10)+75.0000))+1.00000);
    ALGEBRAIC(:,48) = 2000.00./( 320.000.*exp(  - 0.100000.*(ALGEBRAIC(:,10)+75.0000-CONSTANTS(:,126)))+1.00000);
    ALGEBRAIC(:,54) = 1.00000./(ALGEBRAIC(:,32)+ALGEBRAIC(:,48));
    RATES(:,6) = (ALGEBRAIC(:,13) - STATES(:,6))./ALGEBRAIC(:,54);
%     ALGEBRAIC(:,28) = power((1.00000./(1.00000+exp( - ((ALGEBRAIC(:,10)+0.638300) - CONSTANTS(:,115))./10.7071))), 1.0 ./ 2);
    ALGEBRAIC(:,28) = power((1.00000./(1.00000+exp( - ((ALGEBRAIC(:,10)+0.638300) - CONSTANTS(:,115)-CONSTANTS(:,129))./(10.7071*CONSTANTS(:,131))))), 1.0 ./ 2);
%     ALGEBRAIC(:,45) = 28.0000./(1.00000+exp( - ((ALGEBRAIC(:,10) - 40.0000) - CONSTANTS(:,115))./3.00000));
    ALGEBRAIC(:,45) = 28.0000./(1.00000+exp( - ((ALGEBRAIC(:,10) - 40.0000) - CONSTANTS(:,115)-CONSTANTS(:,130))./3.00000));
%     ALGEBRAIC(:,51) =  1.00000.*exp( - ((ALGEBRAIC(:,10) - CONSTANTS(:,115)) - 5.00000)./25.0000);
    ALGEBRAIC(:,51) =  1.00000.*exp( - ((ALGEBRAIC(:,10) - CONSTANTS(:,115)-CONSTANTS(:,130)) - 5.00000)./25.0000);
    ALGEBRAIC(:,57) = 1.00000./(ALGEBRAIC(:,45)+ALGEBRAIC(:,51));
    RATES(:,32) = (ALGEBRAIC(:,28) - STATES(:,32))./ALGEBRAIC(:,57);
%     ALGEBRAIC(:,12) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+42.0504)./8.31060));
    ALGEBRAIC(:,12) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+42.0504-CONSTANTS(:,122))./(8.31060*CONSTANTS(:,124))));
    ALGEBRAIC(:,31) = ALGEBRAIC(:,10)+41.0000-CONSTANTS(:,123);
    ALGEBRAIC(:,47) = piecewise({abs(ALGEBRAIC(:,31))<CONSTANTS(:,41), 2000.00 }, ( 200.000.*ALGEBRAIC(:,31))./(1.00000 - exp(  - 0.100000.*ALGEBRAIC(:,31))));
%     ALGEBRAIC(:,53) =  8000.00.*exp(  - 0.0560000.*(ALGEBRAIC(:,10)+66.0000));
    ALGEBRAIC(:,53) =  8000.00.*exp(  - 0.0560000.*(ALGEBRAIC(:,10)+66.0000-CONSTANTS(:,123)));
    ALGEBRAIC(:,58) = 1.00000./(ALGEBRAIC(:,47)+ALGEBRAIC(:,53));
    RATES(:,5) = (ALGEBRAIC(:,12) - STATES(:,5))./ALGEBRAIC(:,58);
    ALGEBRAIC(:,14) = 1.00000./(1.00000+exp( - ((ALGEBRAIC(:,10) - CONSTANTS(:,44)) - CONSTANTS(:,108))./( CONSTANTS(:,43).*(1.00000+CONSTANTS(:,109)./100.000))));
    ALGEBRAIC(:,33) = piecewise({ALGEBRAIC(:,10)== - 41.8000,  - 41.8000 , ALGEBRAIC(:,10)==0.00000, 0.00000 , ALGEBRAIC(:,10)== - 6.80000,  - 6.80001 }, ALGEBRAIC(:,10));
%     ALGEBRAIC(:,49) = (  - 0.0283900.*(ALGEBRAIC(:,33)+41.8000))./(exp( - (ALGEBRAIC(:,33)+41.8000)./2.50000) - 1.00000) - ( 0.0849000.*(ALGEBRAIC(:,33)+6.80000))./(exp( - (ALGEBRAIC(:,33)+6.80000)./4.80000) - 1.00000);
    ALGEBRAIC(:,49) = (  - 0.0283900.*(ALGEBRAIC(:,33)+41.8000-CONSTANTS(:,108)))./(exp( - (ALGEBRAIC(:,33)+41.8000-CONSTANTS(:,108))./2.50000) - 1.00000) - ( 0.0849000.*(ALGEBRAIC(:,33)+6.80000-CONSTANTS(:,108)))./(exp( - (ALGEBRAIC(:,33)+6.80000-CONSTANTS(:,108))./4.80000) - 1.00000);
    ALGEBRAIC(:,55) = piecewise({ALGEBRAIC(:,10)== - 1.80000,  - 1.80001 }, ALGEBRAIC(:,10));
%     ALGEBRAIC(:,59) = ( 0.0114300.*(ALGEBRAIC(:,55)+1.80000))./(exp((ALGEBRAIC(:,55)+1.80000)./2.50000) - 1.00000);
    ALGEBRAIC(:,59) = ( 0.0114300.*(ALGEBRAIC(:,55)+1.80000-CONSTANTS(:,108)))./(exp((ALGEBRAIC(:,55)+1.80000-CONSTANTS(:,108))./2.50000) - 1.00000);
    ALGEBRAIC(:,61) = 0.00100000./(ALGEBRAIC(:,49)+ALGEBRAIC(:,59));
    RATES(:,7) = (ALGEBRAIC(:,14) - STATES(:,7))./ALGEBRAIC(:,61);
    ALGEBRAIC(:,19) = STATES(:,3);
    ALGEBRAIC(:,37) =  CONSTANTS(:,92).*log(CONSTANTS(:,12)./ALGEBRAIC(:,19));
    ALGEBRAIC(:,62) =  CONSTANTS(:,103).*CONSTANTS(:,24).*power(1.00000+power(CONSTANTS(:,22)./CONSTANTS(:,14), 1.20000),  - 1.00000).*power(1.00000+power(CONSTANTS(:,23)./ALGEBRAIC(:,19), 1.30000),  - 1.00000).*power(1.00000+exp( - ((ALGEBRAIC(:,10) - ALGEBRAIC(:,37))+110.000)./20.0000),  - 1.00000);
    ALGEBRAIC(:,64) = exp((  - CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92)));
    ALGEBRAIC(:,70) = 1.00000+ (CONSTANTS(:,15)./CONSTANTS(:,37)).*(1.00000+exp(( CONSTANTS(:,28).*ALGEBRAIC(:,10))./CONSTANTS(:,92)))+ (CONSTANTS(:,12)./CONSTANTS(:,35)).*(1.00000+ (CONSTANTS(:,12)./CONSTANTS(:,36)).*(1.00000+CONSTANTS(:,12)./CONSTANTS(:,34)));
    ALGEBRAIC(:,72) = ( (( (CONSTANTS(:,12)./CONSTANTS(:,35)).*CONSTANTS(:,12))./CONSTANTS(:,36)).*(1.00000+CONSTANTS(:,12)./CONSTANTS(:,34)).*exp((  - CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92))))./ALGEBRAIC(:,70);
    ALGEBRAIC(:,71) = ( (CONSTANTS(:,15)./CONSTANTS(:,37)).*exp(( CONSTANTS(:,28).*ALGEBRAIC(:,10))./CONSTANTS(:,92)))./ALGEBRAIC(:,70);
    ALGEBRAIC(:,68) = exp(( CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92)));
    ALGEBRAIC(:,63) = ALGEBRAIC(:,19)./(CONSTANTS(:,29)+ALGEBRAIC(:,19));
    ALGEBRAIC(:,73) =  ALGEBRAIC(:,64).*CONSTANTS(:,105).*(ALGEBRAIC(:,72)+ALGEBRAIC(:,71))+ ALGEBRAIC(:,71).*ALGEBRAIC(:,68).*(ALGEBRAIC(:,63)+ALGEBRAIC(:,64));
    ALGEBRAIC(:,65) = 1.00000+ (STATES(:,2)./CONSTANTS(:,30)).*(1.00000+exp((  - CONSTANTS(:,26).*ALGEBRAIC(:,10))./CONSTANTS(:,92))+ALGEBRAIC(:,19)./CONSTANTS(:,33))+ (ALGEBRAIC(:,19)./CONSTANTS(:,31)).*(1.00000+ (ALGEBRAIC(:,19)./CONSTANTS(:,32)).*(1.00000+ALGEBRAIC(:,19)./CONSTANTS(:,29)));
    ALGEBRAIC(:,66) = ( (STATES(:,2)./CONSTANTS(:,30)).*exp((  - CONSTANTS(:,26).*ALGEBRAIC(:,10))./CONSTANTS(:,92)))./ALGEBRAIC(:,65);
    ALGEBRAIC(:,67) = ( (( (ALGEBRAIC(:,19)./CONSTANTS(:,31)).*ALGEBRAIC(:,19))./CONSTANTS(:,32)).*(1.00000+ALGEBRAIC(:,19)./CONSTANTS(:,29)).*exp(( CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92))))./ALGEBRAIC(:,65);
    ALGEBRAIC(:,69) =  ALGEBRAIC(:,68).*ALGEBRAIC(:,63).*(ALGEBRAIC(:,67)+ALGEBRAIC(:,66))+ ALGEBRAIC(:,64).*ALGEBRAIC(:,66).*(CONSTANTS(:,105)+ALGEBRAIC(:,68));
    ALGEBRAIC(:,74) =  ALGEBRAIC(:,67).*ALGEBRAIC(:,63).*(ALGEBRAIC(:,72)+ALGEBRAIC(:,71))+ ALGEBRAIC(:,66).*ALGEBRAIC(:,72).*(ALGEBRAIC(:,63)+ALGEBRAIC(:,64));
    ALGEBRAIC(:,75) =  ALGEBRAIC(:,72).*CONSTANTS(:,105).*(ALGEBRAIC(:,67)+ALGEBRAIC(:,66))+ ALGEBRAIC(:,67).*ALGEBRAIC(:,71).*(CONSTANTS(:,105)+ALGEBRAIC(:,68));
    ALGEBRAIC(:,76) = ( (1.00000 - CONSTANTS(:,38)).*CONSTANTS(:,25).*( ALGEBRAIC(:,69).*ALGEBRAIC(:,71) -  ALGEBRAIC(:,73).*ALGEBRAIC(:,66)))./(ALGEBRAIC(:,73)+ALGEBRAIC(:,69)+ALGEBRAIC(:,74)+ALGEBRAIC(:,75));
    ALGEBRAIC(:,77) =  CONSTANTS(:,92).*log((CONSTANTS(:,12)+ 0.120000.*CONSTANTS(:,14))./(ALGEBRAIC(:,19)+ 0.120000.*CONSTANTS(:,13)));
%     ALGEBRAIC(:,78) =  CONSTANTS(:,39).*power(STATES(:,5), 3.00000).*STATES(:,6).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,77));
    ALGEBRAIC(:,78) =  CONSTANTS(:,128).*CONSTANTS(:,39).*power(STATES(:,5), 3.00000).*STATES(:,6).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,77));
    ALGEBRAIC(:,79) =  CONSTANTS(:,40).*power(STATES(:,5), 3.00000).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,77));
    ALGEBRAIC(:,80) = ALGEBRAIC(:,78)+ALGEBRAIC(:,79);
    ALGEBRAIC(:,50) =  STATES(:,4).*CONSTANTS(:,104).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,37)).*(1.00000 - CONSTANTS(:,20)).*CONSTANTS(:,121);
    ALGEBRAIC(:,83) =  (( 1.85000e-05.*CONSTANTS(:,42).*(ALGEBRAIC(:,10) - 0.00000))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))))).*(ALGEBRAIC(:,19) -  CONSTANTS(:,12).*exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))).*STATES(:,7).*STATES(:,8).*STATES(:,9);
    RATES(:,3) = ( (1.00000 - CONSTANTS(:,16)).* - 1.00000.*(ALGEBRAIC(:,80)+ALGEBRAIC(:,50)+ALGEBRAIC(:,83)+ 3.00000.*ALGEBRAIC(:,62)+ 3.00000.*ALGEBRAIC(:,76)))./( 1.00000.*(CONSTANTS(:,113)+CONSTANTS(:,111)).*CONSTANTS(:,3));
    ALGEBRAIC(:,91) =  CONSTANTS(:,72).*STATES(:,2).*(1.00000 - STATES(:,21)) -  CONSTANTS(:,77).*STATES(:,21);
    RATES(:,21) = ALGEBRAIC(:,91);
    ALGEBRAIC(:,85) =  (( 2.00000.*CONSTANTS(:,49).*ALGEBRAIC(:,10))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*ALGEBRAIC(:,10).*2.00000)./CONSTANTS(:,92))))).*(STATES(:,2) -  CONSTANTS(:,15).*exp((  - 2.00000.*ALGEBRAIC(:,10))./CONSTANTS(:,92))).*STATES(:,10).*STATES(:,11);
    ALGEBRAIC(:,81) =  (( 2.00000.*CONSTANTS(:,42).*(ALGEBRAIC(:,10) - 0.00000))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000).*2.00000)./CONSTANTS(:,92))))).*(STATES(:,2) -  CONSTANTS(:,15).*exp((  - 2.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))).*STATES(:,7).*STATES(:,8).*STATES(:,9);
    ALGEBRAIC(:,87) =  CONSTANTS(:,51).*STATES(:,13).*(STATES(:,16) - STATES(:,2));
    ALGEBRAIC(:,89) = (STATES(:,2) - STATES(:,18))./CONSTANTS(:,60);
    RATES(:,2) = ( ALGEBRAIC(:,87).*CONSTANTS(:,112))./CONSTANTS(:,111) - (((ALGEBRAIC(:,81)+ALGEBRAIC(:,85)) -  2.00000.*ALGEBRAIC(:,76))./( 2.00000.*CONSTANTS(:,3).*CONSTANTS(:,111))+ALGEBRAIC(:,89)+ CONSTANTS(:,67).*ALGEBRAIC(:,91));
    ALGEBRAIC(:,94) =  CONSTANTS(:,69).*STATES(:,18).*(1.00000 - STATES(:,22)) -  CONSTANTS(:,74).*STATES(:,22);
    RATES(:,22) = ALGEBRAIC(:,94);
    ALGEBRAIC(:,92) = CONSTANTS(:,99)./(1.00000+exp(( - STATES(:,18)+CONSTANTS(:,63))./CONSTANTS(:,64)));
    ALGEBRAIC(:,95) = (STATES(:,17) - STATES(:,16))./CONSTANTS(:,61);
    RATES(:,17) = ALGEBRAIC(:,92) - ( ALGEBRAIC(:,95).*CONSTANTS(:,112))./CONSTANTS(:,114);
    ALGEBRAIC(:,97) =  CONSTANTS(:,71).*STATES(:,18).*(1.00000 - (STATES(:,23)+STATES(:,19))) -  CONSTANTS(:,75).*STATES(:,23);
    RATES(:,23) = ALGEBRAIC(:,97);
    ALGEBRAIC(:,98) =  CONSTANTS(:,73).*STATES(:,16).*(1.00000 - STATES(:,24)) -  CONSTANTS(:,78).*STATES(:,24);
    RATES(:,24) = ALGEBRAIC(:,98);
    RATES(:,16) = ALGEBRAIC(:,95) - (ALGEBRAIC(:,87)+ CONSTANTS(:,68).*ALGEBRAIC(:,98));
    ALGEBRAIC(:,100) =  CONSTANTS(:,72).*STATES(:,18).*(1.00000 - STATES(:,20)) -  CONSTANTS(:,77).*STATES(:,20);
    RATES(:,20) = ALGEBRAIC(:,100);
    RATES(:,18) = ( 1.00000.*( ALGEBRAIC(:,89).*CONSTANTS(:,111) -  ALGEBRAIC(:,92).*CONSTANTS(:,114)))./CONSTANTS(:,113) - ( CONSTANTS(:,67).*ALGEBRAIC(:,100)+ CONSTANTS(:,65).*ALGEBRAIC(:,94)+ CONSTANTS(:,66).*ALGEBRAIC(:,97));
    ALGEBRAIC(:,56) =  STATES(:,4).*CONSTANTS(:,102).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*(1.00000 - CONSTANTS(:,20)).*CONSTANTS(:,121);
    ALGEBRAIC(:,60) = ALGEBRAIC(:,50)+ALGEBRAIC(:,56);
%     ALGEBRAIC(:,90) =  CONSTANTS(:,117).*CONSTANTS(:,88).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*( 0.900000.*STATES(:,30)+ 0.100000.*STATES(:,29)).*STATES(:,31);
    ALGEBRAIC(:,90) =  CONSTANTS(:,117).*CONSTANTS(:,88).*(ALGEBRAIC(:,10) - CONSTANTS(:,118) - CONSTANTS(:,97)).*( 0.900000.*STATES(:,30)+ 0.100000.*STATES(:,29)).*STATES(:,31);
    ALGEBRAIC(:,93) =  CONSTANTS(:,92).*log((CONSTANTS(:,14)+ 0.120000.*CONSTANTS(:,12))./(CONSTANTS(:,13)+ 0.120000.*ALGEBRAIC(:,19)));
%     ALGEBRAIC(:,96) =  CONSTANTS(:,95).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,93)).*power(STATES(:,32), 2.00000);
    ALGEBRAIC(:,96) =  CONSTANTS(:,132).*CONSTANTS(:,95).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,93)).*power(STATES(:,32), 2.00000);
    ALGEBRAIC(:,88) =  CONSTANTS(:,87).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*STATES(:,27).*STATES(:,28);
    ALGEBRAIC(:,82) =  (( 0.000365000.*CONSTANTS(:,42).*(ALGEBRAIC(:,10) - 0.00000))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))))).*(CONSTANTS(:,13) -  CONSTANTS(:,14).*exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))).*STATES(:,7).*STATES(:,8).*STATES(:,9);
    ALGEBRAIC(:,84) =  (ALGEBRAIC(:,81)+ALGEBRAIC(:,82)+ALGEBRAIC(:,83)).*(1.00000 - CONSTANTS(:,107)).*1.00000.*CONSTANTS(:,106);
    ALGEBRAIC(:,99) = piecewise({CONSTANTS(:,10)>0.00000,  CONSTANTS(:,91).*CONSTANTS(:,90).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*(1.00000+exp((ALGEBRAIC(:,10)+20.0000)./20.0000)).*STATES(:,33) }, 0.00000);
    ALGEBRAIC(:,86) =  CONSTANTS(:,86).*STATES(:,25).*STATES(:,26).*(ALGEBRAIC(:,10) - CONSTANTS(:,97));
    ALGEBRAIC(:,101) = ALGEBRAIC(:,60)+ALGEBRAIC(:,90)+ALGEBRAIC(:,96)+ALGEBRAIC(:,88)+ALGEBRAIC(:,62)+ALGEBRAIC(:,76)+ALGEBRAIC(:,80)+ALGEBRAIC(:,84)+ALGEBRAIC(:,85)+ALGEBRAIC(:,99)+ALGEBRAIC(:,86);
    RATES(:,1) =  - ALGEBRAIC(:,101)./CONSTANTS(:,4);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)

    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,7) =  CONSTANTS(:,70).*CONSTANTS(:,79).*(1.00000 - (STATES(:,23)+STATES(:,19))) -  CONSTANTS(:,76).*STATES(:,19);
    ALGEBRAIC(:,2) = CONSTANTS(:,48)./(CONSTANTS(:,48)+STATES(:,2));
    ALGEBRAIC(:,8) = ( 0.00100000.*ALGEBRAIC(:,2))./CONSTANTS(:,47);
    ALGEBRAIC(:,3) = CONSTANTS(:,52) - (CONSTANTS(:,52) - CONSTANTS(:,53))./(1.00000+power(CONSTANTS(:,54)./STATES(:,16), CONSTANTS(:,55)));
    ALGEBRAIC(:,9) = CONSTANTS(:,56)./ALGEBRAIC(:,3);
    ALGEBRAIC(:,18) =  CONSTANTS(:,57).*ALGEBRAIC(:,3);
    ALGEBRAIC(:,6) = piecewise({VOI>CONSTANTS(:,6)&VOI<CONSTANTS(:,6)+CONSTANTS(:,7), CONSTANTS(:,8) }, CONSTANTS(:,9));
    ALGEBRAIC(:,10) = piecewise({CONSTANTS(:,5)>=1.00000, ALGEBRAIC(:,6) }, STATES(:,1));
%     ALGEBRAIC(:,11) = 1.00000./(( 0.360000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100)))./(exp( 0.0660000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100))) - 1.00000)+( 0.100000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100)))./(1.00000 - exp(  - 0.200000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100))))) - 0.0540000;
    ALGEBRAIC(:,11) = 1.00000./(( 0.360000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120)))./(exp( 0.0660000.*(((ALGEBRAIC(:,10)+148.800) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120))) - 1.00000)+( 0.100000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120)))./(1.00000 - exp(  - 0.200000.*(((ALGEBRAIC(:,10)+87.3000) - CONSTANTS(:,96)) - CONSTANTS(:,100)- CONSTANTS(:,21)- CONSTANTS(:,120))))) - 0.0540000;
%     ALGEBRAIC(:,30) = piecewise({ALGEBRAIC(:,10)< - (((80.0000 - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21)), 0.0132900+0.999210./(1.00000+exp(((((ALGEBRAIC(:,10)+97.1340) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21))./8.17520)) },  0.000250100.*exp( - (((ALGEBRAIC(:,10) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21))./12.8610));
    ALGEBRAIC(:,30) = piecewise({ALGEBRAIC(:,10)< - (((80.0000 - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21) - CONSTANTS(:,119)), 0.0132900+0.999210./(1.00000+exp(((((ALGEBRAIC(:,10)+97.1340) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21) - CONSTANTS(:,119))./8.17520)) },  0.000250100.*exp( - (((ALGEBRAIC(:,10) - CONSTANTS(:,96)) - CONSTANTS(:,100)) - CONSTANTS(:,21) - CONSTANTS(:,119))./12.8610));
    ALGEBRAIC(:,15) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+37.4000+CONSTANTS(:,45))./(5.30000+CONSTANTS(:,46))));
    ALGEBRAIC(:,34) =  0.00100000.*(44.3000+ 230.000.*exp( - power((ALGEBRAIC(:,10)+36.0000)./10.0000, 2.00000)));
    ALGEBRAIC(:,16) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+38.3000)./5.50000));
    ALGEBRAIC(:,35) = 0.00100000./( 1.06800.*exp((ALGEBRAIC(:,10)+38.3000)./30.0000)+ 1.06800.*exp( - (ALGEBRAIC(:,10)+38.3000)./30.0000));
    ALGEBRAIC(:,17) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+58.7000)./3.80000));
    ALGEBRAIC(:,36) = 1.00000./( 16.6700.*exp( - (ALGEBRAIC(:,10)+75.0000)./83.3000)+ 16.6700.*exp((ALGEBRAIC(:,10)+75.0000)./15.3800))+CONSTANTS(:,50);
    ALGEBRAIC(:,38) = 0.00900000./(1.00000+exp((ALGEBRAIC(:,10)+5.00000)./12.0000))+0.000500000;
    ALGEBRAIC(:,20) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+6.00000)./ - 8.60000));
    ALGEBRAIC(:,39) = 0.590000./(1.00000+exp((ALGEBRAIC(:,10)+60.0000)./10.0000))+3.05000;
    ALGEBRAIC(:,21) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+7.50000)./10.0000));
    ALGEBRAIC(:,22) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+49.0000)./13.0000));
    ALGEBRAIC(:,40) =  0.00100000.*0.600000.*(65.1700./( 0.570000.*exp(  - 0.0800000.*(ALGEBRAIC(:,10)+44.0000))+ 0.0650000.*exp( 0.100000.*(ALGEBRAIC(:,10)+45.9300)))+10.1000);
    ALGEBRAIC(:,23) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10) - 19.3000)./15.0000));
    ALGEBRAIC(:,41) =  0.00100000.*0.660000.*1.40000.*(15.5900./( 1.03700.*exp( 0.0900000.*(ALGEBRAIC(:,10)+30.6100))+ 0.369000.*exp(  - 0.120000.*(ALGEBRAIC(:,10)+23.8400)))+2.98000);
%     ALGEBRAIC(:,24) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+10.0144)./7.66070));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,24) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+10.0144-CONSTANTS(:,118))./7.66070));
%     ALGEBRAIC(:,42) = 0.846554./( 4.20000.*exp(ALGEBRAIC(:,10)./17.0000)+ 0.150000.*exp( - ALGEBRAIC(:,10)./21.6000));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,42) = 0.846554./( 4.20000.*exp((ALGEBRAIC(:,10)-CONSTANTS(:,118))./17.0000)+ 0.150000.*exp( - (ALGEBRAIC(:,10)-CONSTANTS(:,118))./21.6000));
%     ALGEBRAIC(:,43) = 1.00000./( 30.0000.*exp(ALGEBRAIC(:,10)./10.0000)+exp( - ALGEBRAIC(:,10)./12.0000));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,43) = 1.00000./( 30.0000.*exp((ALGEBRAIC(:,10)-CONSTANTS(:,118))./10.0000)+exp( - (ALGEBRAIC(:,10)-CONSTANTS(:,118))./12.0000));
%     ALGEBRAIC(:,44) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+28.6000)./17.1000));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,44) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+28.6000-CONSTANTS(:,133))./17.1000));
%     ALGEBRAIC(:,27) = 1.00000./( 100.000.*exp( - ALGEBRAIC(:,10)./54.6450)+ 656.000.*exp(ALGEBRAIC(:,10)./106.157));
    %%% added voltage shift due to lowered [glu]o
    ALGEBRAIC(:,27) = 1.00000./( 100.000.*exp( - (ALGEBRAIC(:,10)-CONSTANTS(:,133))./54.6450)+ 656.000.*exp((ALGEBRAIC(:,10)-CONSTANTS(:,133))./106.157));
    ALGEBRAIC(:,29) =  10.0000.*exp( 0.0133000.*(ALGEBRAIC(:,10)+40.0000));
    ALGEBRAIC(:,46) = CONSTANTS(:,116)./(CONSTANTS(:,116)+ALGEBRAIC(:,29));
    ALGEBRAIC(:,52) = 1.00000./(CONSTANTS(:,116)+ALGEBRAIC(:,29));
%     ALGEBRAIC(:,13) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+69.8040)./4.45650));
    ALGEBRAIC(:,13) = 1.00000./(1.00000+exp((ALGEBRAIC(:,10)+69.8040-CONSTANTS(:,125))./(4.45650*CONSTANTS(:,127))));
%     ALGEBRAIC(:,32) =  20.0000.*exp(  - 0.125000.*(ALGEBRAIC(:,10)+75.0000));
    ALGEBRAIC(:,32) =  20.0000.*exp(  - 0.125000.*(ALGEBRAIC(:,10)+75.0000-CONSTANTS(:,126)));
%     ALGEBRAIC(:,48) = 2000.00./( 320.000.*exp(  - 0.100000.*(ALGEBRAIC(:,10)+75.0000))+1.00000);
    ALGEBRAIC(:,48) = 2000.00./( 320.000.*exp(  - 0.100000.*(ALGEBRAIC(:,10)+75.0000-CONSTANTS(:,126)))+1.00000);
    ALGEBRAIC(:,54) = 1.00000./(ALGEBRAIC(:,32)+ALGEBRAIC(:,48));
%     ALGEBRAIC(:,28) = power((1.00000./(1.00000+exp( - ((ALGEBRAIC(:,10)+0.638300) - CONSTANTS(:,115))./10.7071))), 1.0 ./ 2);
    ALGEBRAIC(:,28) = power((1.00000./(1.00000+exp( - ((ALGEBRAIC(:,10)+0.638300) - CONSTANTS(:,115)-CONSTANTS(:,129))./(10.7071*CONSTANTS(:,131))))), 1.0 ./ 2);
%     ALGEBRAIC(:,45) = 28.0000./(1.00000+exp( - ((ALGEBRAIC(:,10) - 40.0000) - CONSTANTS(:,115))./3.00000));
    ALGEBRAIC(:,45) = 28.0000./(1.00000+exp( - ((ALGEBRAIC(:,10) - 40.0000) - CONSTANTS(:,115)-CONSTANTS(:,130))./3.00000));
%     ALGEBRAIC(:,51) =  1.00000.*exp( - ((ALGEBRAIC(:,10) - CONSTANTS(:,115)) - 5.00000)./25.0000);
    ALGEBRAIC(:,51) =  1.00000.*exp( - ((ALGEBRAIC(:,10) - CONSTANTS(:,115)-CONSTANTS(:,130)) - 5.00000)./25.0000);
    ALGEBRAIC(:,57) = 1.00000./(ALGEBRAIC(:,45)+ALGEBRAIC(:,51));
%     ALGEBRAIC(:,12) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+42.0504)./8.31060));
    ALGEBRAIC(:,12) = 1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+42.0504-CONSTANTS(:,122))./(8.31060*CONSTANTS(:,124))));
    ALGEBRAIC(:,31) = ALGEBRAIC(:,10)+41.0000-CONSTANTS(:,123);
    ALGEBRAIC(:,47) = piecewise({abs(ALGEBRAIC(:,31))<CONSTANTS(:,41), 2000.00 }, ( 200.000.*ALGEBRAIC(:,31))./(1.00000 - exp(  - 0.100000.*ALGEBRAIC(:,31))));
%     ALGEBRAIC(:,53) =  8000.00.*exp(  - 0.0560000.*(ALGEBRAIC(:,10)+66.0000));
    ALGEBRAIC(:,53) =  8000.00.*exp(  - 0.0560000.*(ALGEBRAIC(:,10)+66.0000-CONSTANTS(:,123)));
    ALGEBRAIC(:,58) = 1.00000./(ALGEBRAIC(:,47)+ALGEBRAIC(:,53));
    ALGEBRAIC(:,14) = 1.00000./(1.00000+exp( - ((ALGEBRAIC(:,10) - CONSTANTS(:,44)) - CONSTANTS(:,108))./( CONSTANTS(:,43).*(1.00000+CONSTANTS(:,109)./100.000))));
    ALGEBRAIC(:,33) = piecewise({ALGEBRAIC(:,10)== - 41.8000,  - 41.8000 , ALGEBRAIC(:,10)==0.00000, 0.00000 , ALGEBRAIC(:,10)== - 6.80000,  - 6.80001 }, ALGEBRAIC(:,10));
%     ALGEBRAIC(:,49) = (  - 0.0283900.*(ALGEBRAIC(:,33)+41.8000))./(exp( - (ALGEBRAIC(:,33)+41.8000)./2.50000) - 1.00000) - ( 0.0849000.*(ALGEBRAIC(:,33)+6.80000))./(exp( - (ALGEBRAIC(:,33)+6.80000)./4.80000) - 1.00000);
    ALGEBRAIC(:,49) = (  - 0.0283900.*(ALGEBRAIC(:,33)+41.8000-CONSTANTS(:,108)))./(exp( - (ALGEBRAIC(:,33)+41.8000-CONSTANTS(:,108))./2.50000) - 1.00000) - ( 0.0849000.*(ALGEBRAIC(:,33)+6.80000-CONSTANTS(:,108)))./(exp( - (ALGEBRAIC(:,33)+6.80000-CONSTANTS(:,108))./4.80000) - 1.00000);
    ALGEBRAIC(:,55) = piecewise({ALGEBRAIC(:,10)== - 1.80000,  - 1.80001 }, ALGEBRAIC(:,10));
%     ALGEBRAIC(:,59) = ( 0.0114300.*(ALGEBRAIC(:,55)+1.80000))./(exp((ALGEBRAIC(:,55)+1.80000)./2.50000) - 1.00000);
    ALGEBRAIC(:,59) = ( 0.0114300.*(ALGEBRAIC(:,55)+1.80000-CONSTANTS(:,108)))./(exp((ALGEBRAIC(:,55)+1.80000-CONSTANTS(:,108))./2.50000) - 1.00000);
    ALGEBRAIC(:,61) = 0.00100000./(ALGEBRAIC(:,49)+ALGEBRAIC(:,59));
    ALGEBRAIC(:,19) = STATES(:,3);
    ALGEBRAIC(:,37) =  CONSTANTS(:,92).*log(CONSTANTS(:,12)./ALGEBRAIC(:,19));
    ALGEBRAIC(:,62) =  CONSTANTS(:,103).*CONSTANTS(:,24).*power(1.00000+power(CONSTANTS(:,22)./CONSTANTS(:,14), 1.20000),  - 1.00000).*power(1.00000+power(CONSTANTS(:,23)./ALGEBRAIC(:,19), 1.30000),  - 1.00000).*power(1.00000+exp( - ((ALGEBRAIC(:,10) - ALGEBRAIC(:,37))+110.000)./20.0000),  - 1.00000);
    ALGEBRAIC(:,64) = exp((  - CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92)));
    ALGEBRAIC(:,70) = 1.00000+ (CONSTANTS(:,15)./CONSTANTS(:,37)).*(1.00000+exp(( CONSTANTS(:,28).*ALGEBRAIC(:,10))./CONSTANTS(:,92)))+ (CONSTANTS(:,12)./CONSTANTS(:,35)).*(1.00000+ (CONSTANTS(:,12)./CONSTANTS(:,36)).*(1.00000+CONSTANTS(:,12)./CONSTANTS(:,34)));
    ALGEBRAIC(:,72) = ( (( (CONSTANTS(:,12)./CONSTANTS(:,35)).*CONSTANTS(:,12))./CONSTANTS(:,36)).*(1.00000+CONSTANTS(:,12)./CONSTANTS(:,34)).*exp((  - CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92))))./ALGEBRAIC(:,70);
    ALGEBRAIC(:,71) = ( (CONSTANTS(:,15)./CONSTANTS(:,37)).*exp(( CONSTANTS(:,28).*ALGEBRAIC(:,10))./CONSTANTS(:,92)))./ALGEBRAIC(:,70);
    ALGEBRAIC(:,68) = exp(( CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92)));
    ALGEBRAIC(:,63) = ALGEBRAIC(:,19)./(CONSTANTS(:,29)+ALGEBRAIC(:,19));
    ALGEBRAIC(:,73) =  ALGEBRAIC(:,64).*CONSTANTS(:,105).*(ALGEBRAIC(:,72)+ALGEBRAIC(:,71))+ ALGEBRAIC(:,71).*ALGEBRAIC(:,68).*(ALGEBRAIC(:,63)+ALGEBRAIC(:,64));
    ALGEBRAIC(:,65) = 1.00000+ (STATES(:,2)./CONSTANTS(:,30)).*(1.00000+exp((  - CONSTANTS(:,26).*ALGEBRAIC(:,10))./CONSTANTS(:,92))+ALGEBRAIC(:,19)./CONSTANTS(:,33))+ (ALGEBRAIC(:,19)./CONSTANTS(:,31)).*(1.00000+ (ALGEBRAIC(:,19)./CONSTANTS(:,32)).*(1.00000+ALGEBRAIC(:,19)./CONSTANTS(:,29)));
    ALGEBRAIC(:,66) = ( (STATES(:,2)./CONSTANTS(:,30)).*exp((  - CONSTANTS(:,26).*ALGEBRAIC(:,10))./CONSTANTS(:,92)))./ALGEBRAIC(:,65);
    ALGEBRAIC(:,67) = ( (( (ALGEBRAIC(:,19)./CONSTANTS(:,31)).*ALGEBRAIC(:,19))./CONSTANTS(:,32)).*(1.00000+ALGEBRAIC(:,19)./CONSTANTS(:,29)).*exp(( CONSTANTS(:,27).*ALGEBRAIC(:,10))./( 2.00000.*CONSTANTS(:,92))))./ALGEBRAIC(:,65);
    ALGEBRAIC(:,69) =  ALGEBRAIC(:,68).*ALGEBRAIC(:,63).*(ALGEBRAIC(:,67)+ALGEBRAIC(:,66))+ ALGEBRAIC(:,64).*ALGEBRAIC(:,66).*(CONSTANTS(:,105)+ALGEBRAIC(:,68));
    ALGEBRAIC(:,74) =  ALGEBRAIC(:,67).*ALGEBRAIC(:,63).*(ALGEBRAIC(:,72)+ALGEBRAIC(:,71))+ ALGEBRAIC(:,66).*ALGEBRAIC(:,72).*(ALGEBRAIC(:,63)+ALGEBRAIC(:,64));
    ALGEBRAIC(:,75) =  ALGEBRAIC(:,72).*CONSTANTS(:,105).*(ALGEBRAIC(:,67)+ALGEBRAIC(:,66))+ ALGEBRAIC(:,67).*ALGEBRAIC(:,71).*(CONSTANTS(:,105)+ALGEBRAIC(:,68));
    ALGEBRAIC(:,76) = ( (1.00000 - CONSTANTS(:,38)).*CONSTANTS(:,25).*( ALGEBRAIC(:,69).*ALGEBRAIC(:,71) -  ALGEBRAIC(:,73).*ALGEBRAIC(:,66)))./(ALGEBRAIC(:,73)+ALGEBRAIC(:,69)+ALGEBRAIC(:,74)+ALGEBRAIC(:,75));
    ALGEBRAIC(:,77) =  CONSTANTS(:,92).*log((CONSTANTS(:,12)+ 0.120000.*CONSTANTS(:,14))./(ALGEBRAIC(:,19)+ 0.120000.*CONSTANTS(:,13)));
%     ALGEBRAIC(:,78) =  CONSTANTS(:,39).*power(STATES(:,5), 3.00000).*STATES(:,6).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,77));
    ALGEBRAIC(:,78) =  CONSTANTS(:,128).*CONSTANTS(:,39).*power(STATES(:,5), 3.00000).*STATES(:,6).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,77));
    ALGEBRAIC(:,79) =  CONSTANTS(:,40).*power(STATES(:,5), 3.00000).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,77));
    ALGEBRAIC(:,80) = ALGEBRAIC(:,78)+ALGEBRAIC(:,79);
    ALGEBRAIC(:,50) =  STATES(:,4).*CONSTANTS(:,104).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,37)).*(1.00000 - CONSTANTS(:,20)).*CONSTANTS(:,121);
    ALGEBRAIC(:,83) =  (( 1.85000e-05.*CONSTANTS(:,42).*(ALGEBRAIC(:,10) - 0.00000))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))))).*(ALGEBRAIC(:,19) -  CONSTANTS(:,12).*exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))).*STATES(:,7).*STATES(:,8).*STATES(:,9);
    ALGEBRAIC(:,91) =  CONSTANTS(:,72).*STATES(:,2).*(1.00000 - STATES(:,21)) -  CONSTANTS(:,77).*STATES(:,21);
    ALGEBRAIC(:,85) =  (( 2.00000.*CONSTANTS(:,49).*ALGEBRAIC(:,10))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*ALGEBRAIC(:,10).*2.00000)./CONSTANTS(:,92))))).*(STATES(:,2) -  CONSTANTS(:,15).*exp((  - 2.00000.*ALGEBRAIC(:,10))./CONSTANTS(:,92))).*STATES(:,10).*STATES(:,11);
    ALGEBRAIC(:,81) =  (( 2.00000.*CONSTANTS(:,42).*(ALGEBRAIC(:,10) - 0.00000))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000).*2.00000)./CONSTANTS(:,92))))).*(STATES(:,2) -  CONSTANTS(:,15).*exp((  - 2.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))).*STATES(:,7).*STATES(:,8).*STATES(:,9);
    ALGEBRAIC(:,87) =  CONSTANTS(:,51).*STATES(:,13).*(STATES(:,16) - STATES(:,2));
    ALGEBRAIC(:,89) = (STATES(:,2) - STATES(:,18))./CONSTANTS(:,60);
    ALGEBRAIC(:,94) =  CONSTANTS(:,69).*STATES(:,18).*(1.00000 - STATES(:,22)) -  CONSTANTS(:,74).*STATES(:,22);
    ALGEBRAIC(:,92) = CONSTANTS(:,99)./(1.00000+exp(( - STATES(:,18)+CONSTANTS(:,63))./CONSTANTS(:,64)));
    ALGEBRAIC(:,95) = (STATES(:,17) - STATES(:,16))./CONSTANTS(:,61);
    ALGEBRAIC(:,97) =  CONSTANTS(:,71).*STATES(:,18).*(1.00000 - (STATES(:,23)+STATES(:,19))) -  CONSTANTS(:,75).*STATES(:,23);
    ALGEBRAIC(:,98) =  CONSTANTS(:,73).*STATES(:,16).*(1.00000 - STATES(:,24)) -  CONSTANTS(:,78).*STATES(:,24);
    ALGEBRAIC(:,100) =  CONSTANTS(:,72).*STATES(:,18).*(1.00000 - STATES(:,20)) -  CONSTANTS(:,77).*STATES(:,20);
    ALGEBRAIC(:,56) =  STATES(:,4).*CONSTANTS(:,102).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*(1.00000 - CONSTANTS(:,20)).*CONSTANTS(:,121);
    ALGEBRAIC(:,60) = ALGEBRAIC(:,50)+ALGEBRAIC(:,56);
%     ALGEBRAIC(:,90) =  CONSTANTS(:,88).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*( 0.900000.*STATES(:,30)+ 0.100000.*STATES(:,29)).*STATES(:,31);
    ALGEBRAIC(:,90) =  CONSTANTS(:,117).*CONSTANTS(:,88).*(ALGEBRAIC(:,10) - CONSTANTS(:,118) - CONSTANTS(:,97)).*( 0.900000.*STATES(:,30)+ 0.100000.*STATES(:,29)).*STATES(:,31);
    ALGEBRAIC(:,93) =  CONSTANTS(:,92).*log((CONSTANTS(:,14)+ 0.120000.*CONSTANTS(:,12))./(CONSTANTS(:,13)+ 0.120000.*ALGEBRAIC(:,19)));
%     ALGEBRAIC(:,96) =  CONSTANTS(:,95).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,93)).*power(STATES(:,32), 2.00000);
    ALGEBRAIC(:,96) =  CONSTANTS(:,132).*CONSTANTS(:,95).*(ALGEBRAIC(:,10) - ALGEBRAIC(:,93)).*power(STATES(:,32), 2.00000);
    ALGEBRAIC(:,88) =  CONSTANTS(:,87).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*STATES(:,27).*STATES(:,28);
    ALGEBRAIC(:,82) =  (( 0.000365000.*CONSTANTS(:,42).*(ALGEBRAIC(:,10) - 0.00000))./( CONSTANTS(:,92).*(1.00000 - exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))))).*(CONSTANTS(:,13) -  CONSTANTS(:,14).*exp((  - 1.00000.*(ALGEBRAIC(:,10) - 0.00000))./CONSTANTS(:,92))).*STATES(:,7).*STATES(:,8).*STATES(:,9);
    ALGEBRAIC(:,84) =  (ALGEBRAIC(:,81)+ALGEBRAIC(:,82)+ALGEBRAIC(:,83)).*(1.00000 - CONSTANTS(:,107)).*1.00000.*CONSTANTS(:,106);
    ALGEBRAIC(:,99) = piecewise({CONSTANTS(:,10)>0.00000,  CONSTANTS(:,91).*CONSTANTS(:,90).*(ALGEBRAIC(:,10) - CONSTANTS(:,97)).*(1.00000+exp((ALGEBRAIC(:,10)+20.0000)./20.0000)).*STATES(:,33) }, 0.00000);
    ALGEBRAIC(:,86) =  CONSTANTS(:,86).*STATES(:,25).*STATES(:,26).*(ALGEBRAIC(:,10) - CONSTANTS(:,97));
    ALGEBRAIC(:,101) = ALGEBRAIC(:,60)+ALGEBRAIC(:,90)+ALGEBRAIC(:,96)+ALGEBRAIC(:,88)+ALGEBRAIC(:,62)+ALGEBRAIC(:,76)+ALGEBRAIC(:,80)+ALGEBRAIC(:,84)+ALGEBRAIC(:,85)+ALGEBRAIC(:,99)+ALGEBRAIC(:,86);
    ALGEBRAIC(:,1) =  0.500000.*CONSTANTS(:,92).*log(CONSTANTS(:,15)./STATES(:,2));
    ALGEBRAIC(:,4) = STATES(:,16) - STATES(:,2);
    ALGEBRAIC(:,5) = STATES(:,12)+STATES(:,13)+STATES(:,14)+STATES(:,15);
    ALGEBRAIC(:,25) = (1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+23.2000)./6.60000)))./(0.846554./( 37.2000.*exp(ALGEBRAIC(:,10)./11.9000)+ 0.960000.*exp( - ALGEBRAIC(:,10)./18.5000)));
    ALGEBRAIC(:,26) =  4.00000.*(( 37.2000.*exp(ALGEBRAIC(:,10)./15.9000)+ 0.960000.*exp( - ALGEBRAIC(:,10)./22.5000))./0.846554 - (1.00000./(1.00000+exp( - (ALGEBRAIC(:,10)+23.2000)./10.6000)))./(0.846554./( 37.2000.*exp(ALGEBRAIC(:,10)./15.9000)+ 0.960000.*exp( - ALGEBRAIC(:,10)./22.5000))));
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end
