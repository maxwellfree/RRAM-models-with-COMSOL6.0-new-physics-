function [con]= GetReset(n,Ea,Rtop,Rbottom,Amplitude,threeD)
% FindResetTrans
% --------------------------------------------------------------
% This function performs repeated short-time simulations on a COMSOL model
% to capture transient changes between solution states.
%
% The model represents a resistive switching cell with a conductive 
% filament.
% It tracks:
%   - Maximum filament temperature over time (Tmax)
%   - Voltage evolution across the cell (Vme)
%   - Current evolution through the cell (Ime)
%   - When filament rupture occurs (due to a sequence of electrical pulses)
%
% Inputs:
%   n         - Number of simulation time steps to run
%   Ea        - Activation energy [eV]
%   Rtop      - Filament radius at the top electrode [nm]
%   Rbottom   - Filament radius at the bottom electrode [nm]
%   Amplitude - Pulse amplitude [nm] (applied change to model parameter)
%   threeD    - Boolean flag to export 3D temperature field
%
% Outputs:
%   con       - Current simulation step counter when finished
% --------------------------------------------------------------

% Path where exported COMSOL data files will be saved
pathCal = '/home/moreno/Desktop/TODO/Transitorio/cal/';

% Load the COMSOL model file
model=mphload('Ni_HfO2_Si.mph');

% Set physical and geometrical parameters in the COMSOL model
model.param.set('Ea',strcat(num2str(Ea),'[eV]') , 'Activation energy');
model.param.set('Rtop',strcat(num2str(Rtop),'[nm]'));
model.param.set('Rbottom',strcat(num2str(Rbottom),'[nm]'));
model.param.set('Amplitude',strcat(num2str(Amplitude),'[nm]'));

% Initial export of right filament shape from the theoretical truncated 
% cone (step 0)
model.result.export('plot5').set('filename', strcat(pathCal, ...
    'profileIvsTime',num2str(0),'.txt'));
model.result.export('plot5').set('header', false);
model.result.export('plot5').run;

% Load this initial profile into COMSOL function 'int4'
model.func('int4').set('filename', strcat(pathCal, ...
    'profileIvsTime',num2str(0),'.txt'));
model.func('int4').refresh;

% Initialization of loop variables
logicalreset = true;  % Will remain true until rupture condition or n steps 
                      % reached.

con   = 1;            % Step counter.
time  = 0;            % Current simulation time.
timeSample = 5.0e-8;  % Time increment per step [s].

% =========================================================================
% Main transient simulation loop
% Each iteration simulates a small time step, updates parameters,
% and exports results. Stops when filament ruptures or n steps reached.
% =========================================================================
while (logicalreset)

    % Set simulation start and end times for this time slice
    model.param.set('TimeStart',strcat(num2str(time),'[s]'), ...
        'Inicio de tiempo');

    time = time + timeSample;

    model.param.set('TimeMax',strcat(num2str(time),'[s]'), ...
        'Inicio de tiempo');

    % Update COMSOL interpolation functions with latest exported data. We 
    % will use the previous profiles as a geometric starting point to 
    % generate the new shapes the filament will have. The deformation the 
    % filament will undergo will depend on the differential equation that 
    % governs its shape, where the roles in matter subtraction depend 
    % strongly on the temperature and activation energy.
    model.func('int4').set('filename',strcat(pathCal,'profileIvsTime', ...
        num2str(con-1),'.txt'));
    model.func('int4').refresh;

    % The following code statement is simply information about the 
    % temperature-file being exported
    strcat(pathCal,'T_IvsV', num2str(con-1),'.txt')
    % The solution derived from the final geometry for the temperature 
    % field, from the axisymmetric cut of revolution which is the simulation
    % domain, is used as a starting point or initial condition to carry out
    % the simulation of the subsequent time period.
    model.func('int2').set('filename',strcat(pathCal,'T_IvsV', ...
        num2str(con-1),'.txt'));
    model.func('int2').refresh;

    % The following code statement is simply information about the 
    % voltage-file being exported
    strcat(pathCal,'V_IvsV', num2str(con-1),'.txt')
    % The solution derived from the final geometry for the electric 
    % potential field, from the axisymmetric section of revolution which 
    % is the simulation domain, is used as a starting point or initial 
    % condition to carry out the simulation of the subsequent time period.
    model.func('int3').set('filename',strcat(pathCal,'V_IvsV', ...
        num2str(con-1),'.txt'));
    model.func('int3').refresh;

    % Display model parameters and simulation progress
    mphgetexpressions(model.param)
    strcat('Progress...',num2str((con/n)*100),'%')

    % -----------------------------------------------------------
    % We use MATLAB's 'try' statement here to handle the situation
    % where the filament becomes so thin that no current can pass
    % through it.
    %
    % Physically, at this point the geometry has evolved so that
    % only dielectric material remains between the electrodes,
    % meaning the reset process is complete.
    %
    % When this occurs, COMSOL will not be able to achieve
    % convergence for the next time step because there is no
    % conductive path. Without 'try', this would cause MATLAB to
    % throw an error and stop abruptly.
    %
    % By wrapping the simulation step in 'try', we can detect this
    % condition, interrupt the simulation gracefully, and exit the
    % loop even if not all scheduled time steps have been executed.
    % -----------------------------------------------------------
     try
        % Run COMSOL study for this time slice 
        model.study('std1').run;

        % Export updated profiles for filament geometry, voltage, temperature
        model.result.export('plot2').set('filename', strcat(pathCal, ...
            'profileIvsTime',num2str(con),'.txt'));
        model.result.export('plot2').set('header', false);
        model.result.export('plot2').run;

        model.result.export('plot3').set('filename', strcat(pathCal, ...
            'V_IvsV',num2str(con),'.txt'));
        model.result.export('plot3').set('header', false);
        model.result.export('plot3').run;

        model.result.export('plot4').set('filename', strcat(pathCal, ...
            'T_IvsV',num2str(con),'.txt'));
        model.result.export('plot4').set('header', false);
        model.result.export('plot4').run;  

        % Stop simulation after n steps
        if (con>n)
            logicalreset = false;
        end

     catch
        % If COMSOL throws an error, terminate loop. This error may be 
        % caused primarily because the reset has been reached, i.e. the 
        % filament has been completely broken and therefore it is impossible 
        % to reach numerical convergence. 
        logicalreset = false;
        as='ERROR'
    end

     % Increment step counter. Each step taken increases the time by 5.0e-8 
     % seconds. 
    con = con + 1;

    % Collect measurement results if simulation is still running
    if (logicalreset)
        % Maximum temperature from table 'tbl1'. What we do is take the 
        % domain corresponding to the entire filament and find its maximum 
        % temperature. This is done over time, taking ALL temporal 
        % states of each simulation period. By taking this maximum 
        % temperature, we can see how the filament evolves throughout 
        % the transient process that modifies its geometry during the 
        % reset process.
        model.result.table('tbl1').clearTableData;
        model.result.numerical('max1').setResult;
        model.result.export('tbl1').set('filename', strcat(pathCal, ...
            'current/timeTMax.txt'));
        model.result.export('tbl1').set('header', false);
        model.result.export('tbl1').set('ifexists', 'append');
        model.result.export('tbl1').run;

        % Average voltage from table 'tbl3'. What is done is to determine 
        % the average voltage that is right in the metallization that is 
        % being subjected to the source. This is done over time, taking 
        % ALL temporal states of each simulation period.
        model.result.table('tbl3').clearTableData;
        model.result.numerical('av1').setResult;
        model.result.export('tbl3').set('filename', strcat(pathCal, ...
            'current/timeVapp.txt'));
        model.result.export('tbl3').set('header', false);
        model.result.export('tbl3').set('ifexists', 'append');
        model.result.export('tbl3').run;

        % Integrated current from numerical feature 'int1'. To see what 
        % current is passing through the cell and to see the coherence 
        % in the reset process (we hope to reduce this current) we will 
        % record the integral of the current density over the grounded 
        % metallization to see the total intensity of current passing 
        % through the cell for each voltage pulse supplied to it. This is 
        % done over time, taking ALL temporal states of each simulation 
        % period.
        model.result.table('tbl2').clearTableData;
        model.result.numerical('int1').setResult;
        model.result.export('tbl2').set('filename', strcat(pathCal, ...
            'current/timeIntJz.txt'));
        model.result.export('tbl2').set('header', false);
        model.result.export('tbl2').set('ifexists', 'append');
        model.result.export('tbl2').run;

        % Optional: export full 3D temperature field if requested
        if (threeD) 
            model.result.export('data1').set('filename', strcat(pathCal, ...
                '3D/T3D',num2str(con),'.txt'));
            model.result.export('data1').set('location', 'grid');
            model.result.export('data1').set('gridx3', 'range(-Rtop-2,0.25,Rtop+2)');
            model.result.export('data1').set('gridy3', 'range(-Rtop-2,0.25,Rtop+2)');
            model.result.export('data1').set('gridz3', 'range(zbottom-2,0.25,ztop+2)');
            model.result.export('data1').set('header', false);
            model.result.export('data1').run;
        end  
    
    end % End of logicalreset

end  % End of while loop

end % End of function

%% Example call:
% [con]= GetReset(990,1.08,4,1,1.0e7,false)
