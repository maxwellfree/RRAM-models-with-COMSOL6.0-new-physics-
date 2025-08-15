%% Save current working directory so we can return after changing folders
Currentdir=pwd;

%% Change into a project directory (where your transient/“Transitorio” 
%% files live)
cd(['/home/moreno/Desktop/TODO/Transitorio' ...
    '']);  % The empty string is redundant but harmless

%% Start the COMSOL Multiphysics server as a background process on Linux
% - This server is what LiveLink for MATLAB connects to in order to run 
%   COMSOL models.
% - The equivalent command if typed directly into a Linux terminal is:
%       comsol server -port 2026
%   This launches the COMSOL server process, listening on TCP port 2026.
% - In MATLAB we can launch it via system(...). 
% - The & puts it in the background so MATLAB can continue execution.
system('comsolmphserver &');

%% Go back to the original working directory
cd(Currentdir);

%% (Again) save current working directory (optional redundancy)
Currentdir=pwd;
%% Change into COMSOL LiveLink for MATLAB interface directory
% This folder contains mphstart.m and supporting Java classes.
cd('/usr/local/comsol62/multiphysics/mli');

%% Start a LiveLink session and connect MATLAB to the COMSOL server
% - mphstart(2026) tells MATLAB to connect to COMSOL Server running on 
%   port 2026.
% - This must match the port number used when starting the server 
%   (comsol server -port 2026).
mphstart(2026);

%% Return to original directory
cd(Currentdir);

%% Import Java classes to access the COMSOL Model API from MATLAB
% com.comsol.model.*       → Core model classes (geometry, physics, study, 
% etc.)
% com.comsol.model.util.*  → Utility classes (model management, progress, 
% licensing)
import com.comsol.model.*;
import com.comsol.model.util.*

%% Show current model tags (IDs) in the COMSOL session
% Useful for checking which models are loaded.
ModelUtil.tags

%% Remove the default model named 'Model' if it exists
% This clears the workspace so you can load or create new models without 
% conflicts.
ModelUtil.remove('Model')

%% Enable progress display in COMSOL during operations
% When true, COMSOL will show progress windows/messages for long tasks.
ModelUtil.showProgress(true)