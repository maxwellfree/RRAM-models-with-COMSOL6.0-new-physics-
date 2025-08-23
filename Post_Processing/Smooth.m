function [outputInfo] = Smooth(scaleGain, path, inFile, outFile, opts)
% SMOOTH  Wrapper that reads a 2-column text file (time, signal),
%         applies pulse-trend smoothing, and writes a space-separated
%         text file with (time, smoothed_signal).
%
%   [outputInfo] = Smooth(scaleGain, inFile, outFile, opts)
%
% Inputs
%   scaleGain : numeric gain applied to the input signal before smoothing
%               (use 1 if you do not want to scale). Another very attractive
%               option is for this value to be -1. This has the advantage 
%               of inverting the data. Let's imagine that the current, for 
%               example, has a negative value but we want to express all 
%               positive data. By changing this to minus one, we can invert 
%               the data.
%   inFile    : input filename (relative to 'path' or absolute if preferred)
%   outFile   : output filename (relative to 'path' or absolute if preferred)
%   opts      : (optional) struct passed directly to SMOOTH_PULSE_TREND
%               The most relevant fields are:
%                 .medianFrac   (default 0.03)
%                 .sgolayFrac   (default 0.12)
%                 .sgolayOrder  (default 3)
%                 .minWin       (default 5)
%                 .maxWin       (default 501)
%                 .periodGuess  (seconds; override auto period estimation)
%                 .plot         (true/false; default true)
%
% Output
%   outputInfo : diagnostics from SMOOTH_PULSE_TREND (fs, dt, N, Tdom, windows, ...)
%
% Notes
%   - This function expects a variable 'path' in the workspace that points
%     to the directory where files are read/written. If you prefer, replace
%     'fullfile(path, ...)' by absolute paths.
%   - The input file must have exactly two numeric columns: [time value].
%   - The output is written as plain text with a single space delimiter.

    if nargin < 4, opts = struct; end
    if nargin < 1 || isempty(scaleGain), scaleGain = 1; end

    % Build absolute/relative paths for I/O
    infile  = strcat(path, inFile);
    outfile = strcat(path, outFile);

    % Delegate the heavy lifting to the inner helper
    outputInfo = smooth_dataII(infile, outfile, scaleGain, opts);
end

function outputInfo = smooth_dataII(infile, outfile, scaleGain, opts)
% Helper that:
%   1) reads [time value] from 'infile'
%   2) scales 'value' by 'scaleGain'
%   3) calls smooth_pulse_trend(t, x, opts)
%   4) writes [time smoothed] into 'outfile' as space-separated text

    % --- Validate input file exists ---
    if ~isfile(infile)
        error('Input file does not exist: %s', infile);
    end

    % --- Read data (expects 2 columns: time, value) ---
    data = readmatrix(infile);
    if size(data,2) < 2
        error('The file must have at least two numeric columns: [time value].');
    end

    % Ensure column vectors
    t = data(:,1);
    x = data(:,2);

    % Optional scaling (sign/gain). Use scaleGain=1 to keep original scale.
    x = x .* scaleGain;

    try
        % --- Pulse-trend smoothing (median + Savitzky–Golay) ---
        %     opts is forwarded so you can control windows/order/plot/periodGuess, etc.
        [y_smooth, info] = smooth_pulse_trend(t, x, opts);
        outputInfo = info;
    catch ME
        % Graceful fallback while preserving context
        warning('Smoothing failed: %s', ME.message);
        outputInfo = struct('error', ME.message);
        % In case of failure, return original signal to keep output consistent
        y_smooth = x;
    end

    % --- Write output: two columns (time, smoothed), space-separated text ---
    out = [t(:) y_smooth(:)];
    writematrix(out, outfile, 'Delimiter', ' ', 'FileType', 'text');

    fprintf('Smoothed file saved to: %s\n', outfile);
end

% ===============================================================
% ----------------------
% OPTIONS (struct 'opts')
% ----------------------
%   .medianFrac   : Median filter window as fraction of dominant period
%                   (default = 0.03). Larger -> more robust spike removal.
%
%   .sgolayFrac   : Savitzky–Golay window as fraction of dominant period
%                   (default = 0.12). Larger -> stronger smoothing, but may
%                   flatten sharp features.
%
%   .sgolayOrder  : Polynomial order for Savitzky–Golay filter (default = 3).
%                   Higher order preserves curvature but can fit noise.
%
%   .minWin       : Minimum window size in samples (default = 5).
%   .maxWin       : Maximum window size in samples (default = 501).
%
%   .periodGuess  : If known, manually set dominant period (seconds).
%                   Overrides automatic estimation for more stable windows.
%
%   .plot         : Logical flag to show raw vs smoothed plot (default = true).
%
% -----------------
% DEFAULT USAGE
% -----------------
% path      = '/your/path/';
% inFile    = 'time_VS_I.txt';
% outFile   = 'tI1.txt';
% opts.plot = false;  % skip plotting
% info = Smooth(1, path, inFile, outFile, opts);
%   -> uses default opts, gain=1 (no scaling), with automatic period detection.
%
% -----------------
% EXAMPLES
% -----------------
% Example 0:
% path='/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s/current/';
% inFile='time_VS_I.txt';
% outFile='tI1.txt';
% opts.plot = false; 
% info = Smooth(-1, path, inFile, outFile, opts)
% Example 1: Stronger smoothing (wider windows), no plotting:
% opts.medianFrac  = 0.05;   % more aggressive spike removal
% opts.sgolayFrac  = 0.18;   % smoother curve
% opts.plot        = false;  % skip plotting
% info = Smooth(1, 'time_VS_I.txt', 'tI1_smooth.txt', opts);
%
% Example 2: Known pulse period (1e-6 s) and signal inversion:
% opts.periodGuess = 1e-6;   % known dominant period
% info = Smooth(-1, 'time_VS_I.txt', 'tI1_inverted.txt', opts);
%
% ===============================================================