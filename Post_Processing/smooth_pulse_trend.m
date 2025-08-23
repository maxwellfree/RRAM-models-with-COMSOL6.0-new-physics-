function [y_smooth, info] = smooth_pulse_trend(t, x, opts)
% SMOOTH_PULSE_TREND  Denoise a pulsed signal while preserving its trend/shape.
%   [y_smooth, info] = smooth_pulse_trend(t, x, opts)
%
% Inputs
%   t    : time vector (monotonic; not necessarily uniform but recommended)
%   x    : signal vector (temperature, voltage, current, etc.)
%   opts : (optional) struct with fields:
%          .medianFrac   : median filter window as fraction of period (default 0.03)
%          .sgolayFrac   : sgolay window as fraction of period (default 0.12)
%          .sgolayOrder  : polynomial order for sgolay (default 3)
%          .minWin       : minimum window in samples (default 5)
%          .maxWin       : maximum window in samples (default 501)
%          .periodGuess  : manual period (in seconds) if you want to override auto
%          .plot         : logical, plot original vs smoothed (default true)
%
% Outputs
%   y_smooth : smoothed/denoised signal preserving pulse trend
%   info     : struct with diagnostics (estimated period, window sizes, etc.)
%
% Method
%   1) Spike/noise suppression via median filter (robust to outliers).
%   2) Edge-preserving smoothing via Savitzky–Golay.
%   Windows are auto-sized from the dominant period estimated by autocorrelation.
%
% Notes
%   - If your data are very irregularly sampled, consider resampling to a
%     near-uniform grid first for best results.
%   - For purely noise-like high-frequency ripples on top of broad pulses,
%     sgolay is a good compromise between smoothing and edge preservation.

    % ----------- Input checks -----------
    if nargin < 3, opts = struct; end
    mustBeVector(t); mustBeVector(x);
    t = t(:); x = x(:);
    if numel(t) ~= numel(x)
        error('t and x must have the same length.');
    end

    % Sort by time if needed
    [t, sortIdx] = sort(t, 'ascend');
    x = x(sortIdx);

    % Defaults
    medianFrac  = getfielddef(opts, 'medianFrac', 0.03);
    sgolayFrac  = getfielddef(opts, 'sgolayFrac', 0.12);
    sgOrder     = getfielddef(opts, 'sgolayOrder', 3);
    minWin      = getfielddef(opts, 'minWin', 5);
    maxWin      = getfielddef(opts, 'maxWin', 501);
    doPlot      = getfielddef(opts, 'plot', true);

    % Estimate sampling interval (robust)
    dt  = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        error('Non-positive or invalid time step detected.');
    end
    fs  = 1/dt;

    % ----------- Estimate dominant period -----------
    if isfield(opts, 'periodGuess') && ~isempty(opts.periodGuess) && opts.periodGuess > 0
        Tdom = opts.periodGuess;
    else
        Tdom = estimate_dominant_period(t, x);
        if ~isfinite(Tdom)
            % Fallback if estimation fails: use 10% of signal duration as "period"
            Tdom = 0.1*(t(end)-t(1));
        end
    end
    Nper = max(1, round(Tdom*fs));  % samples per dominant period

    % ----------- Build windows (odd lengths for filters) -----------
    % Median filter: short window for spike suppression
    Nmed = clamp_odd(round(max(minWin, medianFrac * Nper)), minWin, maxWin);
    % Savitzky–Golay: wider window to smooth ripples but preserve edges
    Nsg  = clamp_odd(round(max(minWin, sgolayFrac * Nper)), max(Nmed+2, minWin), maxWin);
    % Ensure sgolay window > polynomial order
    if Nsg <= sgOrder
        Nsg = sgOrder + 2 + mod(sgOrder,2); % make it odd and > order
    end

    % ----------- Filtering pipeline -----------
    % 1) Robust spike removal
    x_med = medfilt1(x, Nmed, 'truncate');  % 'truncate' keeps window within bounds
    % 2) Edge-preserving smoothing
    y_smooth = sgolayfilt(x_med, sgOrder, Nsg);

    % ----------- Output info -----------
    info.fs       = fs;
    info.dt       = dt;
    info.N        = numel(x);
    info.Tdom     = Tdom;
    info.Nper     = Nper;
    info.Nmed     = Nmed;
    info.Nsg      = Nsg;
    info.sgOrder  = sgOrder;

    % ----------- Optional plot -----------
    if doPlot
        figure; 
        plot(t, x, '.', 'DisplayName','raw'); hold on;
        plot(t, y_smooth, '-', 'LineWidth', 1.5, 'DisplayName','smoothed');
        grid on; xlabel('Time'); ylabel('Signal');
        title(sprintf('Pulse-trend smoothing (N_{med}=%d, N_{sg}=%d, order=%d)', Nmed, Nsg, sgOrder));
        legend('Location','best');
    end
end

% ----------------- Helpers -----------------
function val = getfielddef(s, name, def)
    if isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end

function o = clamp_odd(n, lo, hi)
    n = max(lo, min(hi, n));
    if mod(n,2)==0, n = n+1; end
    o = n;
end

function T = estimate_dominant_period(t, x)
% Estimate dominant period via autocorrelation peak on a uniformly resampled grid.
    % Uniform resample to handle mild jitter
    N   = numel(x);
    dt  = median(diff(t));
    tu  = (t(1):dt:t(end)).';
    xu  = interp1(t, x, tu, 'linear','extrap');

    % Remove mean to emphasize oscillatory structure
    xu = xu - mean(xu);

    % Autocorrelation (biased)
    [ac, lags] = xcorr(xu, 'biased');
    % Keep positive lags only
    ac   = ac(lags >= 0);
    lags = lags(lags >= 0);
    if numel(ac) < 5
        T = NaN; return;
    end

    % Find first significant local maximum beyond lag 0
    % Use a conservative minimum separation to avoid the zero-lag peak
    minSep = max(3, round(0.02*numel(xu))); % 2% of length
    [pks, locs] = findpeaks(ac, 'MinPeakDistance', minSep);
    if isempty(locs)
        T = NaN; return;
    end
    % Choose the strongest peak past lag 0
    [~, idx] = max(pks);
    lagSamples = locs(idx)-1; % adjust because we kept non-negative lags
    if lagSamples <= 0
        T = NaN; return;
    end
    T = lagSamples * dt;
end