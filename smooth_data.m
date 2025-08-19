function smooth_data(infile, outfile, blockSize, method)
% SMOOTH_FILE_BY_MODE
%   smooth_file_by_mode(infile, outfile, blockSize, method)
%
%   infile    : path of the input text file (2 columns: time, value)
%   outfile   : path of the output text file
%   blockSize : (optional) block size, default = 50
%   method    : (optional) 'auto' (default), 'exact', or 'hist'
%
% Example:
%   smooth_file_by_mode('data.txt','data_smoothed.txt',50,'auto')

    if nargin < 3 || isempty(blockSize), blockSize = 50; end
    if nargin < 4 || isempty(method), method = 'auto'; end

    %--- Read data (expects 2 columns: time, variable) ---
    data = readmatrix(infile);
    if size(data,2) < 2
        error('The file must have at least two columns (time and value).');
    end
    t = data(:,1);   % time column
    x = data(:,2);   % variable column

    %--- Preallocate result ---
    n = numel(x);
    y = x; % will be replaced by smoothed values

    %--- Process in blocks ---
    for i = 1:blockSize:n
        idx = i:min(i+blockSize-1, n);
        xi = x(idx);

        % Ignore NaNs in this block
        xi_valid = xi(~isnan(xi));

        if isempty(xi_valid)
            % If all NaNs, leave as is
            continue;
        end

        % Choose strategy to compute "mode"
        blockMode = compute_block_mode(xi_valid, method);

        % Replace all values in this block with the "mode"
        y(idx) = blockMode;
    end

    %--- Save output: [time, smoothed values] ---
    out = [t, y];
    writematrix(out, outfile,'Delimiter',' ');
    fprintf('Smoothed file saved to: %s\n', outfile);
end

function m = compute_block_mode(xv, method)
% Compute the representative value ("mode") for a block

    switch lower(method)
        case 'exact'
            % Strict mode (for repeated discrete values)
            m = mode(xv);

        case 'hist'
            % Estimate mode using histogram (for continuous values)
            m = hist_mode(xv);

        case 'auto'
            % Decide automatically: if many repeated values -> mode,
            % otherwise -> histogram mode
            uq = unique(xv);
            if numel(uq) <= max(5, round(numel(xv)/5))
                m = mode(xv);
            else
                m = hist_mode(xv);
            end

        otherwise
            error('Unknown method. Use ''auto'', ''exact'' or ''hist''.');
    end
end

function m = hist_mode(xv)
% Estimate the mode for continuous data:
% - Compute bin width with Freedman–Diaconis rule
% - If invalid, fallback to Scott’s rule or mean

    n = numel(xv);
    if n == 1
        m = xv; return;
    end

    % Freedman–Diaconis rule
    iq = iqr(xv);
    bw = 2*iq / (n^(1/3));

    % Fallback if FD fails
    if ~isfinite(bw) || bw <= 0
        s = std(xv);
        bw = 3.5*s / (n^(1/3));
    end
    if ~isfinite(bw) || bw <= 0
        % All values equal
        m = xv(1);
        return;
    end

    xmin = min(xv);
    xmax = max(xv);
    if xmin == xmax
        m = xmin; return;
    end

    % Build histogram edges
    edges = xmin:bw:xmax;
    if numel(edges) < 2
        % Too few bins -> fallback to mean
        m = mean(xv);
        return;
    end

    % Histogram counts
    counts = histcounts(xv, edges);
    [~, k] = max(counts);

    % Center of most populated bin
    left = edges(k);
    right = edges(k+1);
    m = (left + right)/2;
end

%% Default block size (50), automatic method:
% smooth_data('input.txt','output.txt');
% 
%% Force strict mode (good for discrete/repeated values):
% smooth_data('input.txt','output_exact.txt',50,'exact');
% 
%% Force histogram mode (good for continuous/noisy data):
% smooth_data('input.txt','output_hist.txt',50,'hist');