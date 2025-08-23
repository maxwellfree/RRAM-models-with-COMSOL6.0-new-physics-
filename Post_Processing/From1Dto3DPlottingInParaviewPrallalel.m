function From1Dto3DPlottingInParaviewPrallalel(inputDir, filePattern, outputDir, opts)
% GENERATE_FRUSTUM_VOLUME_SERIES (parallel + progress)
% Build 3D binary volumes (inside=1, outside=0) by revolving 2D profiles R(z)
% around the Z-axis, and export one VTK per time step.
%
% Parallelization:
%   - Each input file (time step) is processed on a separate worker (parfor).
%   - Uses a DataQueue + waitbar to show progress safely from parfor.
%
% Options (opts):
%   .Nx, .Ny, .Nz   : grid size (default 128 x 128 x 128)
%   .zPad, .rPad    : domain padding (default 0)
%   .zRange         : fixed [zmin zmax] (optional)
%   .rMax           : fixed max radius (optional)
%   .insideVal      : default 1 (uint8)
%   .outsideVal     : default 0 (uint8)
%   .ascii          : true (default)
%   .useParallel    : true/false (default true if PCT available)
%   .numWorkers     : desired workers (default = min( feature('numcores'), 16 ))
%   .showProgress   : true/false (default true)
%
% Notes:
%   - All workers write different output files -> no contention.
%   - Precomputed grid/rho are broadcast to workers (ensure memory fits Nx*Ny).

    if nargin < 4, opts = struct; end
    Nx = getOpt(opts,'Nx',128);
    Ny = getOpt(opts,'Ny',128);
    Nz = getOpt(opts,'Nz',128);
    zPad = getOpt(opts,'zPad',0);
    rPad = getOpt(opts,'rPad',0);
    insideVal  = uint8(getOpt(opts,'insideVal',1));
    outsideVal = uint8(getOpt(opts,'outsideVal',0));
    writeAscii = getOpt(opts,'ascii',true);
    showProgress = getOpt(opts,'showProgress', true);

    % Parallel settings
    havePCT     = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    defaultNW   = min( feature('numcores'), 16 );
    useParallel = getOpt(opts,'useParallel', havePCT);
    numWorkers  = getOpt(opts,'numWorkers', defaultNW);

    if ~isfolder(outputDir), mkdir(outputDir); end

    % ---- List & order files by trailing number ----
    files = dir(fullfile(inputDir, filePattern));
    if isempty(files), error('No files found for pattern: %s', filePattern); end
    [files, idxNums] = sortFilesByTrailingNumber(files); %#ok<ASGLU>
    Nfiles = numel(files);

    % ---- Global bounds (z-range & rmax) ----
    if isfield(opts,'zRange') && numel(opts.zRange)==2
        zmin = opts.zRange(1); zmax = opts.zRange(2);
    else
        zAllMin = +inf; zAllMax = -inf;
        for k = 1:Nfiles
            A = readmatrix(fullfile(files(k).folder, files(k).name));
            Z = A(:,1);
            zAllMin = min(zAllMin, min(Z));
            zAllMax = max(zAllMax, max(Z));
        end
        zmin = zAllMin - zPad;
        zmax = zAllMax + zPad;
    end

    if isfield(opts,'rMax')
        rmax = opts.rMax;
    else
        rAllMax = 0;
        for k = 1:Nfiles
            A = readmatrix(fullfile(files(k).folder, files(k).name));
            R = abs(A(:,2));
            rAllMax = max(rAllMax, max(R));
        end
        rmax = rAllMax + rPad;
    end
    if rmax <= 0, error('Non-positive rmax computed. Check input data.'); end

    % ---- Fixed Cartesian grid & rho(x,y) (broadcasted to workers) ----
    x = linspace(-rmax, +rmax, Nx);
    y = linspace(-rmax, +rmax, Ny);
    z = linspace(  zmin,   zmax, Nz);
    dx = (x(end)-x(1)) / max(1,Nx-1);
    dy = (y(end)-y(1)) / max(1,Ny-1);
    dz = (z(end)-z(1)) / max(1,Nz-1);
    [X,Y] = ndgrid(x,y);
    rhoXY = sqrt(X.^2 + Y.^2);  % Nx×Ny

    % ---- Progress indicator setup ----
    h = [];
    cleanupH = onCleanup(@() closeWaitbarIfOpen());
    nDone = 0;

    if showProgress
        % If headless (no display), waitbar may fail -> use textual fallback
        try
            h = waitbar(0, sprintf('Processing 0/%d (0%%)', Nfiles), 'Name','Frustum series');
        catch
            h = []; fprintf('Progress: 0/%d (0%%)\n', Nfiles);
        end
    end

    if useParallel
        % start or resize pool
        p = gcp('nocreate');
        if isempty(p)
            parpool('local', numWorkers);
        elseif p.NumWorkers ~= numWorkers
            delete(p);
            parpool('local', numWorkers);
        end

        % DataQueue to update progress safely from workers
        q = parallel.pool.DataQueue;
        afterEach(q, @updateProgress);

        % Parallel loop over files
        parfor k = 1:Nfiles
            processOneFile(k, files, idxNums, z, rhoXY, insideVal, outsideVal, ...
                           writeAscii, [Nx Ny Nz], [x(1) y(1) z(1)], [dx dy dz], outputDir);
            send(q, 1); % notify progress
        end
    else
        % Sequential fallback with progress updates
        for k = 1:Nfiles
            processOneFile(k, files, idxNums, z, rhoXY, insideVal, outsideVal, ...
                           writeAscii, [Nx Ny Nz], [x(1) y(1) z(1)], [dx dy dz], outputDir);
            updateProgress(); % local update
        end
    end

    if ~isempty(h) && isvalid(h)
        waitbar(1, h, sprintf('Processing %d/%d (100%%)', Nfiles, Nfiles));
    end
    fprintf('Done. Wrote %d VTK files to: %s\n', Nfiles, outputDir);

    % ----- nested helpers for progress -----
    function updateProgress(~)
        nDone = nDone + 1;
        if showProgress
            pct = nDone / Nfiles;
            msg = sprintf('Processing %d/%d (%.0f%%)', nDone, Nfiles, 100*pct);
            if ~isempty(h) && isvalid(h)
                try
                    waitbar(pct, h, msg);
                catch
                    % fallback textual
                    fprintf('%s\n', msg);
                end
            else
                fprintf('%s\n', msg);
            end
        end
    end

    function closeWaitbarIfOpen()
        if ~isempty(h) && isvalid(h), close(h); end
    end
end


% ---------- Per-file worker task ----------
function processOneFile(k, files, idxNums, z, rhoXY, insideVal, outsideVal, ...
                        writeAscii, dims, origin, spacing, outputDir)

    filePath = fullfile(files(k).folder, files(k).name);
    A = readmatrix(filePath);
    if size(A,2) < 2
        error('File %s must have at least two columns [Z R].', files(k).name);
    end

    Z = A(:,1); R = A(:,2);
    [Z, iu] = unique(Z(:),'stable');
    R = abs(R(:)); R = R(iu);          % non-negative radius

    % Interpolate R(z) on common z-grid. Outside -> 0
    Rz = interp1(Z, R, z, 'linear', NaN);
    Rz(~isfinite(Rz)) = 0;
    Rz = max(0, Rz);

    % Build volume: inside if rho <= R(z)
    vol = rhoXY <= reshape(Rz, 1,1,[]);
    V = uint8(vol);
    V(~vol) = outsideVal;
    V(vol)  = insideVal;

    % Output name by numeric suffix if present; else k-1
    outIndex = idxNums(k);
    if isnan(outIndex), outIndex = k-1; end
    outName = sprintf('frustum_%04d.vtk', outIndex);
    outPath = fullfile(outputDir, outName);

    writeVTKStructuredPoints(outPath, V, dims, origin, spacing, writeAscii);
end


% ---------- Small utilities ----------
function val = getOpt(s, f, def)
    if isfield(s,f) && ~isempty(s.(f)), val = s.(f); else, val = def; end
end

function [filesSorted, idxNums] = sortFilesByTrailingNumber(files)
    idxNums = nan(numel(files),1);
    for i = 1:numel(files)
        nm = files(i).name;
        tok = regexp(nm, '(\d+)(?=\.[^.]+$)', 'tokens', 'once');
        if ~isempty(tok)
            idxNums(i) = str2double(tok{1});
        end
    end
    hasNum = ~isnan(idxNums);
    [~, orderNum] = sort(idxNums(hasNum));
    [~, orderStr] = sort({files(~hasNum).name});
    filesSorted = [files(hasNum); files(~hasNum)];
    filesSorted(hasNum) = files(hasNum(orderNum));
    filesSorted(~hasNum) = files(~hasNum(orderStr));
    idxNums = [idxNums(hasNum); idxNums(~hasNum)];
    idxNums(hasNum) = idxNums(orderNum);
    idxNums(~hasNum) = idxNums(orderStr);
end

function writeVTKStructuredPoints(filename, V, dims, origin, spacing, writeAscii)
    Nx=dims(1); Ny=dims(2); Nz=dims(3);
    x0=origin(1); y0=origin(2); z0=origin(3);
    dx=spacing(1); dy=spacing(2); dz=spacing(3);

    fid = fopen(filename,'w');
    if fid==-1, error('Cannot open %s for writing.', filename); end

    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Frustum volume\n');
    if writeAscii
        fprintf(fid, 'ASCII\n');
    else
        error('Binary write not implemented in this snippet.');
    end
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS %d %d %d\n', Nx, Ny, Nz);
    fprintf(fid, 'ORIGIN %.9g %.9g %.9g\n', x0, y0, z0);
    fprintf(fid, 'SPACING %.9g %.9g %.9g\n', dx, dy, dz);
    fprintf(fid, 'POINT_DATA %d\n', Nx*Ny*Nz);
    fprintf(fid, 'SCALARS frustum unsigned_char 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');

    % x-fastest order
    cnt = 0;
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                fprintf(fid, '%d ', V(i,j,k));
                cnt = cnt + 1;
                if mod(cnt, 20)==0, fprintf(fid, '\n'); end
            end
        end
    end
    if mod(cnt,20)~=0, fprintf(fid, '\n'); end
    fclose(fid);
end

%% Example 1: Basic configuration (128^3 grid, no padding)
% inputDir    = '/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s';
% filePattern = 'profileIvsTime*.txt';   % e.g., profileIvsTime0.txt, profileIvsTime1.txt, ...
% outputDir   = '/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s/3D';
% 
% opts = struct(); % defaults
% From1Dto3DPlottingInParaviewPrallalel(inputDir, filePattern, outputDir, opts);

%% Example 2: Finer grid, fixed global domain (good for time series)
% inputDir    = '/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s';
% filePattern = 'profileIvsTime*.txt';   % e.g., profileIvsTime0.txt, profileIvsTime1.txt, ...
% outputDir   = '/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s/3D';
% 
% opts = struct();
% opts.Nx = 160; opts.Ny = 160; opts.Nz = 200;
% opts.useParallel = true;
% opts.numWorkers = min( feature('numcores'), 14 ); 
% opts.showProgress = true;
% opts.zRange = [0, 500e-9];     % example fixed z extent
% opts.rMax   = 200e-9;          % example fixed max radius
% opts.rPad   = 0;               % no extra radial padding
% From1Dto3DPlottingInParaviewPrallalel(inputDir, filePattern, outputDir, opts);