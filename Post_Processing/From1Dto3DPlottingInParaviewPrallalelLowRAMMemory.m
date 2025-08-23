function From1Dto3DPlottingInParaviewPrallalelLowRAMMemory(inputDir, filePattern, outputDir, opts)
% GENERATE_FRUSTUM_VOLUME_SERIES (memory-capped + parallel + progress + streaming)
% - Revolve 2D profiles R(z) to 3D binary volumes (inside=1, outside=0).
% - Stream-write VTK per z-slice (no full 3D array in RAM).
% - Auto-cap parallel workers under a RAM budget (default 20 GB).
% - Progress indicator compatible with parfor.
%
% opts:
%   .Nx,.Ny,.Nz        grid size (default 128,128,128)
%   .zPad,.rPad        padding (default 0)
%   .zRange            fixed [zmin zmax] (optional)
%   .rMax              fixed max radius (optional)
%   .insideVal         uint8 inside value (default 1)
%   .outsideVal        uint8 outside value (default 0)
%   .ascii             true (default)
%   .useParallel       default true if PCT available
%   .numWorkers        desired workers (default = min(feature('numcores'),16))
%   .showProgress      default true
%   .ramBudgetGB       total RAM budget for this job (default 20)
%   .memSafetyFactor   fraction of budget usable (default 0.8)
%   .perWorkerOverheadMB  extra per-worker margin (default 64)
%
% NOTE: Writes VTK legacy STRUCTURED_POINTS in ASCII. Data streamed per slice.

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

    havePCT     = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    defaultNW   = min(feature('numcores'),16);
    useParallel = getOpt(opts,'useParallel', havePCT);
    reqWorkers  = getOpt(opts,'numWorkers', defaultNW);

    ramBudgetGB      = getOpt(opts,'ramBudgetGB', 20);     % total budget
    memSafetyFactor  = getOpt(opts,'memSafetyFactor', 0.8);% use only 80% of it
    perWorkerOverMB  = getOpt(opts,'perWorkerOverheadMB', 64);

    if ~isfolder(outputDir), mkdir(outputDir); end

    % ---- List & order files ----
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

    % ---- Grid and rho(x,y) as single (reduce memory) ----
    x = linspace(-rmax, +rmax, Nx);
    y = linspace(-rmax, +rmax, Ny);
    z = linspace(  zmin,   zmax, Nz);
    dx = (x(end)-x(1)) / max(1,Nx-1);
    dy = (y(end)-y(1)) / max(1,Ny-1);
    dz = (z(end)-z(1)) / max(1,Nz-1);

    [X,Y] = ndgrid(x,y);
    rhoXY = single(hypot(X,Y));  % Nx×Ny single precision

    dims   = [Nx Ny Nz];
    origin = [x(1) y(1) z(1)];
    spc    = [dx dy dz];

    % ---- Memory-aware worker cap ----
    % Per worker footprint ~ rhoXY(single): 4*Nx*Ny bytes
    % + one slice mask (logical/uint8): 1*Nx*Ny bytes
    % + overhead margin perWorkerOverMB.
    perWorkerBytes = (4 + 1) * Nx * Ny; % ~5 bytes per voxel on a slice
    perWorkerBytes = perWorkerBytes + perWorkerOverMB*1024*1024;
    totalBudgetBytes = ramBudgetGB * 1024^3 * memSafetyFactor;
    maxWorkersByMem = max(1, floor(double(totalBudgetBytes) / double(perWorkerBytes)));
    numWorkers = min([reqWorkers, maxWorkersByMem, defaultNW]);

    if useParallel
        % start/resize pool accordingly
        p = gcp('nocreate');
        if isempty(p)
            parpool('local', numWorkers);
        elseif p.NumWorkers ~= numWorkers
            delete(p); parpool('local', numWorkers);
        end
        fprintf('Parallel mode: %d workers (capped by memory budget ~ %.1f GB).\n', ...
            gcp().NumWorkers, totalBudgetBytes/1024^3);
    else
        fprintf('Sequential mode enabled.\n');
    end

    % ---- Progress indicator ----
    h = []; nDone = 0;
    cleanupH = onCleanup(@() closeWaitbarIfOpen());
    if showProgress
        try
            h = waitbar(0, sprintf('Processing 0/%d (0%%)', Nfiles), 'Name','RUN');
        catch
            h = []; fprintf('Progress: 0/%d (0%%)\n', Nfiles);
        end
    end
    function updateProgressLocal()
        nDone = nDone + 1;
        if showProgress
            pct = nDone / Nfiles;
            msg = sprintf('Processing %d/%d (%.0f%%)', nDone, Nfiles, 100*pct);
            if ~isempty(h) && isvalid(h)
                try, waitbar(pct, h, msg); catch, fprintf('%s\n', msg); end
            else
                fprintf('%s\n', msg);
            end
        end
    end
    function closeWaitbarIfOpen()
        if ~isempty(h) && isvalid(h), close(h); end
    end

    if useParallel
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) updateProgressLocal());

        parfor k = 1:Nfiles
            processOneFileStream(k, files, idxNums, z, rhoXY, insideVal, outsideVal, ...
                                 writeAscii, dims, origin, spc, outputDir);
            send(q, 1);
        end
    else
        for k = 1:Nfiles
            processOneFileStream(k, files, idxNums, z, rhoXY, insideVal, outsideVal, ...
                                 writeAscii, dims, origin, spc, outputDir);
            updateProgressLocal();
        end
    end

    if ~isempty(h) && isvalid(h)
        waitbar(1, h, sprintf('Processing %d/%d (100%%)', Nfiles, Nfiles));
    end
    fprintf('Done. Wrote %d VTK files to: %s\n', Nfiles, outputDir);
end


% ---------- Per-file task (STREAMING writer, low-memory) ----------
function processOneFileStream(k, files, idxNums, z, rhoXY, insideVal, outsideVal, ...
                              writeAscii, dims, origin, spacing, outputDir)

    Nx=dims(1); Ny=dims(2); Nz=dims(3);

    filePath = fullfile(files(k).folder, files(k).name);
    A = readmatrix(filePath);
    if size(A,2) < 2
        error('File %s must have at least two columns [Z R].', files(k).name);
    end

    Z = A(:,1); R = A(:,2);
    [Z, iu] = unique(Z(:),'stable');
    R = abs(R(:)); R = R(iu);

    % Interpolate R(z) on common z-grid. Outside -> 0
    Rz = interp1(Z, R, z, 'linear', NaN);
    Rz(~isfinite(Rz)) = 0;
    Rz = max(0, Rz);

    % Prepare output file and write header
    outIndex = idxNums(k); if isnan(outIndex), outIndex = k-1; end
    outName  = sprintf('frustum_%04d.vtk', outIndex);
    outPath  = fullfile(outputDir, outName);

    fid = fopen(outPath,'w');
    if fid==-1, error('Cannot open %s for writing.', outPath); end

    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Frustum volume (streaming)\n');
    if writeAscii
        fprintf(fid, 'ASCII\n');
    else
        error('Binary write not implemented.');
    end
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS %d %d %d\n', Nx, Ny, Nz);
    fprintf(fid, 'ORIGIN %.9g %.9g %.9g\n', origin(1), origin(2), origin(3));
    fprintf(fid, 'SPACING %.9g %.9g %.9g\n', spacing(1), spacing(2), spacing(3));
    fprintf(fid, 'POINT_DATA %d\n', Nx*Ny*Nz);
    fprintf(fid, 'SCALARS frustum unsigned_char 1\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');

    % Stream per z-slice: compute mask and write as ASCII rows
    cnt = 0;
    rhoXYs = rhoXY; % local ref
    for kz = 1:Nz
        thr = single(Rz(kz));              % radius at this z (single)
        mask = rhoXYs <= thr;              % Nx×Ny logical
        % write in x-fastest order: loop y (rows) then x
        for j = 1:Ny
            % convert row to uint8 with inside/outside mapping
            row = uint8(mask(:,j)) * insideVal;   % insideVal (1 or chosen)
            % set outsides explicitly to outsideVal if insideVal~=1
            if outsideVal ~= 0
                row(~mask(:,j)) = outsideVal;
                row(mask(:,j))  = insideVal;
            end
            % print line
            fprintf(fid, '%d ', row);
            cnt = cnt + numel(row);
            fprintf(fid, '\n');
        end
    end
    % ensure trailing newline
    if mod(cnt,20)~=0, fprintf(fid, '\n'); end
    fclose(fid);
end


% ---------- Utilities ----------
function val = getOpt(s, f, def)
    if isfield(s,f) && ~isempty(s.(f)), val = s.(f); else, val = def; end
end

function [filesSorted, idxNums] = sortFilesByTrailingNumber(files)
    idxNums = nan(numel(files),1);
    for i = 1:numel(files)
        nm = files(i).name;
        tok = regexp(nm, '(\d+)(?=\.[^.]+$)', 'tokens', 'once');
        if ~isempty(tok), idxNums(i) = str2double(tok{1}); end
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

%% Example 1: Basic configuration (400*256^3 grid, no padding)
% inputDir    = '/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s';
% filePattern = 'profileIvsTime*.txt';   % e.g., profileIvsTime0.txt, profileIvsTime1.txt, ...
% outputDir   = '/home/moreno/Desktop/TODO/Transitorio/results/non-lineal/cal7.5E-6s/3D';
% 
% opts = struct(); % defaults
% opts.useParallel   = true;    % paraleliza por archivo
% opts.numWorkers    = 16;      % solicita 16 (se ajustará si no cabe en RAM)
% opts.ramBudgetGB   = 20;      % presupuesto total de RAM
% opts.memSafetyFactor = 0.8;   % usa solo el 80% del presupuesto
% opts.perWorkerOverheadMB = 64;% margen por worker
% opts.Nx = 256; opts.Ny = 256; opts.Nz = 400; % ejemplo de grid
% 
% From1Dto3DPlottingInParaviewPrallalelLowRAMMemory(inputDir, filePattern, outputDir, opts)
%% Clean pool
% delete(gcp('nocreate'));   % close parallel pool
% clearvars                  % clear all workspace variables
% close all                  % close figures
% java.lang.System.gc();     % ask JVM to free memory
