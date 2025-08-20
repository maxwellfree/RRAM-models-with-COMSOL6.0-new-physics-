% generar_volumenes_superficie.m
% Toma superficies 2D (r, z, valor) o (z, r, valor) y las revoluciona
% alrededor del eje Z para producir un volumen 3D escalar listo para ParaView.
%
% Campo 3D resultante:
%   V(x,y,z) = f( r = sqrt(x^2 + y^2), z ), independiente del ángulo theta.
%
% Salida: archivos VTK secuenciales (filamento_0000.vtk, filamento_0001.vtk, ...)
% y opcionalmente un .pvd con la colección temporal.

%% ======================= CONFIGURACIÓN =======================
inputFolder    = "entrada_superficies";    % Carpeta con las superficies 2D
filePattern    = "*.txt";                  % Patrón de archivos
hasHeader      = false;                    % true si la primera fila es cabecera
delimiter      = "";                       % "" para autodetección

% Orden/semántica de columnas en los archivos de entrada:
%   Elige una de: 'r_z_val'  (col1=r, col2=z, col3=valor)
%                 'z_r_val'  (col1=z, col2=r, col3=valor)
columnOrder    = "r_z_val";

% Nombre del escalar para VTK (ej.: 'voltaje', 'temperatura')
scalarName     = "voltaje";

% Carpeta y nombres de salida
outputFolder   = "salida_volumen_3D";
outputPrefix   = "filamento_";
padNumber      = 4;                        % ceros a la izquierda (0000, 0001, ...)

% Resolución del volumen (x-y-z). nx,ny fijan el dominio radial; nz en altura.
nx = 128;
ny = 128;
nz = 256;

% Margen radial sobre r_max de los datos (para encuadre)
radialPadding  = 0.05;    % 5%

% Interpolación del campo f(r,z):
%   'auto'    -> usa griddedInterpolant si la malla es regular, si no scatteredInterpolant
%   'linear'  -> fuerza scatteredInterpolant('linear','none')
%   'natural' -> scatteredInterpolant('natural','none')
interpMode     = "auto";

% Valor fuera del dominio de datos (r,z) al extrapolar:
%   Usa NaN para poder "enmascarar" fuera; o un valor constante (p.ej. 0)
outsideValue   = NaN;

% Generar archivo .pvd con la serie temporal
emitirPVD      = true;
pvdName        = "serie_superficie.pvd";
%% =============================================================

if ~isfolder(inputFolder), error("No existe carpeta de entrada: %s", inputFolder); end
if ~isfolder(outputFolder), mkdir(outputFolder); end

files = dir(fullfile(inputFolder, filePattern));
if isempty(files), error("No se encontraron archivos con patrón %s", filePattern); end

% Ordenar de forma robusta por sufijo numérico (si existe) y luego alfabético
names = string({files.name});
idxNums = extractNumericSuffix(names);
[~, order] = sortrows([idxNums, names.'], [1 2]);
files = files(order);

% --------- PASO 1: Explorar dominios globales (r,z) ----------
RminG = +inf; RmaxG = -inf; ZminG = +inf; ZmaxG = -inf;
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    M = readmatrixSmart(fpath, delimiter, hasHeader);
    [r, z, val] = parseColumns(M, columnOrder);
    [r, z, val] = cleanTriplets(r, z, val);
    RminG = min(RminG, min(r));   RmaxG = max(RmaxG, max(r));
    ZminG = min(ZminG, min(z));   ZmaxG = max(ZmaxG, max(z));
end

% Dominios de la grilla 3D
Rmax = RmaxG * (1 + radialPadding);
x = linspace(-Rmax, Rmax, nx);
y = linspace(-Rmax, Rmax, ny);
zGrid = linspace(ZminG, ZmaxG, nz);

% Malla XY y radio
[X, Y] = meshgrid(x, y);          % (ny x nx)
Rxy = sqrt(X.^2 + Y.^2);

% Parámetros VTK (ImageData / STRUCTURED_POINTS)
dx = (nx>1) * (x(2)-x(1)) + (nx==1)*1.0;
dy = (ny>1) * (y(2)-y(1)) + (ny==1)*1.0;
dz = (nz>1) * (zGrid(2)-zGrid(1)) + (nz==1)*1.0;
origin  = [x(1), y(1), zGrid(1)];
spacing = [dx, dy, dz];

vtkFiles = strings(numel(files),1);
fprintf("Generando %d pasos temporales...\n", numel(files));

% --------- PASO 2: Procesar cada superficie -> volumen 3D ----
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    M = readmatrixSmart(fpath, delimiter, hasHeader);
    [r, z, val] = parseColumns(M, columnOrder);
    [r, z, val] = cleanTriplets(r, z, val);

    % Construir interpolante f(r,z) -> valor
    F = buildInterpolator(r, z, val, interpMode, outsideValue);

    % Volumen 3D: V(x,y,z) = f( sqrt(x^2+y^2), z )
    V = nan(ny, nx, nz, 'double');
    for kk = 1:nz
        rv = Rxy(:);
        zv = zGrid(kk) * ones(numel(rv),1);
        Vk = F(rv, zv);
        % Forzar fuera del cilindro simulado (r>Rmax) al valor externo
        if ~isnan(outsideValue)
            Vk(rv > Rmax) = outsideValue;
        else
            Vk(rv > Rmax) = NaN;
        end
        V(:,:,kk) = reshape(Vk, ny, nx);
    end

    % Reordenar a (x,y,z) y vectorizar para VTK (x rápido, luego y, luego z)
    dataVTK = permute(V, [2 1 3]);            % (nx, ny, nz)
    dataVTK = dataVTK(:);

    % Escribir archivo VTK (ASCII). Si hay NaN, ParaView los soporta.
    fname   = sprintf("%s%0*d.vtk", outputPrefix, padNumber, k-1);
    outPath = fullfile(outputFolder, fname);
    writeVTKStructuredPointsScalar(outPath, nx, ny, nz, origin, spacing, dataVTK, scalarName);
    vtkFiles(k) = fname;

    fprintf("  [%d/%d] -> %s\n", k, numel(files), fname);
end

% --------- PASO 3: Escribir .pvd opcional ---------------------
if emitirPVD
    pvdPath = fullfile(outputFolder, pvdName);
    writePVD(pvdPath, vtkFiles);
    fprintf("Serie PVD escrita: %s\n", pvdPath);
end

fprintf("¡Listo! Serie 3D en: %s\n", outputFolder);

%% ===================== FUNCIONES LOCALES =====================

function M = readmatrixSmart(fpath, delimiter, hasHeader)
    opts = detectImportOptions(fpath, 'NumHeaderLines', 0);
    if ~isempty(delimiter)
        opts = setvaropts(opts, 'WhitespaceRule', 'preserve');
        opts.Delimiter = delimiter;
    end
    if hasHeader
        opts.DataLines = [2 inf];
    end
    M = readmatrix(fpath, opts);
end

function [r, z, val] = parseColumns(M, columnOrder)
    if size(M,2) < 3
        error("Cada archivo debe tener al menos 3 columnas: r, z y valor (o z, r y valor).");
    end
    switch string(columnOrder)
        case "r_z_val"
            r   = M(:,1); z = M(:,2); val = M(:,3);
        case "z_r_val"
            z   = M(:,1); r = M(:,2); val = M(:,3);
        otherwise
            error("columnOrder debe ser 'r_z_val' o 'z_r_val'.");
    end
end

function [r, z, val] = cleanTriplets(r, z, val)
    % Elimina NaN y ordena por (z, r) para mayor estabilidad
    mask = ~(isnan(r) | isnan(z) | isnan(val));
    r = r(mask); z = z(mask); val = val(mask);
    A = [z(:), r(:), val(:)];
    A = sortrows(A, [1 2]);  % por z luego r
    z = A(:,1); r = A(:,2); val = A(:,3);

    % Colapsa duplicados exactos (r,z) haciendo media del valor
    [ZR, ~, ic] = unique([z r], "rows");
    if numel(ic) < size(A,1)
        val = accumarray(ic, val, [], @mean);
        z = ZR(:,1); r = ZR(:,2);
    end
end

function F = buildInterpolator(r, z, val, mode, outsideValue)
    % Determina si los datos están en malla regular (tensor) y crea interpolante.
    ru = unique(r);
    zu = unique(z);
    isRegular = (numel(ru)*numel(zu) == numel(r));
    if isRegular
        % Intento de reordenar a grilla (r,z)
        [Rgrid, Zgrid] = ndgrid(ru, zu);
        % Mapear (r,z)->índices
        [~, ir] = ismember(r, ru);
        [~, iz] = ismember(z, zu);
        Vgrid = nan(numel(ru), numel(zu));
        Vgrid(sub2ind(size(Vgrid), ir, iz)) = val;
        % Completar posibles huecos con interpolación scattered puntual
        if any(isnan(Vgrid(:)))
            S = scatteredInterpolant(r, z, val, 'linear', 'none');
            Vfill = S(Rgrid, Zgrid);
            Vgrid(isnan(Vgrid)) = Vfill(isnan(Vgrid));
        end
        % griddedInterpolant espera orden creciente por dimensión
        Fg = griddedInterpolant({ru, zu}, Vgrid, 'linear', 'none');
        F = @(rq, zq) fillOutside(Fg(rq, zq), rq, zq, ru([1 end]), zu([1 end]), outsideValue);
    else
        % Disperso: scatteredInterpolant
        meth = 'linear';
        if mode == "natural", meth = 'natural'; end
        S = scatteredInterpolant(r, z, val, meth, 'none');
        F = @(rq, zq) fillOutside(S(rq, zq), rq, zq, [min(r) max(r)], [min(z) max(z)], outsideValue);
    end
end

function v = fillOutside(v, rq, zq, rlim, zlim, outsideValue)
    maskOut = (rq < rlim(1)) | (rq > rlim(2)) | (zq < zlim(1)) | (zq > zlim(2));
    if any(maskOut, 'all')
        if isnan(outsideValue)
            v(maskOut) = NaN;
        else
            v(maskOut) = outsideValue;
        end
    end
end

function nums = extractNumericSuffix(names)
    nums = inf(numel(names),1);
    for i=1:numel(names)
        t = regexp(names(i), '(\d+)(?!.*\d)', 'tokens', 'once');
        if ~isempty(t), nums(i) = str2double(t{1}); end
    end
end

function writeVTKStructuredPointsScalar(filename, nx, ny, nz, origin, spacing, dataVec, scalarName)
% Escribe un VTK Legacy ASCII con DATASET STRUCTURED_POINTS y un SCALAR.
    fid = fopen(filename, 'w');
    if fid==-1, error("No se puede escribir %s", filename); end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, "# vtk DataFile Version 3.0\n");
    fprintf(fid, "Volumen 3D por revolucion con campo escalar\n");
    fprintf(fid, "ASCII\n");
    fprintf(fid, "DATASET STRUCTURED_POINTS\n");
    fprintf(fid, "DIMENSIONS %d %d %d\n", nx, ny, nz);
    fprintf(fid, "ORIGIN %.10g %.10g %.10g\n", origin(1), origin(2), origin(3));
    fprintf(fid, "SPACING %.10g %.10g %.10g\n", spacing(1), spacing(2), spacing(3));
    fprintf(fid, "POINT_DATA %d\n", nx*ny*nz);
    fprintf(fid, "SCALARS %s double 1\n", scalarName);
    fprintf(fid, "LOOKUP_TABLE default\n");

    % dataVec es double y puede contener NaN (ParaView los interpreta)
    lineLen = 6; % valores por línea (para legibilidad)
    for i = 1:lineLen:numel(dataVec)
        j = min(i+lineLen-1, numel(dataVec));
        fprintf(fid, '%.15g ', dataVec(i:j));
        fprintf(fid, '\n');
    end
end

function writePVD(pvdPath, vtkFiles)
    fid = fopen(pvdPath, 'w');
    if fid==-1, error("No se puede escribir %s", pvdPath); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '  <Collection>\n');
    for k=1:numel(vtkFiles)
        fprintf(fid, '    <DataSet timestep="%g" group="" part="0" file="%s"/>\n', k-1, vtkFiles(k));
    end
    fprintf(fid, '  </Collection>\n');
    fprintf(fid, '</VTKFile>\n');
end