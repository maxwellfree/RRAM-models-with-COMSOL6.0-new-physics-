function generate3Dfrom1D(inputFolder, outputFolder, nTheta, gridResolution)
% generate3Dfrom1D - Convert 1D curve data into 3D binary volumes for Paraview
%
% INPUTS:
%   inputFolder    - folder containing the 1D files (two columns: [r, z])
%   outputFolder   - folder where 3D volumes will be saved
%   nTheta         - number of angular divisions around Z (e.g. 180)
%   gridResolution - number of points in radial and vertical directions (e.g. 100)
%
% Each file in inputFolder is assumed to represent curve data at a given time.
% The output will be a set of .vti files (VTK ImageData) named sequentially.

    % Create output folder if it does not exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % List all input files
    files = dir(fullfile(inputFolder, 'pro*.txt')); % 
    
    for f = 1:length(files)
        % ---- Load curve data ----
        data = load(fullfile(inputFolder, files(f).name));
        zCurve = data(:,1); % radial coordinate
        rCurve = data(:,2); % vertical coordinate
        
        % ---- Define 3D grid ----
        rMax = max(rCurve);
        zMin = min(zCurve);
        zMax = max(zCurve);
        
        r = linspace(0, rMax, gridResolution);
        z = linspace(zMin, zMax, gridResolution);
        theta = linspace(0, 2*pi, nTheta);
        
        [R,Theta,Z] = ndgrid(r,theta,z);
        X = R .* cos(Theta);
        Y = R .* sin(Theta);
        
        % ---- Interpolate curve to full Z-range ----
        rInterp = interp1(zCurve, rCurve, z, 'linear','extrap');
        
        % ---- Create 3D binary volume ----
        volume = zeros(size(R));
        for k = 1:length(z)
            % For each z-slice, mark all radii below curve as 1
            mask = R(:,:,k) <= rInterp(k);
            volume(:,:,k) = mask;
        end
        
        % ---- Write to VTK ImageData (.vti) ----
        outputFile = fullfile(outputFolder, sprintf('frame_%04d.vti', f));
        writeVTKvolume(volume, X, Y, Z, outputFile);
        
        fprintf('Saved %s\n', outputFile);
    end
end


function writeVTKvolume(V, X, Y, Z, filename)
% writeVTKvolume - Save a 3D matrix as a VTK ImageData (.vti)
% This function uses simple ASCII XML format.

    dims = size(V);
    spacing = [X(2,1,1)-X(1,1,1), Y(1,2,1)-Y(1,1,1), Z(1,1,2)-Z(1,1,1)];
    origin = [min(X(:)), min(Y(:)), min(Z(:))];
    
    fid = fopen(filename, 'w');
    fprintf(fid, '<?xml version="1.0"?>\n');
    fprintf(fid, '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid, '<ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="%f %f %f" Spacing="%f %f %f">\n', ...
        dims(1)-1, dims(2)-1, dims(3)-1, origin(1), origin(2), origin(3), spacing(1), spacing(2), spacing(3));
    fprintf(fid, '<Piece Extent="0 %d 0 %d 0 %d">\n', dims(1)-1, dims(2)-1, dims(3)-1);
    fprintf(fid, '<PointData Scalars="scalars">\n');
    fprintf(fid, '<DataArray type="Int32" Name="scalars" format="ascii">\n');
    fprintf(fid, '%d ', V(:));
    fprintf(fid, '\n</DataArray>\n');
    fprintf(fid, '</PointData>\n');
    fprintf(fid, '<CellData>\n</CellData>\n');
    fprintf(fid, '</Piece>\n');
    fprintf(fid, '</ImageData>\n');
    fprintf(fid, '</VTKFile>\n');
    fclose(fid);
end

% Sample 
% generate3Dfrom1D('/home/moreno/Desktop/TODO/Transitorio/cal', 'cal/3D', 180, 100)