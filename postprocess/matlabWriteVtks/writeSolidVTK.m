function [] = writeSolidVTK(solid,writeFile)
% file path
fid = fopen(writeFile, 'w');

% Write VTK file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK file generated by MATLAB\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');

% Write grid dimensions
fprintf(fid, 'DIMENSIONS %d %d %d\n', solid.nx, solid.ny, solid.nz);

% Write point coordinates
fprintf(fid, 'POINTS %d float\n', solid.nx * solid.ny * solid.nz);

% Write coordinates 
for i = 1:solid.nx
    fprintf(fid, '%f %f %f\n', solid.x(i), solid.y(i), solid.zl(i));
end
for k = 1:solid.nx
    fprintf(fid, '%f %f %f\n', solid.x(k), solid.y(k), solid.zr(k));
end
% Write point data
fprintf(fid, 'POINT_DATA %d\n', solid.nx * solid.ny * solid.nz);

% Close the file
fprintf('Writing ready : %s\n', writeFile)
fclose all;
end
