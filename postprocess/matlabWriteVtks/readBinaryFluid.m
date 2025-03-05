function [mesh] = readBinaryFluid(filePath,time,inflowVelocity,Uref,Lref,Tref)
% Open file
fileID = fopen(filePath,'r');
if fileID == -1
    error('Can not found fluid mesh files! : %s',filePath);
end
% Read paramters
mesh.nx = fread(fileID, 1, 'int32');
mesh.ny = fread(fileID, 1, 'int32');
mesh.nz = fread(fileID, 1, 'int32');
mesh.id = fread(fileID, 1, 'int32');
mesh.xmin = fread(fileID, 1, 'double') - inflowVelocity(1) / Lref * time * Tref;  % dimensionless
mesh.ymin = fread(fileID, 1, 'double') - inflowVelocity(2) / Lref * time * Tref;
mesh.zmin = fread(fileID, 1, 'double') - inflowVelocity(3) / Lref * time * Tref;
mesh.dh = fread(fileID, 1, 'double');
% Error
if mod(mesh.nx,2) == 0 || mod(mesh.ny,2) == 0 || mod(mesh.nz,2) == 0
    error('The number of grids in %s is odd!', filePath)
end
% Calculate max coordinates
mesh.xmax = mesh.xmin + (mesh.nx - 1) * mesh.dh;
mesh.ymax = mesh.ymin + (mesh.ny - 1) * mesh.dh;
mesh.zmax = mesh.zmin + (mesh.nz - 1) * mesh.dh;
% Read velocities
u = fread(fileID, mesh.nz * mesh.ny * mesh.nx, 'float32') - inflowVelocity(1) / Uref; % k,j,i
v = fread(fileID, mesh.nz * mesh.ny * mesh.nx, 'float32') - inflowVelocity(2) / Uref; % dimensionless
w = fread(fileID, mesh.nz * mesh.ny * mesh.nx, 'float32') - inflowVelocity(3) / Uref;
mesh.u = permute(reshape(u, [mesh.nz mesh.ny mesh.nx]),[3 2 1]); % i,j,k
mesh.v = permute(reshape(v, [mesh.nz mesh.ny mesh.nx]),[3 2 1]);
mesh.w = permute(reshape(w, [mesh.nz mesh.ny mesh.nx]),[3 2 1]);
% Calculate coordinates
array.x = mesh.xmin : mesh.dh : mesh.xmax;
array.y = mesh.ymin : mesh.dh : mesh.ymax;
array.z = mesh.zmin : mesh.dh : mesh.zmax;
[mesh.x, mesh.y, mesh.z] = ndgrid(array.x, array.y, array.z);
fclose all;
end
