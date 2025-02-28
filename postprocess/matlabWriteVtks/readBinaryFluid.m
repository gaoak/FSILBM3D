function [mesh] = readBinaryFluid(filePath,extraVelocity,time)
% Open file
fileID = fopen(filePath,'r');
% Read paramters
mesh.nx = fread(fileID, 1, 'int32');
mesh.ny = fread(fileID, 1, 'int32');
mesh.nz = fread(fileID, 1, 'int32');
mesh.id = fread(fileID, 1, 'int32');
mesh.xmin = fread(fileID, 1, 'double');
mesh.ymin = fread(fileID, 1, 'double');
mesh.zmin = fread(fileID, 1, 'double');
mesh.dh = fread(fileID, 1, 'double');
% Read velocities
u = fread(fileID, mesh.nz * mesh.ny * mesh.nx, 'float32'); % k,j,i
v = fread(fileID, mesh.nz * mesh.ny * mesh.nx, 'float32');
w = fread(fileID, mesh.nz * mesh.ny * mesh.nx, 'float32');
mesh.u = permute(reshape(u, [mesh.nz mesh.ny mesh.nx]),[3 2 1]) + extraVelocity(1); % i,j,k
mesh.v = permute(reshape(v, [mesh.nz mesh.ny mesh.nx]),[3 2 1]) + extraVelocity(2);
mesh.w = permute(reshape(w, [mesh.nz mesh.ny mesh.nx]),[3 2 1]) + extraVelocity(3);
% Calculate coordiantes
[mesh.x, mesh.y, mesh.z] = ndgrid(mesh.xmin : mesh.dh : (mesh.xmin + (mesh.nx - 1) * mesh.dh), ...
                                  mesh.ymin : mesh.dh : (mesh.ymin + (mesh.ny - 1) * mesh.dh), ...
                                  mesh.zmin : mesh.dh : (mesh.zmin + (mesh.nz - 1) * mesh.dh));
mesh.x = mesh.x - extraVelocity(1) * time;
mesh.y = mesh.x - extraVelocity(2) * time;
mesh.z = mesh.x - extraVelocity(3) * time;
end