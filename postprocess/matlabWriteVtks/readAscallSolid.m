function [solid] = readAscallSolid(filePath,time,initialVelocity,Lref,Tref)
if ~isfile(filePath)
    error('Can not found solid mesh files! : %s',filePath);
end
% Read paramters
solid.data = importdata(filePath).data;
solid.length = floor((size(solid.data,1)/2)) + 1;
% Calculate coordinates
solid.x  = solid.data(1:2:solid.length,1) - initialVelocity(1) / Lref * time * Tref;  % dimensionless
solid.yl = solid.data(1:2:solid.length,2) - initialVelocity(2) / Lref * time * Tref;
solid.yr = solid.data(2:2:solid.length,2) - initialVelocity(2) / Lref * time * Tref;
solid.zl = solid.data(1:2:solid.length,3) - initialVelocity(3) / Lref * time * Tref;
solid.zr = solid.data(2:2:solid.length,3) - initialVelocity(3) / Lref * time * Tref;
% Get index numbers
solid.nx = size(solid.x,1);
solid.ny = 1;
solid.nz = 2;
fclose all;
end
