function [solid] = readAscallSolid(filePath,extraVelocity,time)
if ~isfile(filePath)
    error('Can not found solid mesh files! : %s',filePath);
end
% Read paramters
solid.data = importdata(filePath).data;
solid.length = floor((size(solid.data,1)/2)) + 1;
% Calculate coordinates
solid.x  = solid.data(1:2:solid.length,1) - extraVelocity(1) * time;
solid.yl = solid.data(1:2:solid.length,2) - extraVelocity(1) * time;
solid.yr = solid.data(2:2:solid.length,2) - extraVelocity(1) * time;
solid.zl = solid.data(1:2:solid.length,3) - extraVelocity(1) * time;
solid.zr = solid.data(2:2:solid.length,3) - extraVelocity(1) * time;
% Get index numbers
solid.nx = size(solid.x,1);
solid.ny = 1;
solid.nz = 2;
fclose(filePath);
end
