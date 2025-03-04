function [solid] = readAscallSolid(filePath,extraVelocity,time)
% Read paramters
solid.data = importdata(filePath).data;
solid.length = floor((size(solid.data,1)/2)) + 1;
% Calculate coordinates
solid.x = solid.data(1:2:solid.length,1);
solid.y = solid.data(1:2:solid.length,2);
solid.zl = solid.data(1:2:solid.length,3);
solid.zr = solid.data(2:2:solid.length,3);
solid.x = solid.x - extraVelocity(1) * time * 1e-5;
solid.y = solid.y - extraVelocity(2) * time * 1e-5;
solid.zl = solid.zl - extraVelocity(3) * time * 1e-5;
solid.zr = solid.zr - extraVelocity(3) * time * 1e-5;
% Get index numbers
solid.nx = size(solid.x,1);
solid.ny = 1;
solid.nz = 2;
end
