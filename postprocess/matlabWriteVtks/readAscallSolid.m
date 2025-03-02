function [solid] = readAscallSolid(filePath,extraVelocity,time)
% Read paramters
solid.data = importdata(filePath).data;
length = floor((size(solid.data,1)/2));
% Calculate coordinates
solid.x = solid.data(1:2:length,1);
solid.y = solid.data(1:2:length,2);
solid.z = [solid.data(1,3); solid.data(2,3)];
solid.x = solid.x - extraVelocity(1) * time * 1e-5;
solid.y = solid.y - extraVelocity(2) * time * 1e-5;
solid.z = solid.z - extraVelocity(3) * time * 1e-5;
% Get index numbers
solid.nx = size(solid.x,1);
solid.ny = 1;
solid.nz = size(solid.z,1);
end
