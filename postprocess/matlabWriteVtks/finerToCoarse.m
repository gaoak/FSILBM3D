function [meshFather] = finerToCoarse(meshson,meshFather)
% get indexs
nz1 = round((meshson.zmin - meshFather.zmin) / meshFather.dh) + 1;
nz2 = nz1 + (meshson.nz - 1)/2;
ny1 = round((meshson.ymin - meshFather.ymin) / meshFather.dh) + 1;
ny2 = ny1 + (meshson.ny - 1)/2;
nx1 = round((meshson.xmin - meshFather.xmin) / meshFather.dh) + 1;
nx2 = nx1 + (meshson.nx - 1)/2;
% update date in coarse mesh
meshFather.u(nx1:nx2,ny1:ny2,nz1:nz2) = meshson.u(1:2:meshson.nx,1:2:meshson.ny,1:2:meshson.nz);
meshFather.v(nx1:nx2,ny1:ny2,nz1:nz2) = meshson.v(1:2:meshson.nx,1:2:meshson.ny,1:2:meshson.nz);
meshFather.w(nx1:nx2,ny1:ny2,nz1:nz2) = meshson.w(1:2:meshson.nx,1:2:meshson.ny,1:2:meshson.nz);
end
