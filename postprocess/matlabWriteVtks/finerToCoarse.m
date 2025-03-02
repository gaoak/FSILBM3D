function [meshCoarse] = finerToCoarse(meshFiner,meshCoarse)
% get indexs
nz1 = (meshFiner.zmin - meshCoarse.zmin) / meshCoarse.dh + 1;
nz2 = nz1 + (meshFiner.nz - 1)/2;
ny1 = (meshFiner.ymin - meshCoarse.ymin) / meshCoarse.dh + 1;
ny2 = ny1 + (meshFiner.ny - 1)/2;
nx1 = (meshFiner.xmin - meshCoarse.xmin) / meshCoarse.dh + 1;
nx2 = nx1 + (meshFiner.nx - 1)/2;
% update date in coarse mesh
meshCoarse.u(nx1:nx2,ny1:ny2,nz1:nz2) = meshFiner.u(1:2:meshFiner.nx,1:2:meshFiner.ny,1:2:meshFiner.nz);
meshCoarse.v(nx1:nx2,ny1:ny2,nz1:nz2) = meshFiner.v(1:2:meshFiner.nx,1:2:meshFiner.ny,1:2:meshFiner.nz);
meshCoarse.w(nx1:nx2,ny1:ny2,nz1:nz2) = meshFiner.w(1:2:meshFiner.nx,1:2:meshFiner.ny,1:2:meshFiner.nz);
end
