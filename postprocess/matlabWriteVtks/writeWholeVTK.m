function [] = writeWholeVTK(mesh,inflowVelocity,time,Lref,Tref,writePath)
if norm(inflowVelocity) ~= 0
    % Update whole block parameters
    exMesh = mesh;
    exDisplacemnt = time * Tref * inflowVelocity / Lref;
    exMesh.u = zeros(size(exMesh.u));
    exMesh.v = zeros(size(exMesh.v));
    exMesh.w = zeros(size(exMesh.w));
    % Calculate enlarge coefficients
    enlarge.x = (exMesh.xmax - exMesh.xmin + exDisplacemnt(1)) / (exMesh.xmax - exMesh.xmin);
    enlarge.y = (exMesh.ymax - exMesh.ymin + exDisplacemnt(2)) / (exMesh.ymax - exMesh.ymin);
    enlarge.z = (exMesh.zmax - exMesh.zmin + exDisplacemnt(3)) / (exMesh.zmax - exMesh.zmin);
    % Enlarge whole block
    exMesh.x = (exMesh.x - exMesh.xmax) * enlarge.x + exMesh.xmax + exDisplacemnt(1);
    exMesh.y = (exMesh.y - exMesh.ymax) * enlarge.y + exMesh.ymax + exDisplacemnt(2);
    exMesh.z = (exMesh.z - exMesh.zmax) * enlarge.z + exMesh.zmax + exDisplacemnt(3);
    % Write whole block
    fprintf("Writing whole block for moving mesh\n")
    writeFluidVTK(exMesh,[writePath 'b000Flow0000000000.vtk'])
end
end
