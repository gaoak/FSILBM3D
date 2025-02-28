run globPalarameters.m
%% File serises
readPath  = [casePath '\DatFlow\'];
writePath = [casePath '\DatFlow\'];
for n = 1:nfile
    % Generate reading and writing file paths
    time = sTime + (n - 1) * dTime;
    [readName2, writeName2] = generateFluidPath(readPath,writePath,time,2);
    [readName1, writeName1] = generateFluidPath(readPath,writePath,time,1);
    % Read mesh data
    mesh2 = readBinaryFluid(readName2,exVel,time);  % finer mesh
    mesh1 = readBinaryFluid(readName1,exVel,time);  % coarse mesh
    % Update mesh data in coarse mesh
    mesh1 = finerToCoarse(mesh2,mesh1);
    % Write mesh data in vtks
    writeFluidVTK(mesh2,writeName2)
    writeFluidVTK(mesh1,writeName1)
end
