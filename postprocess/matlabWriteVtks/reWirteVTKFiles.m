clear;clc;close all
%% File serises
nfile = 81;
sTime = 00000;
dTime = 10000;
readPath  = 'G:\TandemPlates\Demos\FlexiblePlate\DatFlow\';
writePath = 'G:\TandemPlates\Demos\FlexiblePlate\DatFlow\';
for n = 1:nfile
    % Generate reading and writing file paths
    time = sTime + (n - 1) * dTime;
    [readName2, writeName2] = generateFilePath(readPath,writePath,time,2);
    [readName1, writeName1] = generateFilePath(readPath,writePath,time,1);
    % Read data
    mesh2 = readBinaryFile(readName2);  % finer mesh
    mesh1 = readBinaryFile(readName1);  % coarse mesh
    % Update data in coarse mesh
    mesh1 = finerToCoarse(mesh2,mesh1);
    % Write tecplot file
    writeParaviewVTK(mesh2,writeName2)
    writeParaviewVTK(mesh1,writeName1)
end