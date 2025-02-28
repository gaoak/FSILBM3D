run globalParameters.m
%% File serises
readPath  = [casePath '\DatBodySpan\'];
writePath = [casePath '\DatBodySpan\'];
for n = 1:nfile
    % Generate reading and writing file paths
    time = sTime + (n - 1) * dTime;
    [readName2, writeName2] = generateSolidPath(readPath,writePath,time,2);
    [readName1, writeName1] = generateSolidPath(readPath,writePath,time,1);
    % Read mesh data
    % solid2 = readAscallSolid(readName2,exVel,time);
    solid1 = readAscallSolid(readName1,exVel,time);
    % Write mesh data in vtks
    % writeSolidVTK(solid2,writeName2)
    writeSolidVTK(solid1,writeName1)
end
