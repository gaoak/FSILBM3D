run getParameters.m
%% File serises
readPath  = [casePath '\DatBodySpan\'];
writePath = [casePath '\DatBodySpan\'];
for n = 1:nfile
    time = sTime + (n - 1) * dTime;
    for i = 1:nSolid
        % Generate reading and writing file paths
        [readName, writeName] = generateSolidPath(readPath,writePath,time,i);
        % Read mesh data
        solid = readAscallSolid(readName,time,UVW,Lref,Tref);
        % Write mesh data in vtks
        writeSolidVTK(solid,writeName)
    end
end
