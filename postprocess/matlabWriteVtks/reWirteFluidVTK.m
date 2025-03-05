run getParameters.m
%% File serises
readPath  = [casePath '\DatFlow\'];
writePath = [casePath '\DatFlow\'];
for n = 1:nfile
    time = sTime + (n - 1) * dTime;
    for i = nBlock:-1:1
        sonBlocks = length(meshContain{i});
        % Read father mesh data
        [readNameFather, writeNameFather] = generateFluidPath(readPath,writePath,time,i);
        meshFather = readBinaryFluid(readNameFather,time,UVW,Uref,Lref,Tref);  
        % Read son mesh data
        for j = 1:sonBlocks
            [readNameSon, writeNameSon] = generateFluidPath(readPath,writePath,time,meshContain{i}(j));
            meshSon = readBinaryFluid(readNameSon,time,UVW,Uref,Lref,Tref); 
            % Update mesh data in coarse mesh
            meshFather = finerToCoarse(meshSon,meshFather);
            fprintf('deliver data from block %d to block %d\n',meshContain{i}(j),i)
        end
        % Write mesh data in vtks
        if ~isOnlyWriteRootBlock 
            writeFluidVTK(meshFather,writeNameFather)
        end
    end
    if isOnlyWriteRootBlock 
        writeFluidVTK(meshFather,writeNameFather)
    end
end
% write whole block contains the calculation domain at all times
writeWholeVTK(meshFather,time,UVW,Lref,Tref,writePath)
