run getParameters.m
% Read begin point data
for i = 1:nSolid
    fileStruct    = importdata([casePath '\DatInfo\FishNodeBegin_00' num2str(i) '.plt']);
    fileData      = fileStruct.data;
    % Update velocity
    fileData(:,5) = fileData(:,5)  - UVW(1) / Uref;
    fileData(:,6) = fileData(:,6)  - UVW(2) / Uref;
    fileData(:,7) = fileData(:,7)  - UVW(3) / Uref;
    % Write data
    fileID = fopen([casePath '\DatInfo\FishNodeBegin_00' num2str(i) '_n.plt'], 'w');
    fprintf(fileID, '%s\n',fileStruct.textdata{1});
    for j = 1:size(fileData,1)
        time = fileData(j,1);
        fileData(j,2) = fileData(j,2)  - UVW(1) / Lref * time * Tref;
        fileData(j,3) = fileData(j,3)  - UVW(2) / Lref * time * Tref;
        fileData(j,4) = fileData(j,4)  - UVW(3) / Lref * time * Tref;
        fprintf(fileID, [repmat('%.8g ', 1, size(fileData,2)), '\n'], fileData(j,:));
    end
    fprintf('Writing file read : %s\n', [casePath '\DatInfo\FishNodeBegin_00' num2str(i) '_n.plt'])
end
% Read center point data
for i = 1:nSolid
    fileStruct    = importdata([casePath '\DatInfo\FishNodeCenter_00' num2str(i) '.plt']);
    fileData      = fileStruct.data;
    % Update velocity
    fileData(:,5) = fileData(:,5)  - UVW(1) / Uref;
    fileData(:,6) = fileData(:,6)  - UVW(2) / Uref;
    fileData(:,7) = fileData(:,7)  - UVW(3) / Uref;
    % Write data
    fileID = fopen([casePath '\DatInfo\FishNodeCenter_00' num2str(i) '_n.plt'], 'w');
    fprintf(fileID, '%s\n',fileStruct.textdata{1});
    for j = 1:size(fileData,1)
        time = fileData(j,1);
        fileData(j,2) = fileData(j,2)  - UVW(1) / Lref * time * Tref;
        fileData(j,3) = fileData(j,3)  - UVW(2) / Lref * time * Tref;
        fileData(j,4) = fileData(j,4)  - UVW(3) / Lref * time * Tref;
        fprintf(fileID, [repmat('%.8g ', 1, size(fileData,2)), '\n'], fileData(j,:));
    end
    fprintf('Writing file read : %s\n', [casePath '\DatInfo\FishNodeCenter_00' num2str(i) '_n.plt'])
end
% Read end point data
for i = 1:nSolid
    fileStruct    = importdata([casePath '\DatInfo\FishNodeEnd_00' num2str(i) '.plt']);
    fileData      = fileStruct.data;
    % Update velocity
    fileData(:,5) = fileData(:,5)  - UVW(1) / Uref;
    fileData(:,6) = fileData(:,6)  - UVW(2) / Uref;
    fileData(:,7) = fileData(:,7)  - UVW(3) / Uref;
    % Write data
    fileID = fopen([casePath '\DatInfo\FishNodeEnd_00' num2str(i) '_n.plt'], 'w');
    fprintf(fileID, '%s\n',fileStruct.textdata{1});
    for j = 1:size(fileData,1)
        time = fileData(j,1);
        fileData(j,2) = fileData(j,2)  - UVW(1) / Lref * time * Tref;
        fileData(j,3) = fileData(j,3)  - UVW(2) / Lref * time * Tref;
        fileData(j,4) = fileData(j,4)  - UVW(3) / Lref * time * Tref;
        fprintf(fileID, [repmat('%.8g ', 1, size(fileData,2)), '\n'], fileData(j,:));
    end
    fprintf('Writing file read : %s\n', [casePath '\DatInfo\FishNodeEnd_00' num2str(i) '_n.plt'])
end
% Read averaging data
for i = 1:nSolid
    fileStruct    = importdata([casePath '\DatInfo\FishNodeMean_00' num2str(i) '.plt']);
    fileData      = fileStruct.data;
    % Update velocity
    fileData(:,5) = fileData(:,5)  - UVW(1) / Uref;
    fileData(:,6) = fileData(:,6)  - UVW(2) / Uref;
    fileData(:,7) = fileData(:,7)  - UVW(3) / Uref;
    % Write data
    fileID = fopen([casePath '\DatInfo\FishNodeMean_00' num2str(i) '_n.plt'], 'w');
    fprintf(fileID, '%s\n',fileStruct.textdata{1});
    for j = 1:size(fileData,1)
        time = fileData(j,1);
        fileData(j,2) = fileData(j,2)  - UVW(1) / Lref * time * Tref;
        fileData(j,3) = fileData(j,3)  - UVW(2) / Lref * time * Tref;
        fileData(j,4) = fileData(j,4)  - UVW(3) / Lref * time * Tref;
        fprintf(fileID, [repmat('%.8g ', 1, size(fileData,2)), '\n'], fileData(j,:));
    end
    fprintf('Writing file read : %s\n', [casePath '\DatInfo\FishNodeMean_00' num2str(i) '_n.plt'])
end
