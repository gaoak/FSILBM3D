clear;clc;close all;format long
%% Case path
casePath  = 'G:\TandemPlates\Examples\TwoMeshBlocks-1';
%% Read key lines
UrefLine  = readKeyLines([casePath '\check.dat'] ,'Uref',1);
UVWLine   = readKeyLines([casePath '\inFlow.dat'],'uvwIn',1);
TimeLine  = readKeyLines([casePath '\inFlow.dat'],'timeWriteBegin',1);
deltaLine = readKeyLines([casePath '\inFlow.dat'],'timeWriteFlow',1); 
if (str2double(deltaLine{1}) ~= str2double(deltaLine{2}))
    error('solid writing interval do not equals fluid writing interval!')
end
fluidLine = readKeyLines([casePath '\inFlow.dat'],'nblock',1);
solidLine = readKeyLines([casePath '\inFlow.dat'],'nFish',1); % writing intervl
%% Read key parameters
Uref   = str2double(UrefLine{3}); % reference velocity
UVW    = [str2double(UVWLine{1}) str2double(UVWLine{2}) str2double(UVWLine{3})]; % inflow velocity
sTime  = str2double(TimeLine{1}); % begine time for writing
eTime  = str2double(TimeLine{2}); % ending time for writing
dTime  = str2double(deltaLine{1}); % writing intervl
nBlock = floor(str2double(fluidLine{1})); % the number of fluid blocks
nSolid = floor(str2double(solidLine{1})); % the number of solids
meshContain = cell(nBlock); % the contian relationship of the fluid blocks, the space means no son block
for i = 1:nBlock
    blockLine = readKeyLines([casePath '\check.dat'],'sonBlocks',i);
    if length(blockLine) == 2
        meshContain{i} = [ ];
    else
        for j = 1:length(blockLine)-2
        meshContain{i} = [meshContain{i} str2double(blockLine{2+j})];
        end
    end
end
%% Calculate key parameters
nfile  = (eTime - sTime) / dTime + 1;
exVel  = UVW / Uref;