clear;clc;close all;format long
%% Case path
isOnlyWriteRootBlock = true;
casePath  = 'G:\TandemPlates\Validation\Comparison\Case3DL064T06400-double';
%% Read key lines
UrefLine  = readKeyLines([casePath '\check.dat' ],'Uref'          ,1);
TrefLine  = readKeyLines([casePath '\check.dat' ],'Tref'          ,1);
LrefLine  = readKeyLines([casePath '\check.dat' ],'Lref'          ,1);
UVWLine   = readKeyLines([casePath '\inFlow.dat'],'uvwIn'         ,1);
TimeLine  = readKeyLines([casePath '\inFlow.dat'],'timeWriteBegin',1);
deltaLine = readKeyLines([casePath '\inFlow.dat'],'timeWriteFlow' ,1); 
fluidLine = readKeyLines([casePath '\inFlow.dat'],'nblock'        ,1);
solidLine = readKeyLines([casePath '\inFlow.dat'],'nFish'         ,1);
if (str2double(deltaLine{1}) ~= str2double(deltaLine{2}))
    error('solid writing interval does not equal fluid writing interval!')
end
%% Read key parameters
UVW    = [str2double(UVWLine{1}) str2double(UVWLine{2}) str2double(UVWLine{3})]; % inflow velocity
Uref   = str2double(UrefLine{3}); % reference velocity
Tref   = str2double(TrefLine{3}); % reference time
Lref   = str2double(LrefLine{3}); % reference length
sTime  = str2double(TimeLine{1}); % begine time for writing 
eTime  = str2double(TimeLine{2}); % ending time for writing 
dTime  = str2double(deltaLine{1}); % writing intervl 
nBlock = floor(str2double(fluidLine{1})); % the number of fluid blocks
nSolid = floor(str2double(solidLine{1})); % the number of solids
meshContain = cell(nBlock); % the contian relationship of the fluid blocks, the space means no son block
for i = 1:nBlock
    blockLine = readKeyLines([casePath '\check.dat'],'sonBlocks',i);
    if length(blockLine) == 2
        meshContain{i} = [];
    else
        for j = 1:length(blockLine)-2
        meshContain{i} = [meshContain{i} str2double(blockLine{2+j})];
        end
    end
end
%% Calculate key parameters
nfile  = (eTime - sTime) / dTime + 1;
