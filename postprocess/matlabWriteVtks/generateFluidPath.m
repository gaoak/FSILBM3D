function [readName,writeName] = generateFluidPath(readPath,writePath,nTime,id)
% Convert input numbers to strings
numStr = num2str(nTime * 1e5);

% Calculate the number of zeros that need to be filled in the first place
zeroNumbers = 10 - length(numStr);

% If zero is required
if zeroNumbers > 0
    % Make up zeros in front
    newStr = [repmat('0', 1, zeroNumbers), numStr];
else
    % If the length is sufficient, return directly
    newStr = numStr;
end
readName  = [readPath 'Flow' newStr '_b00' num2str(id)];
writeName = [writePath 'b00' num2str(id) 'Flow' newStr '.vtk' ];
end
