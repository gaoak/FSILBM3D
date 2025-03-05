function [keyLine] = readKeyLines(filePath,keyWord,nLines)
% Open file
fileID = fopen(filePath,'r');
if fileID == -1
    error('Cant found file : %s',filePath);
end
ntimes = 0; % The number of times the keyword is read
found = false;
while ~feof(fileID)
    keyLine = fgetl(fileID);
    % check the keyword
    if contains(keyLine, keyWord)
        ntimes = ntimes + 1;
        if ntimes == nLines
            found = true;
            break;
        end
    end
end
if ~found
    error('Can not find keyword "%s" for %d times in file %s!', keyWord, nLines, filePath);
end
% Replace D with E
keyLine = strrep(keyLine, 'D', 'E');
% Remove whitespace characters at the beginning and end
keyLine = strtrim(keyLine);
% Split strings by spaces
keyLine = strsplit(keyLine, ' ');
fclose all;
end