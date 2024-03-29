function [MatFileList] = FindMatFiles()
%Find .mat file to accumulate data
%   자세한 설명 위치
CurrentDir = dir;

MatFileList = {};
k = 0;
for i = 1:length(CurrentDir(:,1))
    if length(CurrentDir(i).name) > 4
        if strcmp(CurrentDir(i).name(end-3:end),'.mat') == 1
            k = k+1;
            MatFileList{k,1} = CurrentDir(i).name;
        end
    end
end
end

