function [FolderList] = FindFolders()
%Find .mat file to accumulate data
%   자세한 설명 위치
CurrentDir = dir;

FolderList = {};
k = 0;
for i = 1:length(CurrentDir(:,1))
    if strcmp(CurrentDir(i).name,'.') == 0 && strcmp(CurrentDir(i).name,'..') == 0 && CurrentDir(i).isdir == 1
        k = k+1;
        FolderList{k,1} = CurrentDir(i).name;
    end
end
end

