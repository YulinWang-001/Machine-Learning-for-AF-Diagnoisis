function [Dataset] = Loadall(folder_name)
matfiles = dir(fullfile(folder_name, '*.mat'));
Dataset = cell(length(matfiles),1);
for i = 1:length(matfiles)
    Dataset{i} = load(matfiles(i).name);
end
end

