function [DS_AP] = Preprocess(DS_Raw, size)
k = 1;
DS_AP = cell(size,1);
for i = 1:length(DS_Raw)
    if DS_Raw{i}.ECG.age > 0 
        DS_AP{k} = DS_Raw{i};
        k = k + 1;
    end 
end  


