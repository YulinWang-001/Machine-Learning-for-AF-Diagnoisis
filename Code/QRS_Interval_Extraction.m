% Number of patients - The number of patients of corresponding dataset 
% dataset - the dataset that will be used, either normal or AF dataset
% data type - two options - 1 or 0; 
% 1 - Normal and 0 - AF
% Depending on different data type, different offset values from peak index
% will be used
function [QRS_interval] = QRS_Interval_Extraction(Number_of_patients,dataset,data_type)
[~,col] = size(dataset); % How many data values per person will have
QRS_interval = zeros(Number_of_patients,1);

if data_type == 1
    offset_QRS = 60; 
    % Since the ARS interval should be <= 120ms for normal people Check if 
    % the peak is valid, as in, if there could potentially be
    MinH = 0.1;
else
    offset_QRS = 90;
    MinH = 0.05;
    % Since the ARS could be > 120ms
    % Check if the peak is valid, as in, if there could potentially be
    % a QRS interval at this peak
end 

for m = 1 : Number_of_patients
    [~,qrs_i] = findpeaks(dataset(m,:),'MinPeakHeight',MinH,'MinPeakDistance',200);
    No_Peaks = size(qrs_i,2);
    QRS_interval_temp = 0;
    No_valid_peak = 0; % Number of Valid Peak, since some peaks are near 
    % the beginning and the end of dataset, which we will ignore
    
    for i = 1 : No_Peaks
        temp_peak = qrs_i(i);
        
        if qrs_i(i) > offset_QRS && qrs_i(i) < col-offset_QRS
            dataset_temp = -dataset(m,temp_peak-offset_QRS:temp_peak+offset_QRS);
            [~,Q_and_S] = findpeaks(dataset_temp); % Find the Q and S values
            Q_and_S = Q_and_S + temp_peak-offset_QRS;
            
            substraction_vector = temp_peak*ones(1,size(Q_and_S,2))-Q_and_S;
            positive_part = substraction_vector(substraction_vector>0);
            negative_part = substraction_vector(substraction_vector<0);
            
            if size(positive_part,2)==0 || size(negative_part,2)==0
                continue
            else 
                Q = positive_part(end);
                S = negative_part(1);
                QRS_interval_temp = QRS_interval_temp + abs(Q-S);
                No_valid_peak = No_valid_peak + 1;
            end 
            Average_QRS_interval_AF = QRS_interval_temp/No_valid_peak;
            QRS_interval(m,1) = Average_QRS_interval_AF;
        end
    end 
        
end

QRS_interval = QRS_interval(~isnan(QRS_interval)); % Features

end