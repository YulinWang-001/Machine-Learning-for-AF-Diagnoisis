% Number of patients - The number of patients of corresponding dataset 
% dataset - the dataset that will be used, either normal or AF dataset
% data type - two options - 1 or 0; 
% 1 - Normal and 0 - AF
% Depending on different data type, different offset values from peak index
% will be used

function [Relative_amp] = PWave_Extraction(Number_of_patients,dataset,data_type)
P_wave_matrix = zeros(Number_of_patients,1);

if data_type == 1
    offset_P = 100; 
    % Since the ARS interval should be <= 120ms for normal people Check if
    % the peak is valid, as in, if there could potentially be
    MinH = 0.1;
else
    offset_P = 100;
    % For patients with AF, the offset should be bigger as the QRS interval
    % might be larger 
    MinH = 0.05;
end 

for m = 1 : Number_of_patients
    [~,qrs_i] = findpeaks(dataset(m,:),'MinPeakHeight',MinH,'MinPeakDistance',200);
    No_Peaks = size(qrs_i,2); 
    P_wave_AMP_temp = [];
    temp_peak_matrix_pp = [];
    
    for i = 1 : No_Peaks
        temp_peak = qrs_i(i);
        
        if qrs_i(i) > offset_P
            dataset_temp_P = dataset(m,temp_peak-offset_P:temp_peak);
            [~,P_wave] = findpeaks(dataset_temp_P); 
            P_wave_amp = dataset(m,P_wave); % Calcuate the amplititude
            P_wave_amp_ave = mean(P_wave_amp); % Find out the average amplititude 
            P_wave_AMP_temp = [P_wave_AMP_temp;P_wave_amp_ave];
            % For all the R points, get the coresponding P waves
            temp_peak_matrix_pp = [temp_peak_matrix_pp; temp_peak];
        else 
            continue 
        end
        temp_peak_ave_pp = mean(temp_peak_matrix_pp);
        % For each person, check the mean p_wave value
        P_wave_amp_pp = mean(P_wave_AMP_temp);
        
    end
    R_amp_matrix(m) = temp_peak_ave_pp;  
    Relative_amp = 10^5*abs(P_wave_matrix./R_amp_matrix);
    P_wave_matrix(m) = P_wave_amp_pp;
    Relative_amp(isnan(Relative_amp)) = 0;
    % =================================================================== %
end

end

