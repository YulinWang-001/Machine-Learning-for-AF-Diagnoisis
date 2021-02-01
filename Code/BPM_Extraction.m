% Inputs - Number of Patients and the dataset, note that the data set
% contains the same row of signal from all the patients, Fs is user
% pre-defined sampling frequency
% data_type - two options --> 1 - Normal and 0 AF
function [Rhythm] = BPM_Extraction(Number_of_patients,dataset,Fs,data_type)
if data_type == 1
    MinH = 0.1;
else 
    MinH = 0.05;
end 
for m = 1:Number_of_patients
    Peak_height = max(dataset(m,:));
    [~,qrs_i] = findpeaks(dataset(m,:),'MinPeakHeight',MinH,'MinPeakDistance',200);
    for i = 1:length(qrs_i)-1
        RR_interval(i) = qrs_i(i+1) - qrs_i(i);
    end
    secperbeat = mean(RR_interval)*1/Fs; % Second per beat (Sampling at 500 Hz)
    bpm = 60/secperbeat;
    Rhythm(m,1) = bpm; % Features
end
end 

