clear;
clc;

%% Load the dataset of normal patient and AF patients
% ====== ATTENTION: PLEASE CHANGE THE FOLDER TO CORRESPONDING PATH ====== %
folder_normal='D:\Sensing, Perception and Neuroergonomics\Coursework2\ECG normal\ECG normal\';
folder_AF = 'D:\Sensing, Perception and Neuroergonomics\Coursework2\ECG AF\ECG AF\';
% folder_normal='H:\MATLAB\SPN CW2\ECG normal\ECG normal';
% folder_AF = 'H:\MATLAB\SPN CW2\ECG AF\ECG AF';
% Loadall is a sefl-defined function 
DS_normal_raw = Loadall(folder_normal);
DS_AF_raw = Loadall(folder_AF);
open('Loadall');
% Build a cell for users to easier access data values
% Use dataset{number}.ECG.age/sex/data to check the data values. Here,
% number one is normal, number two is AF 

%% Proprocess the data values
% Check if all the features makes sense
% Afrer matlab obsevation, there are three 'bad vales' in normal dataset
% and one 'bad value' in AF dataset
DS_normal = Preprocess(DS_normal_raw, 915);
DS_AF = Preprocess(DS_AF_raw, 1097);
open('Preprocess')
DS = {DS_normal; DS_AF};
Total_No_Normal = size(DS{1},1);
Total_No_AF = size(DS{2},1);

%% Predefine some values
startx = 1;
upperxlim = 4000;
endx = startx+upperxlim;
lead_no = 10;
sample_person = 500; % Sample person's data to be plotted

%% Noise Elimination - Low Pass Fiter 
% Use Low-Pass Filter
N = 5; % The size of Low pass filter 
kernel = (1/N)*ones(1,N);
normal_LPfiltered = zeros(Total_No_Normal, endx);
AF_LPfiltered = zeros(Total_No_AF, endx);
for m = 1:Total_No_Normal
    normal_LPfiltered(m,:) = conv(DS{1}{m}.ECG.data(lead_no,startx:endx), kernel, 'same');
end 

for m = 1:Total_No_AF
    AF_LPfiltered(m,:) = conv(DS{2}{m}.ECG.data(lead_no,startx:endx), kernel, 'same');
end

%% Plot the Raw data and data values after low-pass filter
figure
subplot(2,2,1)
plot(DS_normal{sample_person}.ECG.data(lead_no,startx:endx));
title('Normal - Raw');
subplot(2,2,3)
plot(DS_AF{sample_person}.ECG.data(lead_no,startx:endx));
title('AF - Raw')
subplot(2,2,2)
plot(normal_LPfiltered(sample_person,:))
title('Normal - LP filtered');
subplot(2,2,4)
plot(AF_LPfiltered(sample_person,:))
title('AF - LP filtered');

%% Feature Extraction - Find the P points, Rhythm(BPM) and QRS interval 
% In this section we will extract all the features inclduing te duration of
% QRS complexes (threshold 120 Binarise), P waves(no visible P waves for
% AF ---  we will first extract feature - Irregular rhythms 
% Rhythm should be measured the distance between complexes (Neighbouring R)


% After extracting multiple features from one person's data, we will only 
% choose the first five values
% According to online medical research, the QRS ranges from 35-120 ms, 
qrs_i_normal = [];
qrs_i_AF = [];

Fs = 500; % Sample frequency
% QRS_interval_normal = zeros(Total_No_Normal,1);
% QRS_interval_AF = zeros(Total_No_AF,1);
Q_AF = zeros(Total_No_AF,1);
S_AF = zeros(Total_No_AF,1);
P_wave_normal_matrix = zeros(Total_No_Normal,1);
P_wave_AF_matrix = zeros(Total_No_AF,1);
R_normal_amp_matrix = zeros(Total_No_Normal,1);
R_AF_amp_matrix = zeros(Total_No_AF,1);

% ============================= Normal Subjects ========================= %
for m = 1 : Total_No_Normal % Number of all the patients
    [~,qrs_i_normal] = findpeaks(normal_LPfiltered(m,:),'MinPeakHeight',...
        0.1,'MinPeakDistance',200);
    No_Peaks = size(qrs_i_normal,2);
    P_wave_normal_AMP_temp = [];
    temp_peak_matrix_pp = [];
    
    for i = 1 : No_Peaks % Number of peaks of one lead 
        temp_peak = qrs_i_normal(i);
        offset_P = 100; % Detect the existence of P_wave, the data domain we selected 
        % starts from highest peak minus offset and ends at highest peak
               
        % ============================== P Wave ========================= %
        % P wave is will not be counted as a feature, just to filter out
        % some 'bad' data values
        if qrs_i_normal(i) > offset_P
            dataset_temp_P = normal_LPfiltered(m,temp_peak-offset_P:temp_peak);
            [~,P_wave_normal] = findpeaks(dataset_temp_P); 
            P_wave_amp_normal = normal_LPfiltered(m,P_wave_normal); % Calcuate the amplititude
            P_wave_amp_normal_ave = mean(P_wave_amp_normal); % Find out the average amplititude 
            P_wave_normal_AMP_temp = [P_wave_normal_AMP_temp;P_wave_amp_normal_ave];
            % For all the R points, get the coresponding P waves
            temp_peak_matrix_pp = [temp_peak_matrix_pp; temp_peak];
        else 
            continue 
        end
        P_wave_amp_normal_pp = mean(P_wave_normal_AMP_temp); % For each person, check the mean p_wave value
        temp_peak_ave_pp = mean(temp_peak_matrix_pp);
    end
    P_wave_normal_matrix(m) = P_wave_amp_normal_pp;
    R_normal_amp_matrix(m) = temp_peak_ave_pp;
    Relative_amp_normal = 10^5*abs(P_wave_normal_matrix./R_normal_amp_matrix);
    Relative_amp_normal(isnan(Relative_amp_normal)) = 0;
    % =================================================================== %

end

% ===================== QRS Interval Extraction Normal ================== %
QRS_interval_normal = QRS_Interval_Extraction(Total_No_Normal,normal_LPfiltered,1);
QRS_interval_normal = QRS_interval_normal(~isnan(QRS_interval_normal)); % Features
% ======================================================================= %

% ============================ Rhythm (BPM) ============================= %
Rhythm_normal = BPM_Extraction(Total_No_Normal,normal_LPfiltered, Fs, 1);
% ======================================================================= %

% ============================ AF Subjects ============================== %
for m = 1 : Total_No_AF
    [~,qrs_i_AF] = findpeaks(AF_LPfiltered(m,:),'MinPeakHeight',0.05,'MinPeakDistance',200);
    No_Peaks = size(qrs_i_AF,2);
    
    P_wave_AF_AMP_temp = [];
    temp_peak_matrix_pp = [];
    
    for i = 1 : No_Peaks % Number of peaks of one lead 
        temp_peak = qrs_i_AF(i);
       % ===================================== P Wave =================== %
       if qrs_i_AF(i) > offset_P
            dataset_temp_P = AF_LPfiltered(m,temp_peak-offset_P:temp_peak);
            [~,P_wave_AF] = findpeaks(dataset_temp_P); 
            P_wave_amp_AF = AF_LPfiltered(m,P_wave_AF); % Calcuate the amplititude
            P_wave_amp_AF_ave = mean(P_wave_amp_AF); % Find out the average amplititude 
            P_wave_AF_AMP_temp = [P_wave_AF_AMP_temp;P_wave_amp_AF_ave];
            % For all the R points, get the coresponding P waves
            temp_peak_matrix_pp = [temp_peak_matrix_pp; temp_peak];
        else 
            continue 
       end
        temp_peak_ave_pp = mean(temp_peak_matrix_pp);
        P_wave_amp_AF_pp = mean(P_wave_AF_AMP_temp); % For each person, check the mean p_wave value

    end
    R_AF_amp_matrix(m) = temp_peak_ave_pp;
    Relative_amp_AF = 10^5*abs(P_wave_AF_matrix./R_AF_amp_matrix);
    P_wave_AF_matrix(m) = P_wave_amp_AF_pp;
    Relative_amp_AF(isnan(Relative_amp_AF)) = 0;
    % =================================================================== %
end  

% ============================ QRS Interval Extraction AF =============== %
QRS_interval_AF = QRS_Interval_Extraction(Total_No_AF,AF_LPfiltered,0);
QRS_interval_AF = QRS_interval_AF(~isnan(QRS_interval_AF)); % Features
% ======================================================================= %

% ============================ Rhythm (BPM) AF ========================= %
Rhythm_AF = BPM_Extraction(Total_No_AF,AF_LPfiltered, Fs, 0);
% =================================================================== % 

% No zero values in both Rhythm_normal and Rhythm_AF
zeros_row_QRS_normal = find(QRS_interval_normal==0);
zeros_row_QRS_AF = find(QRS_interval_AF==0);
zeros_row_Rel_amp_normal = find(Relative_amp_normal==0);
zeros_row_Rel_amp_AF = find(Relative_amp_AF==0); % This indicates the P wave existence
all_normal = [zeros_row_QRS_normal;zeros_row_Rel_amp_normal];
all_AF = [zeros_row_QRS_AF;zeros_row_Rel_amp_AF];

% =============== All the features extracted and filtered =============== %
% Some values are zeros after feature extraction, so we will need to get
% rid of all the zeros 
Rhythm_normal(all_normal) = [];
Rhythm_AF(all_AF) = [];
QRS_interval_normal(all_normal) = [];
QRS_interval_AF(all_AF) = [];
Relative_amp_normal(all_normal) = [];
Relative_amp_AF(all_AF) = [];

% ======================================================================= %
Rhythm = [Rhythm_normal;Rhythm_AF];
QRS_interval = [QRS_interval_normal;QRS_interval_AF];
Relative_amp = [Relative_amp_normal;Relative_amp_AF];
Real_outcomes = [ones(size(Rhythm_normal,1),1);zeros(size(Rhythm_AF,1),1)];

% =================== End of Feature Extraction ========================= %
%% Classification Training and Cross-Validation
% Normal - 1; AF - 0
Matrix = [Rhythm,QRS_interval,Relative_amp,Real_outcomes];
% Permutation first
idx = randperm(size(Matrix,1));
Matrix = Matrix(idx,:);
ECG_features = Matrix(:,1:2);
Health_condition = Matrix(:,4);
figure;
gscatter(ECG_features(:,1),ECG_features(:,2),Health_condition,'rb','..');
xlim([60,135]);
ylim([10,135]);
xlabel('BPM');
ylabel('QRS Interval');
legend('AF','Normal')
fprintf('The table below illustrates the percentage of Nornal and AF \n')
tabulate(Health_condition); % percentage of normal and AF 
% Check the percentage of 1 (Normal) and 0 (AF)
prior = [0.53 0.47];

% Call the trainClassifier function, which contains classifier and k-fold
% cross-validation
[trainedClassifier, validationAccuracy,classificationNaiveBayes] = trainClassifier(ECG_features, Health_condition);
% open('trainClassifier');
fprintf('The accuracy obtained after applying k-fold cross-validation (k=10) is %.3f%%\n',100*validationAccuracy);

%% Test the result - Training set = Test Set - Train performance
labels = predict(classificationNaiveBayes,ECG_features);
perf=sum(labels==Health_condition)/size(labels,1); % performance in the range of 0 to 1
fprintf('The Accuracy of classification model without using k-fold cross-validation is %f%%\n',100*perf);

%% 2D Plot the classification using K-Fold Cross-Validation
[x_range, y_range] = meshgrid(60:1:140,0:1:120);
XGrid = [x_range(:) y_range(:)];
[~,Posterior,~] = predict(classificationNaiveBayes,XGrid);
sz = size(x_range);
s = max(Posterior,[],2);

figure
hold on
surf(x_range,y_range,reshape(Posterior(:,1),sz),'EdgeColor','none')
surf(x_range,y_range,reshape(Posterior(:,2),sz),'EdgeColor','none')
xlabel('BPM');
ylabel('QRS Interval');
title('Classification Result')
colorbar
hold off

%% 3D Plot of Classification 
figure
hold on
surf(x_range,y_range,reshape(Posterior(:,1),sz),'FaceColor','red','EdgeColor','none')
surf(x_range,y_range,reshape(Posterior(:,2),sz),'FaceColor','blue','EdgeColor','none')
xlabel('BPM');
ylabel('QRS Interval');
zlabel('Probability');
% legend(labels)
title('Classification Probability')
alpha(0.2); % Opacity
view(3)
hold off
