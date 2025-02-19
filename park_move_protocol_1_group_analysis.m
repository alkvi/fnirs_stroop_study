% Group-level analysis
% Ran on R2024a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/phd_study_3
%% Select dataset

clear;

% One of cbsi, hbo, hbr
hemo_measure = "cbsi";

% File locations
subjstats_file = "../Park-MOVE_fnirs_dataset_v2/mat_files/SubjStats_setup_1_" + hemo_measure + ".mat";
results_file = "data/results_table_" + hemo_measure + "_protocol_1.csv";
results_file_contrast = "data/results_table_" + hemo_measure + "_protocol_1_contrast.csv";
folder_figures = "figures";
folder_contrast = "figures";

% Load
SubjStats = importdata(subjstats_file);

%% Set groups in demographics

demographics = nirs.createDemographicsTable(SubjStats);
demographics.SubjectID = erase(demographics.SubjectID, "sub-");
disp(demographics);

ya_ids = readtable("../Park-MOVE_fnirs_dataset_v2/identifiers_YA.csv");
oa_ids = readtable("../Park-MOVE_fnirs_dataset_v2/identifiers_OA.csv");
pd_ids = readtable("../Park-MOVE_fnirs_dataset_v2/identifiers_PD.csv");
subj_ids = cell2table(demographics.SubjectID, "VariableNames", ["id_nummer"]);

ya_idx = ismember(subj_ids, ya_ids);
oa_idx = ismember(subj_ids, oa_ids);
pd_idx = ismember(subj_ids, pd_ids);

demographics.group = repmat("NA", height(demographics), 1);
demographics.group(ya_idx) = "YA";
demographics.group(oa_idx) = "OA";
demographics.group(pd_idx) = "PD";

job=nirs.modules.AddDemographics;
job.varToMatch = 'UUID';
job.demoTable=demographics;
SubjStats=job.run(SubjStats);

demographics = nirs.createDemographicsTable(SubjStats);
disp(demographics);

fprintf('Sum YA: %d\n', sum(ya_idx));
fprintf('Sum OA: %d\n', sum(oa_idx));
fprintf('Sum PD: %d\n', sum(pd_idx));

%% Add covariates

demographics = nirs.createDemographicsTable(SubjStats);

% Get the demographics
balance_data = readtable('../Park-MOVE_fnirs_dataset_v2/REDcap_data/MiniBEST_data.csv'); 
w12_data = readtable('../Park-MOVE_fnirs_dataset_v2/REDcap_data/Walk12_data.csv'); 
hads_data = readtable('../Park-MOVE_fnirs_dataset_v2/REDcap_data/HADS_data.csv'); 
updrs_data = readtable('../Park-MOVE_fnirs_dataset_v2/REDcap_data/UPDRS_data.csv'); 
neuropsych_data = readtable('../Park-MOVE_fnirs_dataset_v2/REDcap_data/Neuropsychological_data.csv'); 
redcap_data = readtable('../Park-MOVE_fnirs_dataset_v2/REDcap_data/All_REDcap_data.csv'); 
measurement_data = readtable('../Park-MOVE_fnirs_dataset_v2/measurement_dates.csv'); 
gait_as_data = readtable('data/mixed_model_parameters.csv'); 

% Prepare the table
demographics.hy = repmat("NA", height(demographics), 1);
demographics.updrs_3_motor = NaN(height(demographics),1);
demographics.balance = NaN(height(demographics),1);
demographics.w12 = NaN(height(demographics),1);
demographics.hads_anxiety = NaN(height(demographics),1);
demographics.st_walk_speed = NaN(height(demographics),1);
demographics.dt_walk_speed = NaN(height(demographics),1);
demographics.dt_cost_walk_speed = NaN(height(demographics),1);
demographics.st_step_time_var = NaN(height(demographics),1);
demographics.dt_step_time_var = NaN(height(demographics),1);
demographics.tmt2 = NaN(height(demographics),1);
demographics.tmt4 = NaN(height(demographics),1);
demographics.ledd = NaN(height(demographics),1);
demographics.dt_cost_stroop_time = NaN(height(demographics),1);
demographics.prio = NaN(height(demographics),1);
demographics.edu = NaN(height(demographics),1);
demographics.disease_dur = NaN(height(demographics),1);

% Add balance data
for idx=1:height(balance_data)
    subj_id_seek = string(balance_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).balance = balance_data(idx,:).mb_total;
end

% Add TMT data
for idx=1:height(neuropsych_data)
    subj_id_seek = string(neuropsych_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).tmt2 = neuropsych_data(idx,:).tmt_2;
    demographics(match_idx,:).tmt4 = neuropsych_data(idx,:).tmt_4;
end

% UPDRS3 motor score and HY
for idx=1:height(updrs_data)
    subj_id_seek = string(updrs_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    updrs_3_scores = table2array(updrs_data(idx,33:65));
    if sum(isnan(updrs_3_scores)) == length(updrs_3_scores)
        continue
    end
    updrs_3 = sum(updrs_3_scores, "omitnan");
    demographics(match_idx,:).updrs_3_motor = updrs_3;
    hy_value = updrs_data(idx,:).mdsupdrs3_hy;
    demographics(match_idx,:).hy(hy_value == 1 | hy_value == 2) = "HY_1_2";
    demographics(match_idx,:).hy(hy_value == 3 | hy_value == 4) = "HY_3_4";
end

% Add LEDD data and education
for idx=1:height(redcap_data)
    subj_id_seek = string(redcap_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).ledd = redcap_data(idx,:).led_total;
    demographics(match_idx,:).edu = redcap_data(idx,:).crf_utbildning_ar;
end

% Disease duration
for idx=1:height(measurement_data)
    subj_id_seek = string(measurement_data.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    match_idx_redcap = strcmp(redcap_data.id_nummer, subj_id_seek);
    if sum(match_idx_redcap) < 1
        continue
    end
    measure_year = measurement_data(idx,:).measurement_date_t1.Year;
    demographics(match_idx,:).disease_dur = measure_year - redcap_data(match_idx_redcap,:).crf_pd_year_phone;
end

% Gait data
for idx=1:height(gait_as_data)
    subj_id_seek = string(gait_as_data.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).st_walk_speed = gait_as_data(idx,:).st_walk_speed;
    demographics(match_idx,:).dt_walk_speed = gait_as_data(idx,:).dt_walk_speed;
    demographics(match_idx,:).dt_cost_walk_speed = gait_as_data(idx,:).dt_cost_walk_speed;
    demographics(match_idx,:).st_step_time_var = gait_as_data(idx,:).st_step_time_var;
    demographics(match_idx,:).dt_step_time_var = gait_as_data(idx,:).dt_step_time_var;
    demographics(match_idx,:).dt_cost_stroop_time = gait_as_data(idx,:).dt_cost_stroop_time;
    demographics(match_idx,:).prio = demographics(match_idx,:).dt_cost_walk_speed - demographics(match_idx,:).dt_cost_stroop_time;
end

% Normalize covariates for comparable betas
demographics(:,37:52) = normalize(demographics(:,37:52), 'zscore');
demographics(:,14) = normalize(demographics(:,14), 'zscore');

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
SubjStats = job.run(SubjStats);

%% Select groups

demographics = nirs.createDemographicsTable(SubjStats);

ya_idx = strcmpi(demographics.group,'YA');
oa_idx = strcmpi(demographics.group,'OA');
pd_idx = strcmpi(demographics.group,'PD');

% Keep all subjects
nan_idx = logical(zeros(height(demographics),1));

% Only YA
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(ya_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_ya = SubjStats(selected_idx);
demographics_ya = nirs.createDemographicsTable(SubjStats_ya);

% Only OA
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(oa_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_oa = SubjStats(selected_idx);
demographics_oa = nirs.createDemographicsTable(SubjStats_oa);

% Only PD
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(pd_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_pd = SubjStats(selected_idx);
demographics_pd = nirs.createDemographicsTable(SubjStats_pd);

fprintf('Sum YA: %d\n', length(SubjStats_ya));
fprintf('Sum OA: %d\n', length(SubjStats_oa));
fprintf('Sum PD: %d\n', length(SubjStats_pd));

%% Set up ROI

% Set up ROI
source = [NaN NaN NaN NaN NaN NaN NaN]';    
detector = [1 2 3 4 5 6 7]';
ROI_PFC = table(source,detector);

%% Outlier removal - optional

% This is how outlier removal can be performed based 
% on the leverage on the group model per subject.
% Here, it will remove 2 from OA and 1 from PD.
% SubjStats_oa_orig = SubjStats_oa;
% SubjStats_pd_orig = SubjStats_pd;
% job=nirs.modules.RemoveOutlierSubjects;
% SubjStats_oa = job.run(SubjStats_oa);
% job=nirs.modules.RemoveOutlierSubjects;
% SubjStats_pd = job.run(SubjStats_pd);

%% Aim 1

% st_step_time_var
formula_oa = 'beta ~ -1 + cond + st_step_time_var + age';
formula_pd = 'beta ~ -1 + cond + st_step_time_var + age + updrs_3_motor';

% OA group
job = nirs.modules.MixedEffects();
job.formula = formula_oa;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_oa = job.run(SubjStats_oa);

%  PD group
job = nirs.modules.MixedEffects();
job.formula = formula_pd;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd = job.run(SubjStats_pd);

% Collect results
roi_result_oa_aim1 = nirs.util.roiAverage(GroupStats_oa, ROI_PFC, 'PFC');
roi_result_pd_aim1 = nirs.util.roiAverage(GroupStats_pd, ROI_PFC, 'PFC');

% Diagnostics and validation
covar = "step_time_var";
[AIC, BIC] = PlotDiagnostics(roi_result_oa_aim1.model{1}, "OA_model_1" , covar);
roi_result_oa_aim1.formula = repmat(string(formula_oa), size(roi_result_oa_aim1,1),1);
roi_result_oa_aim1.AIC = repmat(AIC, size(roi_result_oa_aim1,1),1);
roi_result_oa_aim1.BIC = repmat(AIC, size(roi_result_oa_aim1,1),1);
roi_result_oa_aim1.group = repmat("OA", size(roi_result_oa_aim1,1),1);
roi_result_oa_aim1.comment = repmat("Aim 1", size(roi_result_oa_aim1,1),1);
roi_result_oa_aim1.included_n = repmat(GroupStats_oa.demographics.included_subjects_n, size(roi_result_oa_aim1,1),1);
nan_idx = isnan(demographics_oa.st_step_time_var) | isnan(demographics_oa.age);
roi_result_oa_aim1.included_n = roi_result_oa_aim1.included_n - sum(nan_idx);

[AIC, BIC] = PlotDiagnostics(roi_result_pd_aim1.model{1}, "PD_model_1" , covar);
roi_result_pd_aim1.formula = repmat(string(formula_pd), size(roi_result_pd_aim1,1),1);
roi_result_pd_aim1.AIC = repmat(AIC, size(roi_result_pd_aim1,1),1);
roi_result_pd_aim1.BIC = repmat(AIC, size(roi_result_pd_aim1,1),1);
roi_result_pd_aim1.group = repmat("PD", size(roi_result_pd_aim1,1),1);
roi_result_pd_aim1.comment = repmat("Aim 1", size(roi_result_pd_aim1,1),1);
roi_result_pd_aim1.comment = repmat("Aim 1", size(roi_result_pd_aim1,1),1);
roi_result_pd_aim1.included_n = repmat(GroupStats_pd.demographics.included_subjects_n, size(roi_result_pd_aim1,1),1);
nan_idx = isnan(demographics_pd.st_step_time_var) | isnan(demographics_pd.age) | isnan(demographics_pd.updrs_3_motor);
roi_result_pd_aim1.included_n = roi_result_pd_aim1.included_n - sum(nan_idx);

%% Aim 2 

% dt_cost_walk_speed
% dt_cost_stroop_time
% prio (= dt_cost_walk_speed - dt_cost_stroop_time)

formula_oa = 'beta ~ -1 + cond + age + dt_cost_walk_speed';
formula_pd = 'beta ~ -1 + cond + age + updrs_3_motor + dt_cost_walk_speed';

% OA group
job = nirs.modules.MixedEffects();
job.formula = formula_oa;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_oa = job.run(SubjStats_oa);

%  PD group
job = nirs.modules.MixedEffects();
job.formula = formula_pd;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd = job.run(SubjStats_pd);

% Collect results
roi_result_oa_aim2_1 = nirs.util.roiAverage(GroupStats_oa, ROI_PFC, 'PFC');
roi_result_pd_aim2_1 = nirs.util.roiAverage(GroupStats_pd, ROI_PFC, 'PFC');

% Diagnostics and validation
covar = "cost_walk_speed";
formula_oa = 'beta ~ -1 + cond + age + dt_cost_walk_speed';
formula_pd = 'beta ~ -1 + cond + age + updrs_3_motor + dt_cost_walk_speed';
[AIC, BIC] = PlotDiagnostics(roi_result_oa_aim2_1.model{1}, "OA_model_2_1" , covar);
roi_result_oa_aim2_1.formula = repmat(string(formula_oa), size(roi_result_oa_aim2_1,1),1);
roi_result_oa_aim2_1.AIC = repmat(AIC, size(roi_result_oa_aim2_1,1),1);
roi_result_oa_aim2_1.BIC = repmat(AIC, size(roi_result_oa_aim2_1,1),1);
roi_result_oa_aim2_1.group = repmat("OA", size(roi_result_oa_aim2_1,1),1);
roi_result_oa_aim2_1.comment = repmat("Aim 2_1", size(roi_result_oa_aim2_1,1),1);
roi_result_oa_aim2_1.included_n = repmat(GroupStats_oa.demographics.included_subjects_n, size(roi_result_oa_aim2_1,1),1);
nan_idx = isnan(demographics_oa.age) | isnan(demographics_oa.dt_cost_walk_speed);
roi_result_oa_aim2_1.included_n = roi_result_oa_aim2_1.included_n - sum(nan_idx);

[AIC, BIC] = PlotDiagnostics(roi_result_pd_aim2_1.model{1}, "PD_model_2_1" , covar);
roi_result_pd_aim2_1.formula = repmat(string(formula_pd), size(roi_result_pd_aim2_1,1),1);
roi_result_pd_aim2_1.AIC = repmat(AIC, size(roi_result_pd_aim2_1,1),1);
roi_result_pd_aim2_1.BIC = repmat(AIC, size(roi_result_pd_aim2_1,1),1);
roi_result_pd_aim2_1.group = repmat("PD", size(roi_result_pd_aim2_1,1),1);
roi_result_pd_aim2_1.comment = repmat("Aim 2_1", size(roi_result_pd_aim2_1,1),1);
roi_result_pd_aim2_1.included_n = repmat(GroupStats_pd.demographics.included_subjects_n, size(roi_result_pd_aim2_1,1),1);
nan_idx = isnan(demographics_pd.age) | isnan(demographics_pd.dt_cost_walk_speed) | isnan(demographics_pd.updrs_3_motor);
roi_result_pd_aim2_1.included_n = roi_result_pd_aim2_1.included_n - sum(nan_idx);

formula_oa = 'beta ~ -1 + cond + age + dt_cost_stroop_time';
formula_pd = 'beta ~ -1 + cond + age + updrs_3_motor + dt_cost_stroop_time';

% OA group
job = nirs.modules.MixedEffects();
job.formula = formula_oa;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_oa = job.run(SubjStats_oa);

%  PD group
job = nirs.modules.MixedEffects();
job.formula = formula_pd;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd = job.run(SubjStats_pd);

% Collect results
roi_result_oa_aim2_2= nirs.util.roiAverage(GroupStats_oa, ROI_PFC, 'PFC');
roi_result_pd_aim2_2 = nirs.util.roiAverage(GroupStats_pd, ROI_PFC, 'PFC');

covar = "cost_stroop_time";
formula_oa = 'beta ~ -1 + cond + age + dt_cost_stroop_time';
formula_pd = 'beta ~ -1 + cond + age + updrs_3_motor + dt_cost_stroop_time';
[AIC, BIC] = PlotDiagnostics(roi_result_oa_aim2_2.model{1}, "OA_model_2_2" , covar);
roi_result_oa_aim2_2.formula = repmat(string(formula_oa), size(roi_result_oa_aim2_2,1),1);
roi_result_oa_aim2_2.AIC = repmat(AIC, size(roi_result_oa_aim2_2,1),1);
roi_result_oa_aim2_2.BIC = repmat(AIC, size(roi_result_oa_aim2_2,1),1);
roi_result_oa_aim2_2.group = repmat("OA", size(roi_result_oa_aim2_2,1),1);
roi_result_oa_aim2_2.comment = repmat("Aim 2_2", size(roi_result_oa_aim2_2,1),1);
roi_result_oa_aim2_2.included_n = repmat(GroupStats_oa.demographics.included_subjects_n, size(roi_result_oa_aim2_2,1),1);
nan_idx = isnan(demographics_oa.age) | isnan(demographics_oa.dt_cost_stroop_time);
roi_result_oa_aim2_2.included_n = roi_result_oa_aim2_2.included_n - sum(nan_idx);

[AIC, BIC] = PlotDiagnostics(roi_result_pd_aim2_2.model{1}, "PD_model_2_2" , covar);
roi_result_pd_aim2_2.formula = repmat(string(formula_pd), size(roi_result_pd_aim2_2,1),1);
roi_result_pd_aim2_2.AIC = repmat(AIC, size(roi_result_pd_aim2_2,1),1);
roi_result_pd_aim2_2.BIC = repmat(AIC, size(roi_result_pd_aim2_2,1),1);
roi_result_pd_aim2_2.group = repmat("PD", size(roi_result_pd_aim2_2,1),1);
roi_result_pd_aim2_2.comment = repmat("Aim 2_2", size(roi_result_pd_aim2_2,1),1);
roi_result_pd_aim2_2.included_n = repmat(GroupStats_pd.demographics.included_subjects_n, size(roi_result_pd_aim2_2,1),1);
nan_idx = isnan(demographics_pd.age) | isnan(demographics_pd.dt_cost_stroop_time) | isnan(demographics_pd.updrs_3_motor);
roi_result_pd_aim2_2.included_n = roi_result_pd_aim2_2.included_n - sum(nan_idx);

formula_oa = 'beta ~ -1 + cond + age + prio';
formula_pd = 'beta ~ -1 + cond + age + updrs_3_motor + prio';

% OA group
job = nirs.modules.MixedEffects();
job.formula = formula_oa;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_oa = job.run(SubjStats_oa);

%  PD group
job = nirs.modules.MixedEffects();
job.formula = formula_pd;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd = job.run(SubjStats_pd);

% Collect results
roi_result_oa_aim2_3= nirs.util.roiAverage(GroupStats_oa, ROI_PFC, 'PFC');
roi_result_pd_aim2_3 = nirs.util.roiAverage(GroupStats_pd, ROI_PFC, 'PFC');

covar = "prio";
formula_oa = 'beta ~ -1 + cond + age + prio';
formula_pd = 'beta ~ -1 + cond + age + updrs_3_motor + prio';
[AIC, BIC] = PlotDiagnostics(roi_result_oa_aim2_3.model{1}, "OA_model_2_3" , covar);
roi_result_oa_aim2_3.formula = repmat(string(formula_oa), size(roi_result_oa_aim2_3,1),1);
roi_result_oa_aim2_3.AIC = repmat(AIC, size(roi_result_oa_aim2_3,1),1);
roi_result_oa_aim2_3.BIC = repmat(AIC, size(roi_result_oa_aim2_3,1),1);
roi_result_oa_aim2_3.group = repmat("OA", size(roi_result_oa_aim2_3,1),1);
roi_result_oa_aim2_3.comment = repmat("Aim 2_3", size(roi_result_oa_aim2_3,1),1);
roi_result_oa_aim2_3.included_n = repmat(GroupStats_oa.demographics.included_subjects_n, size(roi_result_oa_aim2_3,1),1);
nan_idx = isnan(demographics_oa.age) | isnan(demographics_oa.prio);
roi_result_oa_aim2_3.included_n = roi_result_oa_aim2_3.included_n - sum(nan_idx);

[AIC, BIC] = PlotDiagnostics(roi_result_pd_aim2_3.model{1}, "PD_model_2_3" , covar);
roi_result_pd_aim2_3.formula = repmat(string(formula_pd), size(roi_result_pd_aim2_3,1),1);
roi_result_pd_aim2_3.AIC = repmat(AIC, size(roi_result_pd_aim2_3,1),1);
roi_result_pd_aim2_3.BIC = repmat(AIC, size(roi_result_pd_aim2_3,1),1);
roi_result_pd_aim2_3.group = repmat("PD", size(roi_result_pd_aim2_3,1),1);
roi_result_pd_aim2_3.comment = repmat("Aim 2_3", size(roi_result_pd_aim2_3,1),1);
roi_result_pd_aim2_3.included_n = repmat(GroupStats_pd.demographics.included_subjects_n, size(roi_result_pd_aim2_3,1),1);
nan_idx = isnan(demographics_pd.age) | isnan(demographics_pd.prio) | isnan(demographics_pd.updrs_3_motor);
roi_result_pd_aim2_3.included_n = roi_result_pd_aim2_3.included_n - sum(nan_idx);

%% Aim 3

% st_step_time_var*tmt4 

formula_oa = 'beta ~ -1 + cond + st_step_time_var*tmt4 + age + edu';
formula_pd = 'beta ~ -1 + cond + st_step_time_var*tmt4 + age + updrs_3_motor + edu';

% OA group
job = nirs.modules.MixedEffects();
job.formula = formula_oa;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_oa = job.run(SubjStats_oa);

%  PD group
job = nirs.modules.MixedEffects();
job.formula = formula_pd;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd = job.run(SubjStats_pd);

% Collect results
roi_result_oa_aim3 = nirs.util.roiAverage(GroupStats_oa, ROI_PFC, 'PFC');
roi_result_pd_aim3 = nirs.util.roiAverage(GroupStats_pd, ROI_PFC, 'PFC');

% Diagnostics and validation
covar = "step_time_var_tmt4";
[AIC, BIC] = PlotDiagnostics(roi_result_oa_aim3.model{1}, "OA_model_3" , covar);
roi_result_oa_aim3.formula = repmat(string(formula_oa), size(roi_result_oa_aim3,1),1);
roi_result_oa_aim3.AIC = repmat(AIC, size(roi_result_oa_aim3,1),1);
roi_result_oa_aim3.BIC = repmat(AIC, size(roi_result_oa_aim3,1),1);
roi_result_oa_aim3.group = repmat("OA", size(roi_result_oa_aim3,1),1);
roi_result_oa_aim3.comment = repmat("Aim 3", size(roi_result_oa_aim3,1),1);
roi_result_oa_aim3.included_n = repmat(GroupStats_oa.demographics.included_subjects_n, size(roi_result_oa_aim3,1),1);
nan_idx = isnan(demographics_oa.tmt4) | isnan(demographics_oa.st_step_time_var) ... 
    | isnan(demographics_oa.age) | isnan(demographics_oa.edu);
roi_result_oa_aim3.included_n = roi_result_oa_aim3.included_n - sum(nan_idx);

[AIC, BIC] = PlotDiagnostics(roi_result_pd_aim3.model{1}, "PD_model_3" , covar);
roi_result_pd_aim3.formula = repmat(string(formula_pd), size(roi_result_pd_aim3,1),1);
roi_result_pd_aim3.AIC = repmat(AIC, size(roi_result_pd_aim3,1),1);
roi_result_pd_aim3.BIC = repmat(AIC, size(roi_result_pd_aim3,1),1);
roi_result_pd_aim3.group = repmat("PD", size(roi_result_pd_aim3,1),1);
roi_result_pd_aim3.comment = repmat("Aim 3", size(roi_result_pd_aim3,1),1);
roi_result_pd_aim3.included_n = repmat(GroupStats_pd.demographics.included_subjects_n, size(roi_result_pd_aim3,1),1);
nan_idx = isnan(demographics_pd.updrs_3_motor)|isnan(demographics_pd.tmt4) | isnan(demographics_pd.st_step_time_var) ... 
    | isnan(demographics_pd.age) | isnan(demographics_pd.edu);
roi_result_pd_aim3.included_n = roi_result_pd_aim3.included_n - sum(nan_idx);

%% Export models for external diagnostics

models = {roi_result_oa_aim1;
    roi_result_pd_aim1;
    roi_result_oa_aim2_1;
    roi_result_pd_aim2_1;
    roi_result_oa_aim2_2;
    roi_result_pd_aim2_2;
    roi_result_oa_aim2_3;
    roi_result_pd_aim2_3;
    roi_result_oa_aim3;
    roi_result_pd_aim3};

model_names = ["model_aim1_oa";
    "model_aim1_pd";
    "model_aim2_1_oa";
    "model_aim2_1_pd";
    "model_aim2_2_oa";
    "model_aim2_2_pd";
    "model_aim2_3_oa";
    "model_aim2_3_pd";
    "model_aim3_oa";
    "model_aim3_pd";];

for model_idx=1:length(model_names)
    
    % Which model
    save_model = models{model_idx}.model{1};
    save_model_name = model_names(model_idx);
    save_model_est_name = "data/model_eval/" + save_model_name + "_estimates.csv";
    save_model_data_name = "data/model_eval/" + save_model_name + "_data.csv";
    save_model_formula_name = "data/model_eval/" + save_model_name + "_formula.txt";
    save_model_formula_description = "data/model_eval/" + save_model_name + "_description.txt";
    
    % Export model coefficients
    coeffs = save_model.Coefficients.Estimate;
    intercept = save_model.Coefficients.Estimate(1);
    predictors = save_model.PredictorNames;
    writetable(table(coeffs, predictors, 'VariableNames', {'Coefficients', 'Predictors'}), save_model_est_name);
    
    % Export residuals, fitted values, and variables
    residuals = save_model.Residuals.Raw;
    fitted_values = save_model.Fitted;
    model_data = table(residuals, fitted_values, 'VariableNames', {'Residuals', 'FittedValues'});
    model_data = [model_data, save_model.Variables];
    writetable(model_data, save_model_data_name);

    % Export formula
    file_txt = fopen(save_model_formula_name,'w');
    fprintf(file_txt,string(save_model.Formula));
    fclose(file_txt);
    file_txt = fopen(save_model_formula_description,'w');
    fprintf(file_txt,models{model_idx}.formula(1));
    fclose(file_txt);

end

%% Summarize in a table

hypothesis_table = [roi_result_oa_aim1; 
    roi_result_pd_aim1;
    roi_result_oa_aim2_1;
    roi_result_pd_aim2_1;
    roi_result_oa_aim2_2;
    roi_result_pd_aim2_2;
    roi_result_oa_aim2_3;
    roi_result_pd_aim2_3;
    roi_result_oa_aim3;
    roi_result_pd_aim3;];

% Drop unneeded columns
hypothesis_table(:,[10,11]) = [];
disp(hypothesis_table);

% Write table 
writetable(hypothesis_table, results_file)

%% Exploratory

% Select both OA and PD
selected_idx = zeros(size(SubjStats,2),1);
selected_idx(oa_idx,1) = 1;
selected_idx(pd_idx,1) = 1;
selected_idx = logical(selected_idx);
SubjStats_diff = SubjStats(selected_idx);
demographics_diff = nirs.createDemographicsTable(SubjStats_diff);

% Settings
formula_diff = 'beta ~ -1 + cond:group';

% Run group model
job = nirs.modules.MixedEffects();
job.formula = formula_diff;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_diff = job.run(SubjStats_diff);

% Draw some figures.
GroupStats_diff.probe.defaultdrawfcn='3D mesh (frontal)';
GroupStats_diff.probe = GroupStats_diff.probe.SetFiducials_Visibility(false);
GroupStats_diff.draw('tstat', [-10 10], 'q < 0.05');
GroupStats_diff.printAll('tstat', [-10 10], 'q < 0.05', folder_figures, 'png')

% Get some contrasts
c = [-1 1 0 0 0 0];
ContrastStats = GroupStats_diff.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder_contrast, 'png')
contrast_roi_1 = nirs.util.roiAverage(ContrastStats, ROI_PFC, 'PFC');

c = [0 0 -1 1 0 0];
ContrastStats = GroupStats_diff.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder_contrast, 'png')
contrast_roi_2 = nirs.util.roiAverage(ContrastStats, ROI_PFC, 'PFC');

c = [0 0 0 0 -1 1];
ContrastStats = GroupStats_diff.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder_contrast, 'png')
contrast_roi_3 = nirs.util.roiAverage(ContrastStats, ROI_PFC, 'PFC');

%% Summarize 

contrast_table = [contrast_roi_1; 
    contrast_roi_2;
    contrast_roi_3;];

% Drop unneeded columns
contrast_table(:,[10]) = [];
disp(contrast_table);

writetable(contrast_table, results_file_contrast)

%% Check assumptions

function [AIC, BIC] = PlotDiagnostics(model, titlestr, covar)

fitted_vals = model.Fitted;
residuals_raw = model.Residuals.Raw;
residuals_standardized = model.Residuals.Standardized;
leverage = model.Diagnostics.Leverage;

% residuals are (approximately) normally distributed (Q-Q plot)
figure;
subplot(3,2,1);
histfit(residuals_raw);
title('Histogram - raw residuals');
subplot(3,2,2);
qqplot(residuals_raw);
title('Q-Q - raw residuals');

% Residual vs fitted to assess possible heteroscedascity
% Plot fitted vs residuals
subplot(3,2,3);
h = plot(fitted_vals, residuals_standardized, 'bx');
xlabel("fitted");
ylabel("residual (standardized)");
hold on;

% Add 0 line and a polyline
p = polyfit(fitted_vals, residuals_standardized, 4);
px = linspace(min(fitted_vals), max(fitted_vals));
py = polyval(p, px);
plot(px, py, 'r-', 'LineWidth', 1);
xlim = [min(fitted_vals) max(fitted_vals)];
line(xlim,[0 0],'LineStyle','--', 'LineWidth', 1);
v = axis; % get current values
lo = min( v(1:2:end) ); % lower limit
up = max( v(2:2:end) ); % uppper limit
axis( [lo up lo up] );
hold off;
title('Fitted-Residual');

% Leverage vs residual
subplot(3,2,4);
plot(leverage, residuals_standardized, 'bx');
hold on;
xlabel("leverage");
ylabel("residual (standardized)");
p = polyfit(leverage, residuals_standardized, 4);
px = linspace(min(leverage), max(leverage));
py = polyval(p, px);
plot(px, py, 'r-', 'LineWidth', 1);
xlim = [min(leverage) max(leverage)];
line(xlim,[0 0],'LineStyle','--', 'LineWidth', 1);
hold off;
title('Leverage-Residual');

% Linearity in covariate
if covar ~= ""
    subplot(3,2,5);
    plot(residuals_raw, table2array(model.Variables(:,covar)), 'bx');
    title('Residual-Covariate');
end

sgtitle(titlestr) 

% Save the figure
fig = gcf;
set(gcf,'Position',[100 100 1000 700])
figname  = titlestr + "_" + covar + ".png";
figname = strrep(figname, " ", "_");
exportgraphics(fig,figname,'Resolution',300);

AIC = model.ModelCriterion.AIC;
BIC = model.ModelCriterion.BIC;

end