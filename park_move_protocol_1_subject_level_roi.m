% Subject-level analysis - ROI beta
% Ran on R2024a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/phd_study_3
%% Select dataset

clear;

% One of cbsi, hbo, hbr
hemo_measure = "cbsi";

% File locations
subjstats_file = "../Park-MOVE_fnirs_dataset_v2/mat_files/SubjStats_setup_1_" + hemo_measure + ".mat";

% Load
SubjStats = importdata(subjstats_file);

%% Set groups in demographics

demographics = nirs.createDemographicsTable(SubjStats);
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

% Get the data
gait_as_data = readtable('data/mixed_model_parameters.csv'); 

% Prepare the table
demographics.st_walk_speed = NaN(height(demographics),1);
demographics.dt_walk_speed = NaN(height(demographics),1);
demographics.dt_cost_walk_speed = NaN(height(demographics),1);
demographics.st_step_time_var = NaN(height(demographics),1);
demographics.dt_step_time_var = NaN(height(demographics),1);
demographics.dt_cost_stroop_time = NaN(height(demographics),1);
demographics.prio = NaN(height(demographics),1);

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

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
SubjStats = job.run(SubjStats);

%% ROI 

% Set up ROI
source = [NaN NaN NaN NaN NaN NaN NaN]';    
detector = [1 2 3 4 5 6 7]';
ROI_PFC = table(source,detector);

%% Write betas

ROIs = {ROI_PFC};
ROI_names = ["PFC"];

for i=1:length(ROIs)
    
    % Get ROI results for each condition
    roi_result_all = nirs.util.roiAverage(SubjStats, ROIs{i}, ROI_names(i), false, false);
    roi_result_all.SubjectID = repmat('FNPXXXX',height(roi_result_all),1);
    for idx=1:height(roi_result_all)
        subj_idx = str2num(roi_result_all.FileIdx{idx});
        roi_result_all.age(idx,:) = SubjStats(subj_idx).demographics.age;
        roi_result_all.SubjectID(idx,:) = SubjStats(subj_idx).demographics.SubjectID;
        roi_result_all.group(idx,:) = SubjStats(subj_idx).demographics.group;
        roi_result_all.st_walk_speed(idx,:) = SubjStats(subj_idx).demographics.st_walk_speed;
        roi_result_all.dt_walk_speed(idx,:) = SubjStats(subj_idx).demographics.dt_walk_speed;
        roi_result_all.st_step_time_var(idx,:) = SubjStats(subj_idx).demographics.st_step_time_var;
        roi_result_all.dt_step_time_var(idx,:) = SubjStats(subj_idx).demographics.dt_step_time_var;
        roi_result_all.dt_cost_walk_speed(idx,:) = SubjStats(subj_idx).demographics.dt_cost_walk_speed;
        roi_result_all.dt_cost_stroop_time(idx,:) = SubjStats(subj_idx).demographics.dt_cost_stroop_time;
    end
    fname = 'data/subject_level_ROI_beta_cbsi_' + ROI_names(i) + '.csv';
    writetable(roi_result_all, fname);
    
end
