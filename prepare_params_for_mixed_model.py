import os 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == "__main__":

    pd.set_option('display.max_rows', None)

    param_file_time = "../Park-MOVE_fnirs_dataset_v2/Task_data/auditory_stroop_answer_time.csv"
    param_file_gait = "../Park-MOVE_fnirs_dataset_v2/IMU_data/imu_gait_parameters.csv"
    param_file_var = "../Park-MOVE_fnirs_dataset_v2/IMU_data/imu_variability_parameters.csv"
    
    ya_ids = pd.read_csv("../Park-MOVE_fnirs_dataset_v2/identifiers_YA.csv")['id_nummer'].to_list()
    oa_ids = pd.read_csv("../Park-MOVE_fnirs_dataset_v2/identifiers_OA.csv")['id_nummer'].to_list()
    pd_ids = pd.read_csv("../Park-MOVE_fnirs_dataset_v2/identifiers_PD.csv")['id_nummer'].to_list()

    gait_param_data = pd.read_csv(param_file_gait)
    variability_param_data = pd.read_csv(param_file_var)
    time_data = pd.read_csv(param_file_time)

    ids = sorted(list(set(gait_param_data['subject'].to_list())))
    for id in ya_ids:
        gait_param_data.loc[gait_param_data.subject == id, 'group'] = "YA"
    for id in oa_ids:
        gait_param_data.loc[gait_param_data.subject == id, 'group'] = "OA"
    for id in pd_ids:
        gait_param_data.loc[gait_param_data.subject == id, 'group'] = "PD"

    # Stroop data
    time_data = time_data[time_data['protocol'] == 'protocol_1']

    # One outlier due to audio error
    time_data = time_data[~((time_data['subject'] == "FNP1039") & (time_data['prompt_number'] == 84))]

    # Get avg answer time per subject
    stroop_time_subject_means = time_data.groupby(['subject', 'protocol', 'block_type']).agg(
        stroop_time_mean_value=('answer_time', 'mean')
    ).reset_index()

    # Make wider so we have one column per trial type
    stroop_time_subject_means = stroop_time_subject_means.pivot(
        index='subject', columns=['block_type'], values='stroop_time_mean_value'
    ).reset_index()

    # Add DT cost
    stroop_time_subject_means['dt_diff'] = stroop_time_subject_means['DT'] - stroop_time_subject_means['ST']
    stroop_time_subject_means['dt_cost_stroop_time'] = (stroop_time_subject_means['dt_diff'] / stroop_time_subject_means['ST']) * 100

    # Gait data
    gait_param_data['Walking Speed LR'] = gait_param_data[["Walking Speed R", "Walking Speed L"]].mean(axis=1)
    gait_param_data['Step Time LR'] = gait_param_data[["Step Time R", "Step Time L"]].mean(axis=1)

    # Protocol 1
    gait_data_p1 = gait_param_data[(gait_param_data['session'] == 'protocol1')]
    variability_param_data_p1 = variability_param_data[(variability_param_data['session'] == 'protocol1')]

    # Get condition variables
    trial = "Straight_walking"
    comp_trial = "Straight_walking_and_Aud_Stroop"
    avg_frames = []
    for id in ids:
        st_data = gait_data_p1[(gait_data_p1['subject'] == id) & (gait_data_p1['trial_type'] == trial)]
        dt_data = gait_data_p1[(gait_data_p1['subject'] == id) & (gait_data_p1['trial_type'] == comp_trial)]
        avg_ws_st = st_data["Walking Speed LR"].mean(axis=0)
        avg_ws_dt = dt_data["Walking Speed LR"].mean(axis=0)
        dt_diff = avg_ws_dt - avg_ws_st
        dt_cost_walk_speed = -(dt_diff / avg_ws_st) * 100 # DTE from Kelly et al 2010

        if st_data.empty:
            continue

        var_data = variability_param_data_p1[(variability_param_data_p1['subject'] == id) & 
                                                     (variability_param_data_p1['trial_type'] == trial)]
        st_step_time_var = var_data['Step Time Variability']
        if st_step_time_var.empty:
            st_step_time_var = np.NaN
        else:
            st_step_time_var = st_step_time_var.values[0]

        var_data = variability_param_data_p1[(variability_param_data_p1['subject'] == id) & 
                                                     (variability_param_data_p1['trial_type'] == comp_trial)]
        dt_step_time_var = var_data['Step Time Variability']
        if dt_step_time_var.empty:
            dt_step_time_var = np.NaN
        else:
            dt_step_time_var = dt_step_time_var.values[0]

        avg_data = {
            'subject': [id],
            'group': [st_data['group'].values[0]],
            'protocol': ['protocol1'],
            'st_walk_speed': [avg_ws_st],
            'dt_walk_speed': [avg_ws_dt],
            'dt_cost_walk_speed': [dt_cost_walk_speed],
            'st_step_time_var': [st_step_time_var],
            'dt_step_time_var': [dt_step_time_var],
        }
        avg_frame = pd.DataFrame(data=avg_data)
        avg_frames.append(avg_frame)

    # Merge with stroop and save
    avg_frame = pd.concat(avg_frames)
    avg_frame = pd.merge(avg_frame, stroop_time_subject_means[['subject', 'dt_cost_stroop_time']], on='subject', how='left')

    print(avg_frame)
    avg_frame.to_csv("data/mixed_model_parameters.csv", index=False)
