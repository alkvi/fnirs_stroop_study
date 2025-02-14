% Subject-level analysis
% Ran on R2024a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/phd_study_3
%% Load

original_dir = pwd();
cd 'C:\Users\alkvis\OneDrive - Karolinska Institutet\Dokument\Project\fNIRS_projects\Park-MOVE_fnirs_dataset_v2\fNIRS_data';

% Load the NIRx probe used
probe_folder = "nirx_probe";
nirx_probe = nirs.io.loadNIRxProbe(probe_folder,true);

% Load BIDS
my_data_dir = 'bids_dataset_snirf';
raw_data = nirs.bids.loadBIDS(my_data_dir, true, nirx_probe);

% Go back to start
cd(original_dir);

% Also save the data to file
%save('data/mat_files/raw_data.mat','raw_data');

% ..or load
%load('data/mat_files/raw_data.mat');

demographics = nirs.createDemographicsTable(raw_data);
demographics.SubjectID = erase(demographics.SubjectID, "sub-");
disp(demographics);

%% Add age to demographics

age_data = readtable('../Park-MOVE_fnirs_dataset_v2/basic_demographics.csv'); 
demographics.age = NaN(height(demographics),1);

% Add age data
for idx=1:height(age_data)
    subj_id_seek = string(age_data.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).age = ones(sum(match_idx),1) * age_data(idx,:).age;
end

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
raw_data = job.run(raw_data);

%% Prune bad channels in data

bad_ch_file = "../Park-MOVE_fnirs_dataset_v2/temp_data/sci_bad_channels.csv";
bad_ch_table = readtable(bad_ch_file);

% Loop through bad channel table for each src-det pair.
for i = 1:size(bad_ch_table,1)
    link = bad_ch_table(i,:);
    
    % Find which subjects have this pair marked as bad.
    subject_list = string(bad_ch_table(i,:).subject).split(',');
    if strcmp(subject_list,"None")
        continue
    end
    
    % raw_data(j).data is TxN where N is number of channels
    % replace raw_data(j).data with NaN in col N for bad channel in N
    for j = 1:length(subject_list)
        subj_parts = subject_list(j).split('-');
        subj = subj_parts(1);
        protocol = subj_parts(2);
        subj_idx = find(strcmp(demographics.SubjectID,subj) & strcmp(demographics.session,protocol));
        idx_to_null = raw_data(subj_idx).probe.link.source == link.src ...
            & raw_data(subj_idx).probe.link.detector == link.det;
        
        fprintf("Nulling src-det %d-%d for subj %s\n", link.src, link.det, subj);
        raw_data(subj_idx).data(:,idx_to_null) = NaN;
    end
end

%% Preprocess 

% Label short-sep
job = nirs.modules.LabelShortSeperation();
job.max_distance = 8;
raw_data = job.run(raw_data);

% Extra short-sep labelling:
% Some short-detectors (virtual detector 8-15) make up links 
% with not just the source they're attached to, but ones further away.
% these are still short-separation (so e.g. 1-14 and 1-15 detect the same
% signal). Mark these also. NOTE: the column is called "ShortSeperation" (sic).
for i = 1:numel(raw_data)
    short_separation_idx =(raw_data(i).probe.link.detector >= 8);
    raw_data(i).probe.link.ShortSeperation(short_separation_idx) = 1;
end

% raw_data(4) had an extra onset.
raw_data(4).stimulus('Straight_walking').onset(1) = [];

% Trim baseline
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
raw_data = job.run(raw_data);

% Set each duration to 20 seconds, make sure amp is 1.
raw_data = nirs.design.change_stimulus_duration(raw_data,[],20);
for subj_idx=1:length(raw_data)
    for key_idx=1:length(raw_data(subj_idx).stimulus.keys)
        key = raw_data(subj_idx).stimulus.keys{key_idx};
        raw_data(subj_idx).stimulus(key).amp(:) = 1;
    end
end

job = nirs.modules.OpticalDensity();
od = job.run(raw_data);

job = nirs.modules.BeerLambertLaw();
hb = job.run(od);

%% Visualize

nirs.viz.TimeSeriesViewer(hb);  

demographics = nirs.createDemographicsTable(hb);
disp(demographics);

figure;
raw_data(1).probe.defaultdrawfcn='3D mesh';
raw_data(1).probe.draw;

%% Get combined measure

% Add CBSI
job = nirs.modules.CalculateCBSI();
hb = job.run(hb);

job=nirs.modules.KeepTypes;
job.types={'cbsi', 'hbo', 'hbr'};
hb=job.run(hb);

%% Run on protocol 1

p1_idx = strcmp(demographics.session, "protocol1");
hb_setup1 = hb(p1_idx);

% Skip 117 which has only a few seconds of data
% Exclude also 83, 103 who did not stand still during rest.
hb_setup1([117,83,103]) = [];

nirs.viz.TimeSeriesViewer(hb_setup1);  

%% Then run GLM 

% Only CBSI
hb_setup1_cbsi = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'cbsi'};
hb_setup1_cbsi=job.run(hb_setup1_cbsi);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_cbsi = job.run(hb_setup1_cbsi); 
save('../Park-MOVE_fnirs_dataset_v2/mat_files/SubjStats_setup_1_cbsi.mat','SubjStats_cbsi');
clear SubjStats_cbsi;
clear hb_setup1_cbsi;

% Only HbO
hb_setup1_hbo = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'hbo'};
hb_setup1_hbo=job.run(hb_setup1_hbo);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_hbo = job.run(hb_setup1_hbo); 
save('../Park-MOVE_fnirs_dataset_v2/mat_files/SubjStats_setup_1_hbo.mat','SubjStats_hbo');
clear SubjStats_hbo;
clear hb_setup1_hbo;

% Only HHb
hb_setup1_hbr = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'hbr'};
hb_setup1_hbr=job.run(hb_setup1_hbr);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_hbr = job.run(hb_setup1_hbr); 
save('../Park-MOVE_fnirs_dataset_v2/mat_files/SubjStats_setup_1_hbr.mat','SubjStats_hbr');
clear SubjStats_hbr;
clear hb_setup1_hbr;
