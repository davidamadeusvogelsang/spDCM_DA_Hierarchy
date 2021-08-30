%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% PEB of PEBS approach 
% D. A. Vogelsang 
%%-----------------------------------------------------------------------

clear all; 
clc; 
close all; 

%% Directories
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
dcmbase = [base, 'DCM_models']; 

%output directory 
output_dir = [base, '/Second_Level/']; 

% add SPM 
SPM_path = '/home/despoC/dvogel/Toolbox/spm12_v7487/';
addpath(genpath(SPM_path)); 

%% Subject information 
subFndr = dir([base '/sub*']);
numSub = length(subFndr);

%% Experiment information 
sessions = {'bromo', 'placebo', 'tolcapone'};


% loop to create a PEB for all sessions and subjects 
for subIdx = 1:numSub
    subject = subFndr(subIdx).name; 

    for j = 1:length(sessions)
        sess = sessions(j); 
        ses = char(sess); 
        session = ['session-' ses]; 
        cd(dcmbase); 
        
        clear S1
        %load each of the models 
        S1(j) = load(['/DCM_', subject, '-', session, '-models.mat'], 'DCM');
        %put it into a cell array so spm_dcm_peb can read it 
        GCM{j,1} = S1(j).DCM{1,1}; 

        
    end

    fprintf(1, 'Done with loading DCMs \n');


    % Specify PEB model settings (see batch editor for help on each setting)
    M = struct();
    M.alpha = 1;
    M.beta  = 16;
    M.hE    = 0;
    M.hC    = 1/16;
    M.Q     = 'single';
    % Specify design matrix for N subjects. It should start with a constant column
    M.X = [1 1 1;1 -1 0;0 -1 1]; 
    M.Xnames = {'group', 'brom_vs_plac', 'tolc_vs_plac'}; 
    

    % Choose field
    field = {'A'};
    rng('default') % for reproducibility 
    cd(dcmbase); 
    
    % RUN the first PEB (so second level PEB to compare the sessions for
    % each subject 
    PEB = spm_dcm_peb(GCM,M,field);
    PEBs{subIdx,1} = PEB; %put it into a cell array 

end
clear M 

% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';
M.X = ones(numSub,1); 
M.Xnames = {'group'}; 

% Choose field
field = {'A'};

rng('default') % for reproducibility 

% ------------------------------------------------------------------------------
% Run the PEB of PEBs
% ------------------------------------------------------------------------------
PEB_group = spm_dcm_peb(PEBs, M, field);

% ------------------------------------------------------------------------------
% Run Bayesian Model Comparison 
% ------------------------------------------------------------------------------
BMA = spm_dcm_peb_bmc(PEB_group); 

% ------------------------------------------------------------------------------
% Review results
% ------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCM);

% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEBofPEBS_results.mat'],'GCM', 'PEBs','BMA', 'PEB_group'); 
    

