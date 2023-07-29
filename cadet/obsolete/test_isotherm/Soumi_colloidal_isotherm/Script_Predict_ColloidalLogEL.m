
%% INPUTS
%parpool(12)
% Directory path
matlab_cadet_path = 'C:\Users\Soumi\Documents\MATLAB\matlab';
experiment_folder = 'soumitra\MAbModel_NewCol_MAR21\';

% CADET INPUTS                                              ______________
% Model Parameters are:
%               ke              Bppd            s-diffusion     
initial =  [5  2     log10(2E-13)];


% System Inputs
mol_wt = 146; % Molecular weight of the protein under analysis, kDa


%% SCRIPT START                                             ______________

% Collect csv files and list them
list = dir('*.csv');
listsize = length(list);

for i = 12:12

% Read file
InputFile = ...
 readmatrix(...
    fullfile(...
            matlab_cadet_path, ...
            experiment_folder, ...
            list(i).name));

% Protein Column number - use these to obtain the correct data
col_protein_time    = 1;
col_protein         = 2;
col_salt_time       = 3;
col_salt            = 4;

% 'Sec_UV'  'Protein mg/mL' 'Sec_pH'  'pH'
protein_time    = InputFile(:, col_protein_time); 
protein         = InputFile(:, col_protein); 
pH_time         = InputFile(:, col_salt_time); 
pH              = InputFile(:, col_salt);

% Set parameters and get parameter space
param_1 = initial(1);
param_2 = initial(2);
param_3 = initial(3);


%% FEED INPUTS

% Protein Loads
% We multiply a correction factor (UV295/UV280) to the protein loading 
% concentrations which are obtained from the appropriate CSV file
ProteinLoadCorrectionFactor = InputFile(1,9);
ProteinLoad         = InputFile(1,7);   % From nanodrop
ProteinLoad_295     = InputFile(1,8);   % From uv295 nm correlation
CADET_ProteinLoad   = ProteinLoad;

% Time Sections (remove superfluous NaNs)
TimeSections        = InputFile(:,6);
CADET_TimeSections  = TimeSections(~isnan(TimeSections));

% Flow Sections (remove superfluous NaNs)
FlowSections        = InputFile(:,5);
CADET_FlowSections  = FlowSections(~isnan(FlowSections));


% Protein data (protein_Exp)
    % For protein (Remember that experimental data is in mg/mL)
    % mM = mg/mL divided by Mol. Wt.
    protein_Exp = (1/mol_wt)*protein;
    
    % Apply correction factor for breakthrough and wash regions (for pH>8)
    protein_time_length = length(protein_time);
    [protein_Exp_delta, protein_Exp_index]= min(abs(protein_time-CADET_TimeSections(4)));
    protein_Exp(1:protein_Exp_index) = protein_Exp(1:protein_Exp_index)/ProteinLoadCorrectionFactor;
    
% pH data - Remove superfluous NaNs
    pH = pH(~isnan(pH));
    pH_time = pH_time(~isnan(pH_time));
    
run_name = list(i).name;
%% Solve the model

parameters = [param_1, param_2, param_3];
tic
    
% Solve the model
protein_model = PredictColloidalmodelEL_Rev_00(...
                CADET_ProteinLoad, ...      From nanodrop
                CADET_TimeSections, ...     
                CADET_FlowSections, ...
                pH_time,            ...
                pH,                 ...
                parameters);
toc

%% OUTPUT

% Get sum of squares
sum_of_sq = list(i).name;%rssq(protein_model - protein_Exp);

    
% Plot the solved model
figure %('units','normalized','outerposition',[0 0 1 1])
plot(...
    protein_time, protein_Exp, '.',        ...  Plot the experimental data
    protein_model(:,1), protein_model(:,3),'.', ...  Plot model solution
    'LineWidth', 3 ...
    )
yyaxis right
plot(...
    protein_model(:,1), -log10(protein_model(:,2)), 'LineWidth', 2)
    
    
    % Figure Title
    %title(list(i).name);
    
    % Figure Name
    figname = run_name; 
    figdir = 'C:\Users\Soumi\Documents\MATLAB\matlab\soumitra\MAbModel_NewCol_MAR21\';
    figname = strcat(figname(1:end-3), 'png');
    saveas(gcf,fullfile(figdir, figname), 'png')
    %close all


end

