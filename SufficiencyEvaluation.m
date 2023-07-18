% Henry Burton and Jack Baker
% July 23, 2023

% This code was developed by Henry Burton (UCLA) and Jack Baker (Stanford) 
% to perform and sufficiency-based evaluation of ground motion intensity 
% measures. It is meant to accompany the 
% following paper:

% Burton, H. V. and Baker, J. W. (2023) Evaluating the Effectiveness of 
% Ground Motion Intensity Measures Through the Lens of Causal Inference. 
% Earthquake Engineering Structural Dynamics (accepted for publication)

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           SPECIFY DIRECTORIES                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define path to base directory
BaseDirectory = mfilename('fullpath');
[BaseDirectory,~,~] = fileparts(BaseDirectory);

% Define path to directory with processed (.mat files) data
InputDataDirectoryProcessed = strcat(BaseDirectory,'\Data\',...
    'Processed\matFiles');

% Define path to "sufficiency" results directory
SufficiencyResultsDirectory = strcat(BaseDirectory,'\Results\Sufficiency');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             LOAD DATA                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go to input data directory
cd(InputDataDirectoryProcessed)

% Load EDPs
load allPSDR
load allPFA

%  Load IMs
load allPGA
load allPGV
load allSaavg
load allSaT1

% Load number of stories in each building
load NStories

% Load the upstream parameters. The first, second and third columns are the 
% magnitude, distance and epsilon, respectively.
load upstreamParameters

% Take natural log of distance
upstreamParameters(:,2) = log(upstreamParameters(:,2));

% Define the number of upstream parameters
Nupstream = length(upstreamParameters(1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      DEFINE MISCELLANEOUS VARIABLES                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of buildings
Nbldgs = 5;

% Number of ground motions
Ngms = 240;

% Number of intensity measures
Nims = 4;

% Define the subset of ground motions to run (one for each orthogonal
% direction)
gmsToRun = (1:1:240)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PERFORM BASELINE SUFFICIENCY EVALUATION               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform sufficiency evaluation considering all ground motions
[psdrPValues,pfaPValues,maxPSDRPValues,maxPFAPValues] = ...
    evaluateSufficiency(allPSDR,allPFA,allPGA,allPGV,allSaavg,allSaT1,...
    Nbldgs,NStories,Nims,Nupstream,upstreamParameters,Ngms,gmsToRun);

% Go to directory with .mat input data
cd(SufficiencyResultsDirectory)
save('psdrPValuesBaseline.mat','psdrPValues','-mat')
save('pfaPValuesBaseline.mat','pfaPValues','-mat')
save('maxPSDRPValuesBaseline.mat','maxPSDRPValues','-mat')
save('maxPFAPValuesBaseline.mat','maxPFAPValues','-mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          DEFINE FUNCTION USED TO EVALUATE SUFFICIENCY FOR ALL           %
%                      IMs, BUILDINGS AND EDPs                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function performs sufficiency evaluation for all EDPs, buildings and
% IMs
function [psdrPValues,pfaPValues,maxPSDRPValues,maxPFAPValues] = ...
    evaluateSufficiency(allPSDR,allPFA,allPGA,allPGV,allSaavg,allSaT1,...
    Nbldgs,NStories,Nims,Nupstream,upstreamParameters,Ngms,gmsToRun)

    % INPUTS
    % allPSDR: 1 X Nbldgs cell containing the PSDRs for all buildings, 
    % stories and ground motions  
    
    % allPFA: 1 X Nbldgs cell containing the PFAs for all buildings, 
    % stories and ground motions 
    
    % allPGA: Ngms x 1 array containing the peak ground acceleration (g)
    % for all ground motions
    
    % allPGV: Ngms x 1 array containing the peak ground velocity (cm/s)
    % for all ground motions
    
    % allSaavg: Ngms x Nbldgs array containing the Saavg for all buildings 
    % and for all ground motions
    
    % allSaT1: Ngms x Nbldgs array containing the SaT1 for all buildings 
    % and for all ground motions

    % Nbldgs - number of buildings (scalar)

    % NStories - NStories x 1 vector of the number of stories in each 
    % building

    % Nims - Number of IMs being considered

    % Nupstream - Number of upstream parameters to be considered

    % upstreamParameters - Ngms x Nupstream array of the upstream parameters
    % associated with each ground motion

    % Ngms - number of ground motions (scalar)

    % OUTPUTS
    % psdrPValues: 1 X Nbldgs cell containing the p-values obtained from 
    % the PSDR sufficiency assessment   

    % pfaPValues: 1 X Nbldgs cell containing the p-values obtained from 
    % the PFA sufficiency assessment  

    % Define cell array used to store sufficiency results
    psdrPValues = [];
    pfaPValues = [];
    maxPSDRPValues = zeros(Nbldgs,Nims,Nupstream);
    maxPFAPValues = zeros(Nbldgs,Nims,Nupstream);
    
    % Loop over the number of buildings
    for i = 1:Nbldgs
        % Extract EDPs for current building
        bldgPSDR = allPSDR{i};
        bldgPFA = allPFA{i};
    
        % Extract the number of stories for the currentr building
        nStories = NStories(i,1);

        % Extract maximum EDP over all stories
        for m = 1:Ngms
            maxPSDR(m,1) = (max((bldgPSDR(:,m))))'; 
            maxPFA(m,1) = (max((bldgPFA(:,m))))';
        end
        %maxPSDR = maxPSDR(gmsToRun,:);
        %maxPFA = maxPFA(gmsToRun,:);

        % Loop over the number of upstream parameters
        for j = 1:Nupstream
            % Intialize array used to store p-values
            psdrPValuesCurrent = zeros(nStories,Nims);
            pfaPValuesCurrent = zeros(nStories,Nims);

            % Compute p-value based on maxEDP for pga
            maxPSDRPValues(i,1,j) = computePValue(allPGA(gmsToRun,:),...
                maxPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
            maxPFAPValues(i,1,j) = computePValue(allPGA(gmsToRun,:),...
                maxPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
    
            % Compute p-value based on maxEDP for pgv
            maxPSDRPValues(i,2,j) = computePValue(allPGV(gmsToRun,:),...
                maxPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
            maxPFAPValues(i,2,j) = computePValue(allPGV(gmsToRun,:),...
                maxPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
    
            % Compute p-value based on maxEDP for SaT1
            maxPSDRPValues(i,3,j) = computePValue(allSaT1(gmsToRun,i),...
                maxPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
            maxPFAPValues(i,3,j) = computePValue(allSaT1(gmsToRun,i),...
                maxPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
    
            % Compute p-value based on maxEDP for Saavg
            maxPSDRPValues(i,4,j) = computePValue(allSaavg(gmsToRun,i),...
                maxPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
            maxPFAPValues(i,4,j) = computePValue(allSaavg(gmsToRun,i),...
                maxPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
    
            % Loop over the number of stories
            for k = 1:nStories
                % Extract EDPs for current story
                storyPSDR = (bldgPSDR(k,:))';
                flrPFA = (bldgPFA(k + 1,:))';
        
                % Compute p-value for pga
                psdrPValuesCurrent(k,1) = ...
                    computePValue(allPGA(gmsToRun,:),...
                    storyPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
                pfaPValuesCurrent(k,1) = ...
                    computePValue(allPGA(gmsToRun,:),flrPFA(gmsToRun,:),...
                    upstreamParameters(gmsToRun,j));
        
                % Compute p-value for pgv
                psdrPValuesCurrent(k,2) = ...
                    computePValue(allPGV(gmsToRun,:),...
                    storyPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
                pfaPValuesCurrent(k,2) = ...
                    computePValue(allPGV(gmsToRun,:),...
                    flrPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
        
                % Compute p-value for SaT1
                psdrPValuesCurrent(k,3) = ...
                    computePValue(allSaT1(gmsToRun,i),...
                    storyPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
                pfaPValuesCurrent(k,3) = ...
                    computePValue(allSaT1(gmsToRun,i),...
                    flrPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
        
                % Compute p-value for Saavg
                psdrPValuesCurrent(k,4) = ...
                    computePValue(allSaavg(gmsToRun,i),...
                    storyPSDR(gmsToRun,:),upstreamParameters(gmsToRun,j));
                pfaPValuesCurrent(k,4) = ...
                    computePValue(allSaavg(gmsToRun,i),...
                    flrPFA(gmsToRun,:),upstreamParameters(gmsToRun,j));
            end
            % Add p-values to cell array for current building and upstream
            % parameter
            psdrPValues{i,j} = psdrPValuesCurrent;
            pfaPValues{i,j} = pfaPValuesCurrent;
        end     
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FUNCTION USED TO COMPUTE P-VALUES FOR SUFFICIENCY EVALUATION       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes the p-value used in the sufficiency evaluation
function [pValue] = computePValue(IM,EDP,CausalParameter)

% INPUTS
% IM: Ngms x 1 array with the intensity measure associated with each ground
% motion

% EDP: Ngms x 1 array with the EDP obtained from the analysis with each 
% ground motion

% upstreamParameters - Ngms x 1 array of the upstream parameter being 
% considered

    % Define natural log of IMs
    logIM = log(IM);

    % Define natural log of EDP
    logEDP = log(EDP);
    
    % Fit linear model
    mdl = fitlm(logIM,logEDP);

    % Extract residuals
    residuals = mdl.Residuals.Raw;

    % Fit linear model
    mdl = fitlm(CausalParameter,residuals);
    pValue = coefTest(mdl);
end


%#ok<*AGROW> 