% Henry Burton and Jack Baker
% July 23, 2023

% This code was developed by Henry Burton (UCLA) and Jack Baker (Stanford) 
% to perform and efficiency-based evaluation of ground motion intensity 
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

% Define path to "efficiency" results directory
EfficiencyResultsDirectory = strcat(BaseDirectory,'\Results\Efficiency');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      DEFINE MISCELLANEOUS VARIABLES                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of buildings
Nbldgs = 5;

% Number of ground motions
Ngms = 240;

% Number of intensity measures
Nims = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PERFORM BASELINE EFFICIENCY EVALUATION               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform efficiency evaluation considering all ground motions
[psdrSigmas,pfaSigmas,psdrMeanSigmas,pfaMeanSigmas,maxPSDRSigmas,...
    maxPFASigmas] = evaluateEfficiency(allPSDR,allPFA,allPGA,allPGV,...
    allSaavg,allSaT1,Nbldgs,NStories,Nims,Ngms);

% Go to directory with .mat input data
cd(EfficiencyResultsDirectory)
save('psdrSigmasBaseline.mat','psdrSigmas','-mat')
save('pfaSigmasBaseline.mat','pfaSigmas','-mat')
save('maxPSDRSigmasBaseline.mat','maxPSDRSigmas','-mat')
save('maxPFASigmasBaseline.mat','maxPFASigmas','-mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          DEFINE FUNCTION USED TO EVALUATE EFFICIENCY FOR ALL            %
%                      IMs, BUILDINGS AND EDPs                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function performs efficiency evaluation for all EDPs, buildings and
% IMs
function [psdrSigmas,pfaSigmas,psdrMeanSigmas,pfaMeanSigmas,maxPSDRSigmas,...
    maxPFASigmas] = ...
    evaluateEfficiency(allPSDR,allPFA,allPGA,allPGV,allSaavg,allSaT1,...
    Nbldgs,NStories,Nims,Ngms)

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
    
    % Ngms - Number of grund motions

    % OUTPUTS
    % psdrSigmas: 1 X Nbldgs cell containing the PSDR beta values for all 
    % buildings and stories   
    
    % pfaSigmas: 1 X Nbldgs cell containing the PFA beta values for all 
    % buildings and stories

    % psdrMeanSigmas: Nbldgs x Nims array of mean beta values for PSDR
    % considering all buildings
    
    % pfaMeanSigmas: Nbldgs x Nims array of mean beta values for PFA
    % considering all buildings

    % Define cell array used to store efficiency results
    psdrSigmas = [];
    pfaSigmas = [];
    maxPSDRSigmas = zeros(Nbldgs,Nims);
    maxPFASigmas = zeros(Nbldgs,Nims);

    % Loop over the number of buildings
    for i = 1:Nbldgs
        % Extract EDPs for current building
        bldgPSDR = allPSDR{i};
        bldgPFA = allPFA{i};

        % Extract the number of stories for the currentr building
        nStories = NStories(i,1);
    
        % Initialize array to store sigmas for current building
        psdrSigmasCurrent = zeros(nStories,Nims);
        pfaSigmasCurrent = zeros(nStories,Nims);

        % Extract maximum EDP over all stories
        for m = 1:Ngms
            maxPSDR(m,1) = (max((bldgPSDR(:,m))))'; 
            maxPFA(m,1) = (max((bldgPFA(:,m))))';
        end

        % Compute theta based on maxEDP for pga
        maxPSDRSigmas(i,1) = computeSigma(allPGA,maxPSDR);
        maxPFASigmas(i,1) = computeSigma(allPGA,maxPFA);

        % Compute theta based on maxEDP for pgv
        maxPSDRSigmas(i,2) = computeSigma(allPGV,maxPSDR);
        maxPFASigmas(i,2) = computeSigma(allPGV,maxPFA);

        % Compute theta based on maxEDP for SaT1
        maxPSDRSigmas(i,3) = computeSigma(allSaT1(:,i),maxPSDR);
        maxPFASigmas(i,3) = computeSigma(allSaT1(:,i),maxPFA);

        % Compute theta based on maxEDP for Saavg
        maxPSDRSigmas(i,4) = computeSigma(allSaavg(:,i),maxPSDR);
        maxPFASigmas(i,4) = computeSigma(allSaavg(:,i),maxPFA);


        % Loop over the number of stories
        for j = 1:nStories
            % Extract EDPs for current story
            storyPSDR = (bldgPSDR(j,:))';
            flrPFA = (bldgPFA(j + 1,:))';
    
            % Compute efficiency for pga
            psdrSigmasCurrent(j,1) = computeSigma(allPGA,storyPSDR');
            pfaSigmasCurrent(j,1) = computeSigma(allPGA,flrPFA');
    
            % Compute efficiency for pgv
            psdrSigmasCurrent(j,2) = computeSigma(allPGV,storyPSDR');
            pfaSigmasCurrent(j,2) = computeSigma(allPGV,flrPFA');
    
            % Compute efficiency for SaT1
            psdrSigmasCurrent(j,3) = ...
                computeSigma(allSaT1(:,i),storyPSDR');
            pfaSigmasCurrent(j,3) = ...
                computeSigma(allSaT1(:,i),flrPFA');
    
            % Compute efficiency for Saavg
            psdrSigmasCurrent(j,4) = ...
                computeSigma(allSaavg(:,i),storyPSDR');
            pfaSigmasCurrent(j,4) = ...
                computeSigma(allSaavg(:,i),flrPFA');
        end
    
        % Add sigmas to cell array for current building
        psdrSigmas{i} = psdrSigmasCurrent;
        pfaSigmas{i} = pfaSigmasCurrent;

        % Loop over the number of IMs
        for k = 1:Nims
            % Compute and store mean sigma values
            psdrMeanSigmas(i,k) = mean(psdrSigmasCurrent(:,k));
            pfaMeanSigmas(i,k) = mean(pfaSigmasCurrent(:,k));
        end

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   FUNCTION USED TO COMPUTE IM SIGMA                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes the standard deviation of the residuals (beta)
% from regressing the IM versus the EDP in log space.
function [beta] = computeSigma(IM,EDP)

% INPUTS
% IM: Ngms x 1 array with the intensity measure associated with each ground
% motion
% EDP: Ngms x 1 array with the EDP obtained from the analysis with each 
% ground motion

    % Define natural log of IMs
    logIM = log(IM);

    % Define natural log of EDP
    logEDP = log(EDP);
    
    % Fit linear model
    mdl = fitlm(logIM,logEDP);
    beta = std(mdl.Residuals.Raw);
end


%#ok<*AGROW> 