% Henry Burton and Jack Baker
% July 23, 2023

% This code was developed by Henry Burton (UCLA) and Jack Baker (Stanford). 
% It takes the raw data files in .csv format and converts them to .mat files 
% so that they can be used by the "EfficiencyEvaluation.m" and 
% "SufficiencyEvaluation.m" scripts. It is meant to accompany the 
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


% Define path to directory with raw data (.csv files)
InputDataDirectoryRaw = strcat(BaseDirectory,'\Data\Raw');

% Define path to directory with processed (.mat files) data
InputDataDirectoryProcessed = strcat(BaseDirectory,'\Data\',...
    'Processed\matFiles');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             LOAD DATA                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go to raw input data directory
cd(InputDataDirectoryRaw)

% Factor to convert pga from cm/s^2 to g
pgaFactor = 0.00101971621;

% Ground motion information
gmSpectra = load("groundMotionSpectra.csv");
gmParameters = load("groundMotionParameters.csv");

% Building information
buildingInfo = load("BuildingInfo.csv");

% EDPs - PSDR
allPSDR{1} = load("1StoryPSDR.csv");
allPSDR{2}  = load("5StoryPSDR.csv");
allPSDR{3}  = load("9StoryPSDR.csv");
allPSDR{4}  = load("14StoryPSDR.csv");
allPSDR{5}  = load("19StoryPSDR.csv");

% EDPs - PFA
allPFA{1} = load("1StoryPFA.csv")*pgaFactor;
allPFA{2} = load("5StoryPFA.csv")*pgaFactor;
allPFA{3} = load("9StoryPFA.csv")*pgaFactor;
allPFA{4} = load("14StoryPFA.csv")*pgaFactor;
allPFA{5} = load("19StoryPFA.csv")*pgaFactor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DEFINE GLOBAL VARIABLES                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of ground motions
Ngms = 240;

% Number of buildings 
Nbldgs = 5;

% Number of IMs 
Nims = 4;

% Number of stories in each building
NStories = buildingInfo(:,1);

% Number of upstream parameters
Nupstream = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       EXTRACT CAUSAL PARAMETERS                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intialize upstream parameters array
upstreamParameters = zeros(Ngms,Nupstream);

% Append upstream parameters array
upstreamParameters(:,1) = gmParameters(:,1);
upstreamParameters(:,2) = gmParameters(:,3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       EXTRACT PGA,PGV, SaT1 and Saavg                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract pga and pgv
allPGA = gmParameters(:,4)*pgaFactor;
allPGV = gmParameters(:,5);

% Extract SaT1 and Saavg
TStartRatio = 0.2;
TEndRatio = 1.5;
allT1 = buildingInfo(:,2);
[allSaT1,allSaavg] = extractSaT1Saavg(allT1,TStartRatio,TEndRatio,...
    Ngms,Nbldgs,gmSpectra);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       EXTRACT maxPSDR and maxPFA                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract maximum EDP over all stories
for i = 1:Nbldgs
    % Extract EDPs for current building
    bldgPSDR = allPSDR{i};
    bldgPFA = allPFA{i};    
    for m = 1:Ngms
        maxPSDR(m,i) = (max((bldgPSDR(:,m))))'; 
        maxPFA(m,i) = (max((bldgPFA(:,m))))';
    end
end

% Go to directory with processed input data and store relevant arrays
cd(InputDataDirectoryProcessed)
save('allPGA.mat','allPGA','-mat')
save('allPGV.mat','allPGV','-mat')
save('allSaT1.mat','allSaT1','-mat')
save('allSaavg.mat','allSaavg','-mat')
save('allSaT1.mat','allSaT1','-mat')
save('allSaavg.mat','allSaavg','-mat')
save('allPSDR.mat','allPSDR','-mat')
save('allPFA.mat','allPFA','-mat')
save('upstreamParameters.mat','upstreamParameters','-mat')
save('NStories.mat','NStories','-mat')
save('maxPSDR.mat','maxPSDR','-mat')
save('maxPFA.mat','maxPFA','-mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  FUNCTION USED TO EXTRACT SaT1 and Saavg                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function extracts the first mode spectral acceleration (SaT1) 
% and geomean  of spectral acceleration over a range of periods (Saavg)
function [allSaT1,allSaavg] = extractSaT1Saavg(allT1,TStartRatio,...
    TEndRatio,Ngms,Nbldgs,gmSpectra)

% INPUTS
% allT1 - Nbldgs x 1 vector of 1st mode period for all buildings

% TStartRatio - lower bound period used to compute Saavg defined as a
% fraction of T1  (scalar)

% TEndRatio - upper bound period used to compute Saavg defined as a
% fraction of T1  (scalar)

% Ngms - number of ground motions (scalar)

% Nbldgs - number of buildings (scalar)

% gmSpectra - 800 x Ngms + 1 array of ground motion spectra

% OUTPUTS
% allSaT1 - Ngms x Nbldgs array of first mode spectral acceleration for
% all ground motions and the periods corresponding to all buildings

% allSaavg - Ngms x Nbldgs array of geomean  of spectral acceleration over
% % a range of periods for all ground motions and the periods corresponding 
% % to all buildings

    % Initialize array used to store allSaT1 and allSaavg
    allSaT1 = zeros(Ngms,Nbldgs);
    allSaavg = zeros(Ngms,Nbldgs);

    
    % Loop over the number of buildings
    for i = 1:Nbldgs
        % Compute SaT1 for building i and ground motion j
        T1 = allT1(i,1);
        [~,idx] = min(abs(gmSpectra(:,1) - T1));
        SaT1 = gmSpectra(idx,(2:Ngms + 1));
        allSaT1(:,i) = SaT1';
    
        % Compute Saavg parameters for building i 
        TStart = TStartRatio*T1;
        TEnd = TEndRatio*T1;
        [~,idxStart] = min(abs(gmSpectra(:,1) - TStart));
        [~,idxEnd] = min(abs(gmSpectra(:,1) - TEnd));
        for j = 1:Ngms
            % Compute Saavg for building i and ground motion j
            allSaavg(j,i) = geomean(gmSpectra(idxStart:idxEnd,j + 1));
        end
    end
end

%#ok<*SAGROW> 