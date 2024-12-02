%% Input

mainFolder          = 'A';
saveFolderResults   = 'B';
paperFigFolder      = 'C';

load('filteredMSASKatsu.mat')
katsuMAINSHOCKS = selectedMAINSHOCKS; katsuMSspectra = selectedMSspectra;
katsuAFTERSHOCKS = selectedAFTERSHOCKS; katsuASspectra = selectedASspectra;
clear selectedMAINSHOCKS selectedMSspectra selectedAFTERSHOCKS selectedASspectra

load('simbad150.mat', 'GMs', 'SPECTRA')
simbadGMs = GMs; simbadSPECTRA = SPECTRA;
clear GMs SPECTRA

nGMfromSimbad   = 100;
SFforSimbad     = 1;
SFforKatsu      = 1;

periodsForAvgSA = linspace(0.18, 1.5*1.19, 10); % For Aftershock Eh study

nSamplesLHS = 1000;
nSamplesRND = nSamplesLHS;

FreeVibrationMS = 20;
FreeVibrationAS = 0; % not needed for ruaumoko

limitForSignificantDuration = 0.9;

%% Get the highest GMs from simbad

IMsimbad = zeros(size(simbadGMs,1),1);
for gm = 1 : size(simbadGMs,1)
    simbadGMs{gm}{1}(:,2)       = SFforSimbad * simbadGMs{gm}{1}(:,2);
    simbadSPECTRA{gm}{1}(:,2)   = SFforSimbad * simbadSPECTRA{gm}{1}(:,2);
    
    IMsimbad(gm) = geomean( interp1(simbadSPECTRA{gm}{1}(:,1), simbadSPECTRA{gm}{1}(:,2), periodsForAvgSA) );
    
    cumulativeSquaredAcceleration = cumsum( simbadGMs{gm}{1}(:,2).^2 );
    firstStep = find(cumulativeSquaredAcceleration >= 0.01*cumulativeSquaredAcceleration(end) );
    lastStep = find(cumulativeSquaredAcceleration >= limitForSignificantDuration*cumulativeSquaredAcceleration(end) );
    significantDurationSIMBAD(gm,1) = interp1( ...
        cumulativeSquaredAcceleration(firstStep(1):lastStep(1)), ...
        simbadGMs{gm}{1}(firstStep(1):lastStep(1),1), ...
        limitForSignificantDuration*cumulativeSquaredAcceleration(end) );
end

[ ~, indicesSimbad ]    = sort(IMsimbad, 'descend');

selectedSimbadGMs       = simbadGMs( indicesSimbad(1:nGMfromSimbad) );
selectedSimbadSPECTRA   = simbadSPECTRA( indicesSimbad(1:nGMfromSimbad) );
selectedSimbadIM        = IMsimbad( indicesSimbad(1:nGMfromSimbad) );
selectedSimbadSD        = significantDurationSIMBAD( indicesSimbad(1:nGMfromSimbad) );

clear simbadGMs simbadSPECTRA IMsimbad

%% Calculate intensity measure for the Katsu goda sequences (and scale them)

selectedMAINSHOCKS = katsuMAINSHOCKS;
selectedAFTERSHOCKS = katsuAFTERSHOCKS;
selectedMSspectra = katsuMSspectra;
selectedASspectra = katsuASspectra;
IMmsas = zeros(size(selectedMAINSHOCKS,1) , 2);
for gm = 1 : size(selectedMAINSHOCKS,1)
    selectedMAINSHOCKS{gm}{1}(:,2)  = SFforKatsu * katsuMAINSHOCKS{gm}{1}(:,2);
    selectedMSspectra{gm}{1}(:,2)   = SFforKatsu * katsuMSspectra{gm}{1}(:,2);
    selectedAFTERSHOCKS{gm}{1}(:,2) = SFforKatsu * katsuAFTERSHOCKS{gm}{1}(:,2);
    selectedASspectra{gm}{1}(:,2)   = SFforKatsu * katsuASspectra{gm}{1}(:,2);
    
    IMmsas(gm,:) = [  geomean( interp1(selectedMSspectra{gm}{1}(:,1), ...
        selectedMSspectra{gm}{1}(:,2), periodsForAvgSA) ) ...
        geomean( interp1(selectedASspectra{gm}{1}(:,1), ...
        selectedASspectra{gm}{1}(:,2), periodsForAvgSA ) )  ];
    
    cumulativeSquaredAcceleration = cumsum( selectedMAINSHOCKS{gm}{1}(:,2).^2 );
    firstStep = find(cumulativeSquaredAcceleration >= 0.04*cumulativeSquaredAcceleration(end) );
    lastStep = find(cumulativeSquaredAcceleration >= limitForSignificantDuration*cumulativeSquaredAcceleration(end) );
    significantDurationMS(gm,1) = interp1( ...
        cumulativeSquaredAcceleration(firstStep(1):lastStep(1)), ...
        selectedMAINSHOCKS{gm}{1}(firstStep(1):lastStep(1),1), ...
        limitForSignificantDuration*cumulativeSquaredAcceleration(end) );
    
end
significantDurationMS(:,2) = significantDurationMS(:,1);
clear katsuMAINSHOCKS katsuAFTERSHOCKS katsuMSspectra katsuASspectra

% add simbad IMs at the end
IMmsas(size(selectedMAINSHOCKS,1)+1 : size(selectedMAINSHOCKS,1)+nGMfromSimbad, :) = [ selectedSimbadIM selectedSimbadIM ];
significantDurationMS(size(selectedMAINSHOCKS,1)+1 : size(selectedMAINSHOCKS,1)+nGMfromSimbad, :) = [ selectedSimbadSD selectedSimbadSD ];


combinedMAINSHOCKS = selectedMAINSHOCKS;
combinedMAINSHOCKS(size(selectedMAINSHOCKS,1)+1 : size(selectedMAINSHOCKS,1)+nGMfromSimbad) = selectedSimbadGMs;

combinedMSspectra = selectedMSspectra;
combinedMSspectra(size(selectedMAINSHOCKS,1)+1 : size(selectedMAINSHOCKS,1)+nGMfromSimbad)  = selectedSimbadSPECTRA;

%% Compose new sequences

% latin hypercube sampling
lhsSamples      = ceil(lhsdesign(nSamplesLHS, 2) * size(IMmsas,1));
lhsFromSimbad   = any(lhsSamples > size(selectedMAINSHOCKS,1), 2);

% back to back combinations
[p,q]           = meshgrid(1:size(IMmsas,1), 1:size(IMmsas,1));
B2Bsamples      = [p(:) q(:)];
B2BFromSimbad   = any(B2Bsamples > size(selectedMAINSHOCKS,1), 2);

h1 = histogram2(B2Bsamples(:,1), B2Bsamples(:,2), 'Normalization', 'pdf', 'FaceAlpha', 0.5);
hold on
h2 = histogram2(lhsSamples(:,1), lhsSamples(:,2), h1.XBinEdges, h1.YBinEdges, 'Normalization', 'pdf', 'FaceAlpha', 0.5);

%% Check optimal LHS sampling

sampleSize = [100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 20000 30000 40000];

error = nan(size(sampleSize));
for s = 1 : numel(sampleSize)
    
    lhsSamp      = ceil(lhsdesign(sampleSize(s), 2) * size(IMmsas,1));
    
    h1 = histogram2(B2Bsamples(:,1), B2Bsamples(:,2), 'Normalization', 'pdf', 'Visible', 'off');
    hold on
    h2 = histogram2(lhsSamp(:,1), lhsSamp(:,2), h1.XBinEdges, h1.YBinEdges, 'Normalization', 'pdf', 'Visible', 'off');
    
    residual    = reshape(h2.Values, numel(h2.Values), 1) - reshape(h1.Values, numel(h1.Values), 1);
    error(s)    = sqrt(sum(residual .* residual));
end

close
figure; hold on
plot(sampleSize, error, '-o')
scatter(nSamplesLHS, interp1(sampleSize, error, nSamplesLHS), 100, 'filled')
legend('Trend', 'Chosen')
xlabel('LHS sample size')
ylabel('Misfit from Back2Back')
set(gca, 'FontSize', 16)

%% Plot

close all

maxIMplot   = max(IMmsas(:,2));
colour1     = 0.7 * [1 1 1];
colour2     = 0.0 * [1 1 1];
sizeMarker  = 15;
font        = 18;

% scatter plots

figure('Position', [1    34   560   420]); hold on
scatter(IMmsas(:,1), IMmsas(:,2),           2*sizeMarker, colour1, 'filled')
scatter(selectedSimbadIM, selectedSimbadIM, 2*sizeMarker, colour2, 'filled')
xlabel(sprintf('avgSA_{MS}(%1.1fs-%1.1fs)', periodsForAvgSA(1), periodsForAvgSA(end)))
ylabel(sprintf('avgSA_{AS}(%1.1fs-%1.1fs)', periodsForAvgSA(1), periodsForAvgSA(end)))
title(sprintf('Natural MS-AS (%d seqs.)',size(selectedMAINSHOCKS,1)))
set(gca,'FontSize',font)
axis([0 maxIMplot 0 maxIMplot])

figure('Position', [558    34   560   420]); hold on
scatter(IMmsas(B2Bsamples(B2BFromSimbad==0,1),1), IMmsas(B2Bsamples(B2BFromSimbad==0,2),1), sizeMarker/3, colour1, 'filled')
scatter(IMmsas(B2Bsamples(B2BFromSimbad==1,1),1), IMmsas(B2Bsamples(B2BFromSimbad==1,2),1), sizeMarker/3, colour1, 'filled')
xlabel(sprintf('avgSA_{MS}(%1.1fs-%1.1fs)', periodsForAvgSA(1), periodsForAvgSA(end)))
ylabel(sprintf('avgSA_{AS}(%1.1fs-%1.1fs)', periodsForAvgSA(1), periodsForAvgSA(end)))
%title(sprintf('Back2Back MS-MS (%d seqs.)', size(B2Bsamples,1)))
set(gca,'FontSize',font)
axis([0 maxIMplot 0 maxIMplot])

annotation(gcf,'textbox', [0 0.883333333333338 0.0579285714285718 0.0738095238095277], 'String','a)',...
'LineStyle','none', 'FontSize',18, 'FontName','Helvetica Neue', 'FitBoxToText','off');

figure('Position', [1121    34   560   420]); hold on
scatter(IMmsas(lhsSamples(lhsFromSimbad==0,1),1), IMmsas(lhsSamples(lhsFromSimbad==0,2),1), 2*sizeMarker, colour1, 'filled', 'MarkerEdgeColor', 'k')
scatter(IMmsas(lhsSamples(lhsFromSimbad==1,1),1), IMmsas(lhsSamples(lhsFromSimbad==1,2),1), 2*sizeMarker, colour1, 'filled', 'MarkerEdgeColor', 'k')
xlabel(sprintf('avgSA_{MS}(%1.1fs-%1.1fs)', periodsForAvgSA(1), periodsForAvgSA(end)))
ylabel(sprintf('avgSA_{AS}(%1.1fs-%1.1fs)', periodsForAvgSA(1), periodsForAvgSA(end)))
%legend('From original database', 'MS&/AS from SIMBAD')
%title(sprintf('LHS MS-MS (%d seqs.)', nSamplesLHS))
set(gca,'FontSize',font)
axis([0 maxIMplot 0 maxIMplot])

ax = axes('Parent',gcf,...
    'Position',[0.571428571428571 0.569047619047619 0.326785714285714 0.347619047619049]);
scatter(ax, significantDurationMS(lhsSamples(lhsFromSimbad==0,1),1), significantDurationMS(lhsSamples(lhsFromSimbad==0,2),1), 5, colour1, 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(ax, significantDurationMS(lhsSamples(lhsFromSimbad==1,1),1), significantDurationMS(lhsSamples(lhsFromSimbad==1,2),1), 5, colour1, 'filled', 'MarkerEdgeColor', 'k')
xlabel('SD5-95_{MS}'); ylabel('SD5-95_{AS}')
%title(sprintf('Back2Back MS-MS (%d seqs.)', size(B2Bsamples,1)))
set(gca,'FontSize',font)

annotation(gcf,'textbox', [0 0.883333333333338 0.0579285714285718 0.0738095238095277], 'String','b)',...
'LineStyle','none', 'FontSize',18, 'FontName','Helvetica Neue', 'FitBoxToText','off');

%% Create MSAS sequences

forAnalysisMAINSHOCKS                                   = combinedMAINSHOCKS(lhsSamples(lhsFromSimbad==0,1),:); % only from Katsu
forAnalysisMAINSHOCKS(end+1:end+sum(lhsFromSimbad==1))  = combinedMAINSHOCKS(lhsSamples(lhsFromSimbad==1,1),:); % MS and/or AS from SIMBAD

forAnalysisAFTERSHOCKS                                  = combinedMAINSHOCKS(lhsSamples(lhsFromSimbad==0,2),:);
forAnalysisAFTERSHOCKS(end+1:end+sum(lhsFromSimbad==1)) = combinedMAINSHOCKS(lhsSamples(lhsFromSimbad==1,2),:);

forAnalysisMSspectra                                    = combinedMSspectra(lhsSamples(lhsFromSimbad==0,1),:);
forAnalysisMSspectra(end+1:end+sum(lhsFromSimbad==1))   = combinedMSspectra(lhsSamples(lhsFromSimbad==1,1),:);

forAnalysisASspectra                                    = combinedMSspectra(lhsSamples(lhsFromSimbad==0,2),:);
forAnalysisASspectra(end+1:end+sum(lhsFromSimbad==1))   = combinedMSspectra(lhsSamples(lhsFromSimbad==1,2),:);


RUAUresultsINTERVAL = 0.02; % [s]
[ forAnalysisMSAS,  importantSTEPS ] = sequenceMSASmaker(forAnalysisMAINSHOCKS, forAnalysisAFTERSHOCKS, RUAUresultsINTERVAL, FreeVibrationMS, FreeVibrationAS);

lastMSASfromKatsu = sum(lhsFromSimbad==0);

for gm = size(forAnalysisMSAS,1) : -1 : 1
    stepForResiduals(gm,1) = ( importantSTEPS{1}(gm,1) * RUAUresultsINTERVAL ) / forAnalysisMSAS{gm}{1}(2,1);
end

% control plot
figure; hold on
for gm = 1 : size(forAnalysisMSAS,1)
    plot(forAnalysisMSAS{gm}{1}(:,1), forAnalysisMSAS{gm}{1}(:,2))
end
xlabel('time [s]'); ylabel('Acceleration [g]')
set(gca,'FontSize', 16)

%% Find non-unique main-shock ground motions

duplicatedMS = zeros(numel(forAnalysisMAINSHOCKS), 1);
for GM = 1 : numel(forAnalysisMAINSHOCKS)
    for gm = GM+1 : numel(forAnalysisMAINSHOCKS)
        if numel(forAnalysisMAINSHOCKS{GM}{1}(:,2))==numel(forAnalysisMAINSHOCKS{gm}{1}(:,2)) && ...
                all(forAnalysisMAINSHOCKS{GM}{1}(:,2) == forAnalysisMAINSHOCKS{gm}{1}(:,2))
            
            % figure; hold on
            % plot(forAnalysisMAINSHOCKS{GM}{1}(:,1), forAnalysisMAINSHOCKS{GM}{1}(:,2))
            % plot(forAnalysisMAINSHOCKS{gm}{1}(:,1), forAnalysisMAINSHOCKS{gm}{1}(:,2))
            
            duplicatedMS(gm) = 1;
            
        end
    end
end

%% Save (the random seed is not set. Watch out before saving)

% save(fullfile(mainFolder, 'forAnalysisMSAS.mat'), 'forAnalysisMAINSHOCKS', ...
%     'forAnalysisAFTERSHOCKS', 'forAnalysisMSspectra', 'forAnalysisASspectra', ...
%     'forAnalysisMSAS', 'importantSTEPS', 'lastMSASfromKatsu', 'duplicatedMS')
% 
% if strcmp(analysisType, 'Karim')
%     save(fullfile(mainFolder, 'forAnalysisKARIM.mat'), 'forAnalysisMSspectra', ...
%         'forAnalysisASspectra', 'forAnalysisMSAS', 'stepForResiduals', ...
%         'lastMSASfromKatsu', 'duplicatedMS')
% end
