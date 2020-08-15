
% Load the gaze data file from my local disk
load('/Users/aguirre/Documents/MATLAB/projects/movieGazeTOMEAnalysis/data/gazeData.mat')

% Extract the data for the first acquisition
vq = gazeData.tfMRI_MOVIE_AP_run01.vq;

% How many subjects do we have
nSubs = size(vq,1);

% How many measures do we have?
nMeasures = size(vq,2);

% Store a record of the location of the nan values
vqNaN = isnan(vq);

% This is the mean across all subjects / time points for each measure
vqMeanVal = squeeze(nanmean(nanmean(vq),3));

% These are the mean centered vectors
vqCentered = vq - nanmean(vq,3);

% This is the mean vector across all subjects for each measure
vqMeanVec = squeeze(nanmean(vqCentered));

% Find the slope that relates each individual subject to the mean vector
for mm = 1:nMeasures
    for ii = 1:nSubs
        goodIdx = ~isnan(vqCentered(ii,mm,:));
        p = polyfit(squeeze(vqCentered(ii,mm,goodIdx)),vqMeanVec(mm,goodIdx)', 1);
        slopes(mm,ii) = p(1);
    end
    % Adjust the mean vector to remove compression
    vqMeanVec(mm,:) = vqMeanVec(mm,:) ./ mean(slopes(mm,:));
end

% Now loop through and adjust each subject's data

vqCenteredScaled = nan(size(vqCentered));
% Find the slope that relates each individual subject to the mean vector
for mm = 1:nMeasures
    for ii = 1:nSubs
        goodIdx = ~isnan(vqCentered(ii,mm,:));
        p = polyfit(squeeze(vqCentered(ii,mm,goodIdx)),vqMeanVec(mm,goodIdx)', 1);
        slopes(mm,ii) = p(1);
        vqCenteredScaled(ii,mm,goodIdx) = vqCentered(ii,mm,goodIdx).*p(1);
    end
    vqCenteredScaledSD(mm,:) = squeeze(nanstd(vqCenteredScaled(:,mm,:)))';
end

% Create a cleaned vq matrix that has "imputed" missing values with the
% mean across subject value, and adds back in the mean
vqCleaned = vqCenteredScaled;
for mm = 1:nMeasures
    for ii = 1:nSubs
        badIdx = isnan(vqCenteredScaled(ii,mm,:));
        vqCleaned(ii,mm,badIdx) = vqMeanVec(mm,badIdx);
        
        % add back in the mean to the whole vector
        vqCleaned(ii,mm,:) = vqCleaned(ii,mm,:) + vqMeanVal(mm);
    end
end

% Create a display of the data
figure
titles={'x gaze','y gaze','stop radius','nans'};
for mm = 1:size(vq,2)
    subplot(nMeasures+1,1,mm)
    imagesc(squeeze(vqCleaned(:,mm,:)));
    axis off
    title(titles{mm});
end
subplot(nMeasures+1,1,mm+1)
imagesc(squeeze(vqNaN(:,1,:)));
axis off
title(titles{mm+1});

