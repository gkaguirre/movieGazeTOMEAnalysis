

% These are the fields to process
fieldNames = {'tfMRI_MOVIE_AP_run01','tfMRI_MOVIE_AP_run02','tfMRI_MOVIE_PA_run03','tfMRI_MOVIE_PA_run04'};

% The movie start times (in seconds) for each of the acquisitions
movieStartTimes = [1880, 2216, 892, 1228];

% Account for a quarter-second phase shift that appears to be present
% between the eye tracking and the movie
phaseCorrect = -0.25;

% How long a trail (in frames) do we leave behind each tracking circle?
nTrail = 0;

% Convert gaze coodinates and stop radius to screen coordinates
screenCoord = @(gazeCoord) (-gazeCoord).*(1080/20.8692) + [1920 1080]/2;
symbolRadius = @(relRad) (relRad - min(mean(vqCleaned(:,3,:)))) * 30;

% Set up the video in
v = VideoReader('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_materials/StimulusFiles/PixarShorts.mov');

% Loop over the fieldNames
for ff = 1:length(fieldNames)
    
    % Get this cleaned matrix
    vqCleaned = gazeData.(fieldNames{ff}).vqCleaned;
        
    % Set up the timebase.
    timebaseSecs = (gazeData.timebase./1000) + movieStartTimes(ff) + phaseCorrect;
    
    % Set up the symbol colors
    nSubs = size(vqCleaned,1);
    colors = getDistinguishableColors(nSubs).*255;
    
    % Set up the video out
    fileNameOut = sprintf('/Users/aguirre/Desktop/GazeTrack_%02d.avi',ff);
    vo = VideoWriter(fileNameOut);
    
    % Set the frame rate for the output
    vo.FrameRate = 1/(timebaseSecs(2)-timebaseSecs(1));
    
    % Open the video out object
    open(vo);
    
    % Loop through the frames
    for tt = 1:length(timebaseSecs)
        v.CurrentTime=timebaseSecs(tt);
        f = readFrame(v);
        for ss = 1:nSubs
            thisCoord = screenCoord(squeeze(vqCleaned(ss,1:2,tt)));
            thisRadius = symbolRadius(vqCleaned(ss,3,tt));
            if ~any(isnan([thisCoord thisRadius]))
                f = insertShape(f,'circle',[thisCoord thisRadius],'LineWidth',2,'Color',colors(ss,:));
                for rr = 1:min([nTrail tt-1])
                    lineCoords = [screenCoord(squeeze(vqCleaned(ss,1:2,tt-rr+1))) ...
                        screenCoord(squeeze(vqCleaned(ss,1:2,tt-rr)))];
                    f = insertShape(f,'line',lineCoords,'LineWidth',3,'Color',colors(ss,:));
                end
            end
        end
        writeVideo(vo,f)
    end
    
    % Close and clear the video objects
    close(vo);
    clear vo
    
end

clear v
