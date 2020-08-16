
% Convert gaze coodinates and stop radius to screen coordinates
screenCoord = @(gazeCoord) (-gazeCoord).*(1080/20.8692) + [1920 1080]/2;
symbolRadius = @(relRad) (relRad - min(mean(vqCleaned(:,3,:)))) * 30;

% Set up the video in and out
v = VideoReader('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_materials/StimulusFiles/PixarShorts.mov');
vo = VideoWriter('/Users/aguirre/Desktop/GazeTrack.mov');

% Set up the timebase. Account for a half-second phase shift that appears
% to be present between the eye tracking and the movie
movieStartTime = 1880;
phaseCorrect = -0.25;
timebaseSecs = (gazeData.timebase./1000) + movieStartTime + phaseCorrect;

nTrail = 0;

% Set up the symbol colors
nSubs = size(vqCleaned,1);
colors = getDistinguishableColors(nSubs).*255;

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
        f = insertShape(f,'circle',[thisCoord thisRadius],'Color',colors(ss,:));
        for rr = 1:min([nTrail tt-1])
            lineCoords = [screenCoord(squeeze(vqCleaned(ss,1:2,tt-rr+1))) ...
                screenCoord(squeeze(vqCleaned(ss,1:2,tt-rr)))];
            f = insertShape(f,'line',lineCoords,'LineWidth',3,'Color',colors(ss,:));            
        end
    end
    writeVideo(vo,f)
end

% Close and clear the video objects
close(vo);
clear v vo