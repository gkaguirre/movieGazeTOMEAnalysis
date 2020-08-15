
screenCoord = @(gazeCoord) (-gazeCoord).*(1080/20.8692) + [1920 1080]/2;
symbolRadius = @(relRad) (relRad - min(mean(vqCleaned(:,3,:)))) * 30;

v = VideoReader('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/TOME_materials/StimulusFiles/PixarShorts.mov');
vo = VideoWriter('/Users/aguirre/Desktop/GazeTrack.mov');

% Set up the timebase and the symbol colors
timebase = (gazeData.timebase./1000)+1880;
nSubs = size(vqCleaned,1);
colors = getDistinguishableColors(nSubs).*255;

% Set the frame rate for the output
vo.FrameRate = 1/(timebase(2)-timebase(1));

% Open the video out object
open(vo);

% Loop through the frames
for tt = 1:length(timebase)
    v.CurrentTime=timebase(tt);
    f = readFrame(v);
    for ss = 1:nSubs
        thisCoord = screenCoord(squeeze(vqCleaned(ss,1:2,tt)));
        thisRadius = symbolRadius(vqCleaned(ss,3,tt));
        f = insertShape(f,'circle',[thisCoord thisRadius],'Color',colors(ss,:));
    end
    writeVideo(vo,f)
end

close(vo);
clear v vo