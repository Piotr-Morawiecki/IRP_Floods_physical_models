%% Make some data
%% Initialize video
myVideo = VideoWriter('LowRainfall'); %open video file
myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
open(myVideo)
%% Plot in a loop and grab frames
for i=3:3:996
    plot(z, hData(:, i), 'LineWidth', 1)
    ylim([-1, 0.2])
    xlim([-1, 0])
    xlabel('Depth, z')
    ylabel('Pressure head, h')
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

%plot(hData(:,1000))


plot(z(:), hData(:, 1500))
xlabel('Depth, z')
ylabel('Pressure head, h')
xlim([-0.46 -0.44])