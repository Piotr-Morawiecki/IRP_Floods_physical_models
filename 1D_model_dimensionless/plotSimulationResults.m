function [] = plotSimulationResults(t, z, hData, qData)

    load('DefaultModel.mat', 'mod');

    thetaData = computeTheta(hData, mod);
    surfaceVolume = max(0, hData(1, :));
    subsurfaceVolume = mean(thetaData, 1);
    thetaData = effectiveSaturation(thetaData, mod);
    figure;

    subplot(4,1,1);
    image(t, z, hData,'CDataMapping','scaled');
    set(gca,'YDir','normal')
    title('Hydraulic head');
    colorbar;

    subplot(4,1,2);
    image(t, z, thetaData,'CDataMapping','scaled');
    set(gca,'YDir','normal')
    title('Effective saturation');
    colorbar;

    subplot(4,1,3);
    image(t, [0;z], qData,'CDataMapping','scaled');
    set(gca,'YDir','normal')
    title('Flow');
    colorbar;

    subplot(4,1,4);

    yyaxis left
    plot(t, subsurfaceVolume)
    yyaxis right
    plot(t, surfaceVolume)
    title('water volume per unit area');
    legend('subsuface volume', 'surface volume');
    title('Total water volume');
    
end

