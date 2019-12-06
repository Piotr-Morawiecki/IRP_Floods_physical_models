function [profiles] = plotRegionsForHighRainfall(z, t, hData, epsilon)
%timeMax = find(hData(end - 1, :) < 0, 1, 'last');

timeMax = size(hData, 2) - 1;
profiles = cell(timeMax);
h0 = -1 - z;
region2Boundaries = zeros(timeMax + 1, 2);
region2Boundaries(1, :) = [0; 0];

figure;
hold on
for time = 1:timeMax
    bottom = find(abs(hData(:, time) - h0) < epsilon, 1, 'first');
    if isempty(bottom)
        bottom = length(z);
    end
    top = find(hData(:, time) < 0, 1, 'first');
    profiles(time) = {hData(top:bottom, time)};
    if mod(time, 334) == 0
        plot(hData(top:bottom, time) ./ (-1-z(top:bottom)));
    end
    region2Boundaries(time + 1, :) = [z(bottom); z(top)];
end
legend(string(round(t(335:334:timeMax), 2)));
hold off

figure
subplot(3,1,1);
plot(t(1:(timeMax+1)), region2Boundaries(:, :));
xlabel('Time');
ylabel('Depth');
legend('Top boundary', 'Bottom boundary');
axis([0 t(timeMax+1) -Inf Inf]);

subplot(3,1,2);
plot(t(1:(timeMax+1)), diff(region2Boundaries, 1, 2));
xlabel('Time');
ylabel('Region II size');
axis([0 t(timeMax+1) -Inf Inf]);

subplot(3,1,3);
dt = t(2) - t(1);
averagingPeriod = 100;
speed = (region2Boundaries(1:(end - averagingPeriod), :) - ...
         region2Boundaries((averagingPeriod + 1):end, :)) / ...
         (averagingPeriod * dt);
plot(t(1:(timeMax - averagingPeriod + 1)), speed)
xlabel('Time');
ylabel('Boundary propagation speed');
axis([0 t(timeMax+1) -Inf Inf]);
legend('Top boundary', 'Bottom boundary');

end

