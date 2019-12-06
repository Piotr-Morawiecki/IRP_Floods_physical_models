totalTime = 4;
nTimeSteps = 10000;
dt = totalTime / (nTimeSteps - 1);
t = dt .* (0:(nTimeSteps - 1));

totalLength = 1;
lengthSteps = 31;
dx = totalLength / (lengthSteps - 1);
x = dx .* (0:(lengthSteps - 1));

hData = zeros(lengthSteps, nTimeSteps);
hData(:, 1) = 0;

storeData = hData;
theta = 0.009;
n = 1/2; % Chezy
%n = 2/3; % Manning

qData = zeros(lengthSteps, nTimeSteps - 1);

drainage = zeros(lengthSteps, nTimeSteps);
rainDuration = 2500;
rainIntensity = 1;
drainage(:, 1:rainDuration) = rainIntensity;

for time = 2:nTimeSteps
    h = hData(:, time - 1);
    q = - h .* diff([h.^n; 0]) ./ dx;
    h = h + ( -diff([0; q]) ./ dx + drainage(:, time)) .* dt;
    hData(:, time) = h;
    qData(:, time - 1) = q;
    
    storeData(:, time) = drainage(:, time) .* dt + ...
                         (1 - theta) * storeData(:, time - 1) + ...
                         theta * [0; storeData(1:(end - 1), time - 1)];
end

figure;
subplot(3,2,1);
image(t, x, hData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('Groundwater table height (kinematic wave)');
colorbar;

subplot(3,2,2);
image(t, x, storeData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('Groundwater table height (Box model)');
colorbar;

subplot(3,2,3);
image(t, x, qData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('Water flow (kinematic wave)');
colorbar;

subplot(3,2,4);
image(t, x, theta * storeData,'CDataMapping','scaled');
set(gca,'YDir','normal')
title('Water flow (Box model)');
colorbar;

subplot(3,2,5);
plot(t(1:(end-1)), qData(end, :))
title('River flow (kinematic wave)');

subplot(3,2,6);
plot(t, theta * storeData(end, :) * dx / dt)
title('River flow (Box model)');