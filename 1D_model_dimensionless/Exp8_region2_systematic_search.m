alphaValues = 0.4:0.1:40;
hProfiles = zeros(1700, size(alphaValues, 2));
tStars = zeros(1, size(alphaValues, 2));

for alphaId = 1:size(alphaValues, 2)
    disp(alphaId)
    tic
    [hProfile, tStar] = simulateBoundaryGrowth(alphaValues(alphaId));
    toc
    hProfiles(:, alphaId) = hProfile;
    tStars(1, alphaId) = tStar;
end

%alphaValues = alphaValues(1:313);
%hProfiles = hProfiles(:, 1:313);

save('hProfilesRefined2.mat', 'hProfiles')
%save('z.mat', 'z')

interestingProfiles = [1 7 97];
plot(z, hProfiles(:, interestingProfiles))
axis([-0.04 0 -1 0.1])
xlabel('z')
ylabel('h')
legend('\alpha=0.4', '\alpha=1.0', '\alpha=2.0')

region2Width = zeros(size(alphaValues, 2), 1);
for alphaId = flip(1:size(alphaValues, 2))
    top = z(find(hProfiles(1:(end - 1), alphaId) > 0, 1, 'last'));
    bottom = z(hProfiles(:, alphaId) == min(hProfiles(:, alphaId)));
    region2Width(alphaId) = top - bottom;
end

plot(alphaValues, region2Width)

figure;
loglog(alphaValues, region2Width, '.')
hold on
loglog(alphaValues, 0.0122 * alphaValues .^ (-1), '-')
hold off
grid on
xlabel('Alpha')
ylabel('Region II width')
legend('Numerical results', 'Fitted function y = 0.0122 \alpha^{-1}')


loglog(alphaValues(1:161), tStars(1, 1:161), '.')
hold on
loglog(alphaValues, 0.0032 * alphaValues .^ (-0.872), '-')
hold off
grid on
xlabel('Alpha')
ylabel('Time until saturation at h=0')
legend('Numerical results', 'Fitted function y = 0.0032 \alpha^{-0.872}')
