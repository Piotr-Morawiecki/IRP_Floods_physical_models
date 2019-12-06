z0Values = 1:100;
hProfiles = zeros(1700, size(z0Values, 2));
tStars = zeros(1, size(z0Values, 2));

for z0Id = 1:size(z0Values, 2)
    disp(z0Id)
    tic
    [hProfile, tStar] = simulateBoundaryGrowthGivenZ0(z0Values(z0Id));
    toc
    hProfiles(:, z0Id) = hProfile;
    tStars(1, z0Id) = tStar;
end

%alphaValues = alphaValues(1:313);
%hProfiles = hProfiles(:, 1:313);

save('hProfilesRefinedVaryingZ.mat', 'hProfiles')
%save('z.mat', 'z')