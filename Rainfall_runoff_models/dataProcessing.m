imageData = read(Tiff('\\myfiles\pwm27\data\KS_sl7.tif','r'));
regiona1Data = imageData(3450:3615,3760:3965);
%imshow(regiona1Data);
catchmentLocation = read(Tiff('catch1.tif','r'));
catchmentLocation = catchmentLocation(:,:,1);
catchmentLocation = imresize(catchmentLocation,size(regiona1Data));
%ukData = imageData(2500:3900,2500:4000);
%imshow(ukData);
%catchment1Data = imageData(3480:3650,3800:3850);
%catchment1Data = imageData(3300:3580,3790:4000);
imshow(catchmentLocation);
%imshow(imfuse(catchmentLocation,regiona1Data));
%imshow(catchmentLocation);