function result = getscielab(img1,img2,bit_depth, graphic_on)

%Distance to Display that the image will be displayed
vDist =0.3;

%% Explicit calls.  

% we scale by 2^bit_depth.
% Then we correct for gamma curve.
% The result is a 'linear' RGB between 0 and 1.
scale = 2^bit_depth;
img1_rgb  = dac2rgb(double(img1)/scale);
img2_rgb = dac2rgb(double(img2)/scale);


%% 2.  Load the display calibration information
% coneFile = fullfile(isetRootPath,'data','human','SmithPokornyCones');
% spCones = ieReadSpectra(coneFile,wave);   %plot(wave,spCones)

displayFile = fullfile(isetRootPath,'data','displays','crt');
dsp = displayCreate(displayFile);

% Determine transform matrix and white point.  
% Works with:  imageLinearTransform
rgb2xyz  = displayGet(dsp,'rgb2xyz');       % rowXYZ = rowRGB * rgb2xyz
whiteXYZ = displayGet(dsp,'white point');  

%% Convert the linear RGB data to XYZ values

img1XYZ = imageLinearTransform(img1_rgb,rgb2xyz);
img2XYZ = imageLinearTransform(img2_rgb,rgb2xyz);


%% Run the scielab function.

% The FOV depends on the display dpi, viewing distance, and image size
dots     = size(img1);
imgWidth = dots(2)*displayGet(dsp,'meters per dot');  % Image width (meters)
fov      = ieRad2deg(2*atan2(imgWidth/2,vDist));      % Horizontal fov in deg
sampPerDeg = dots(2)/fov;

% Run CIELAB 2000
params.deltaEversion = '2000';
params.sampPerDeg  = sampPerDeg;
params.imageFormat = 'xyz';
params.filterSize  = sampPerDeg;
params.filters = [];

[errorImage, params] = scielab(img1XYZ, img2XYZ, whiteXYZ, params);
result =mean(errorImage(:));

%% Examining and interpreting the results.
if(graphic_on ==1)
vcNewGraphWin;
imagesc(errorImage);
colorbar('vert');
title('S-CIELAB error map')

vcNewGraphWin;
hist(errorImage(:),100)
title('S-CIELAB delta E histogram')
figure;
imshow(img1);
figure;
imshow(img2);
end
end
