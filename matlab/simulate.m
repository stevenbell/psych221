% Simulate the demosiacking/super-resolution process
clear;

%% Read in the source image
source = im2double(imread('../testimages/PC150634.JPG'));

% Downsample to reduce any artifacts from the original image and to have
% something of reasonable size to work with
sourceSmall = imresize(source, 0.25);

imshow(sourceSmall);

%% Create the artificial Bayer matrix
rawSize = [240, 320];

bayerPattern = [1, 2; 2, 3]; % 1=red, 2=green, 3=blue
bayerFull = false([rawSize, 3]);

for c = 1:3
  bayerFull(:, :, c) = repmat(bayerPattern == c, rawSize / 2);
end

%% Create a number of sample mosiacked images
downFactor = 2;
sliceWidth = rawSize(2) * downFactor;
sliceHeight = rawSize(1) * downFactor;
sliceX = 100;
sliceY = 100;
%%
% TODO: Seed the random number generator so we get repeatable results
nSamples = 2;
x = [0, 10, 0, 1];
y = [0, 9, 1, 1];

rawList = cell(nSamples, 1);

for ii = 1:nSamples
  % This offset is in terms of original-image pixels
  offsetX = x(ii); % Eventually we'll use a random offset
  offsetY = y(ii);

  % Convolve
  xFraction = rem(offsetX, 1);
  yFraction = rem(offsetX, 1);
  f = [xFraction; 1-xFraction] * [yFraction, 1-yFraction];
  subSampled = imfilter(sourceSmall, f);
  
  % Downsample
  subSampled = subSampled((1:downFactor:sliceHeight) + sliceY - fix(offsetY), ...
                          (1:downFactor:sliceWidth) + sliceX - fix(offsetX), :);

  % Mosiac
  raw = zeros([rawSize, 3]);
  raw(bayerFull) = subSampled(bayerFull);
  raw = sum(raw, 3); % Flatten the image

  rawList{ii} = raw;
end

%% Now demosaic these images!


% Calculate the relative alignment
a = rawList{1};
b = rawList{2};

aFreq = fft2(a);
bFreq = fft2(b);

% TODO: drop high frequencies introduced by the Bayer pattern
g = aFreq .* conj(bFreq);
pcorr = fftshift(abs(ifft2(g ./ abs(abs(g)))));

[yy, xx] = find(pcorr == max(pcorr(:)));
imagesc(pcorr);
offset = [yy, xx] - rawSize / 2 - [1 1]

% TODO: subpixel registeration

% imagesc(log(fftshift(abs(aFreq))));
% colormap(gray);

% Splat

% Sharpen?



