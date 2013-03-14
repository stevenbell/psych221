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
noiseVar = 0.001; % Signal has range 0-1

%%
% Seed the random number generator so we get repeatable results
RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', 0));

nSamples = 2;
x = [0, 6, .25, .125];
y = [0, 6, 0, 0];

rawList = cell(nSamples, 1);

for ii = 1:nSamples
  % This offset is in terms of original-image pixels
  offsetX = x(ii); % Eventually we'll use a random offset
  offsetY = y(ii);

  % Convolve
  xFraction = rem(offsetX, 1);
  yFraction = rem(offsetY, 1);
  f = [yFraction; 1-yFraction] * [xFraction, 1-xFraction];
  subSampled = imfilter(sourceSmall, f);
  
  % Downsample
  subSampled = subSampled((1:downFactor:sliceHeight) + sliceY - fix(offsetY), ...
                          (1:downFactor:sliceWidth) + sliceX - fix(offsetX), :);

  % Mosiac
  raw = zeros([rawSize, 3]);
  raw(bayerFull) = subSampled(bayerFull);
  raw = sum(raw, 3); % Flatten the image

  % Add noise
  raw = raw + randn(rawSize) * sqrt(noiseVar);
  
  rawList{ii} = raw;
end

%% Align the images

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

%% Do subpixel registration by fitting a 2D parabola
x = repmat([xx - 1; xx; xx + 1], [3, 1]);
y = reshape(repmat([yy - 1, yy, yy + 1], [3, 1]), [9, 1]);

v = [x.^2, y.^2, x.*y, x, y, ones(9, 1)];
z = reshape(pcorr(yy-1:yy+1, xx-1:xx+1)', [9, 1]);

c = inv(v'*v) * v' * z;

xInterp = ((2*c(2)*c(4)/c(3)) - c(5)) / (c(3) - (4*c(2)*c(1)/c(3)));
yInterp = -(2*c(1)*xInterp + c(4)) / c(3);

intOffset = [yy, xx] - rawSize / 2 - [1 1]
floatOffset = [yInterp, xInterp] - rawSize / 2 - [1 1]

hold on;
plot(xInterp, yInterp, 'wx');


%% Subpixel registration by upsampling

factor = 1;
a = imresize(rawList{1}, factor);
b = imresize(rawList{2}, factor);

% a(size(a, 1) * factor, size(a, 2) * factor) = 0; % Hack to pad with zeros
% b(size(b, 1) * factor, size(b, 2) * factor) = 0;

aFreq = fft2(a);
bFreq = fft2(b);

g = aFreq .* conj(bFreq);

% Drop high frequencies introduced by the Bayer pattern
middleY = size(aFreq, 1) / 2;
middleX = size(aFreq, 2) / 2;
band = 100;

aFreq(middleY-band:middleY+band, :) = 0;
aFreq(:, middleX-band:middleX+band) = 0;

% imagesc(log(abs(g)));
% colormap gray;

% imshow(abs(ifft2(aFreq)));


pcorr = fftshift(abs(ifft2(g ./ abs(abs(g)))));

[yy, xx] = find(pcorr == max(pcorr(:)));
imagesc(pcorr);

offset = [yy, xx] - rawSize * factor / 2 - [1 1]

% imagesc(log(fftshift(abs(aFreq))));
% colormap(gray);

% Splat

% Sharpen?


%% Perform demosiacking

