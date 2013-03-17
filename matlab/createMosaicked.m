function imSet = createMosaicked(imIn, nSamples, baseX, baseY, width, height, offsetX, offsetY)
% Create a set of random mosiacked images from a single RGB image
% imIn - Source RGB image to use, as a matrix of double
% nSamples - Number of mosiacked images to generate
% baseX, baseY - Base offset into the source image
% width, height - Dimensions of the input image region to use
% offsetX, offsetY - Offsets to produce, in pixels of the input image
%
% Returns a cell array with each of the mosaicked images as a 2D matrix.
% The resulting images will be half the size specified by the 'width' and
% 'height' parameters.
%
% The image sizes must be a multiple of 4, since the image is downsampled
% by a factor of 2, and the result must be a complete Bayer pattern.
% This code uses a Bayer 'rggb' pixel layout.

%% Error checking
if(numel(offsetX) ~= nSamples || numel(offsetY) ~= nSamples)
  error('The number of specified offsets must match the number of requested samples');
end
if(mod(width, 4) ~= 0 || mod(height, 4) ~= 0)
  error('width and height must be a multiple of 4.');
end

%% Create the artificial Bayer matrix
rawSize = [height, width] / 2; % Size of the output image

bayerPattern = [1, 2; 2, 3]; % 1=red, 2=green, 3=blue
bayerFull = false([rawSize, 3]);

for c = 1:3
  bayerFull(:, :, c) = repmat(bayerPattern == c, rawSize / 2);
end

%% Create a number of sample mosiacked images
noiseVar = 0.001; % Signal has range 0-1
imSet = cell(nSamples, 1);

for ii = 1:nSamples
  % Offset is specified in input-image-pixels
  dx = offsetX(ii);
  dy = offsetY(ii);

  % Convolve
  xFraction = rem(dx, 1);
  yFraction = rem(dy, 1);
  f = [yFraction; 1-yFraction] * [xFraction, 1-xFraction];
  subSampled = imfilter(imIn, f);
  
  % Downsample
  subSampled = subSampled((1:2:height) + baseY - fix(dy), ...
                          (1:2:width) + baseX - fix(dx), :);

  % Mosiac
  raw = zeros([rawSize, 3]);
  raw(bayerFull) = subSampled(bayerFull);
  raw = sum(raw, 3); % Flatten the image

  % Add noise
  raw = raw + randn(rawSize) * sqrt(noiseVar);
  
  imSet{ii} = raw;
end
