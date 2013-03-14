
lrWidth = size(rawList{1}, 2);
lrHeight = size(rawList{1}, 1);

% Resolution of the "high resolution" output image
gridScale = 1; % Scale mapping new 
hrWidth = 320;
hrHeight = 240;


% Known offsets for each of the images in rawList
offsetX = [0, 3];
offsetY = [0, 3];

% Indices for the Bayer pattern in the low resolution images
% Assumes an R G / G B pattern
[redGridX, redGridY] = meshgrid(1:2:lrWidth, 1:2:lrHeight);
[greenGridX, greenGridY] = meshgrid(1:2:320, 1:1:240);
greenGridX(1:2:lrHeight, :) = greenGridX(1:2:lrHeight, :) + 1;
[blueGridX, blueGridY] = meshgrid(2:2:lrWidth, 2:2:lrHeight);

% Create logical masks for simple demosaicking
isRed = false(lrHeight, lrWidth);
isGreen = isRed;
isBlue = isRed;

isRed  (sub2ind([lrHeight, lrWidth], redGridY(:),   redGridX(:)))   = 1;
isGreen(sub2ind([lrHeight, lrWidth], greenGridY(:), greenGridX(:))) = 1;
isBlue (sub2ind([lrHeight, lrWidth], blueGridY(:),  blueGridX(:)))  = 1;

%%
redSamples = []; redPosX = []; redPosY = [];
greenSamples = []; greenPosX = []; greenPosY = [];
blueSamples = [];  bluePosX = []; bluePosY = [];

% For each color channel, create a vector with all of the samples, and
% corresponding vectors with the positions of the samples
for ii = 1:length(rawList)  
  im = rawList{ii};
  redSamples = cat(1, redSamples, im(isRed));
  redPosX = cat(1, redPosX, redGridX(:) - offsetX(ii));
  redPosY = cat(1, redPosY, redGridY(:) - offsetY(ii));

  % Have to use sub2ind here to keep the orders aligned
  greenSamples = cat(1, greenSamples, ...
    im(sub2ind([lrHeight, lrWidth], greenGridY(:), greenGridX(:))) );
  greenPosX = cat(1, greenPosX, greenGridX(:) - offsetX(ii));
  greenPosY = cat(1, greenPosY, greenGridY(:) - offsetY(ii));
  
  blueSamples = cat(1, blueSamples, im(isBlue));
  bluePosX = cat(1, bluePosX, blueGridX(:) - offsetX(ii));
  bluePosY = cat(1, bluePosY, blueGridY(:) - offsetY(ii));
end

% Draw a pretty plot showing where all the samples lie
% figure(1); clf; hold on;
% plot(redPosX, redPosY, 'rx');
% plot(greenPosX, greenPosY, 'gx');
% plot(bluePosX, bluePosY, 'bx');

%%
upGrid = zeros([hrHeight, hrWidth, 3]);
posX = {redPosX, greenPosX, bluePosX};
posY = {redPosY, greenPosY, bluePosY};
samples = {redSamples, greenSamples, blueSamples};

% Take the nearest sample
for ii = 1:hrHeight
  for jj = 1:hrWidth
    interpPointX = jj / gridScale; % Convert to LR coordinate space
    interpPointY = ii / gridScale;
    for c = 1:3 % Channel
      % Gaussian weighting
      % This really should be much sharper than Gaussian.  We want very close
      % points to be weighted extremely highly, but if two are roughly the
      % same distance but far away, they should be treated equally.
      weight = exp(-((posY{c} - interpPointY).^2 + (posX{c} - interpPointX).^2));

      normval = sum(weight);
      upGrid(ii, jj, c) = samples{c}' * weight / normval;
    end
  end
  fprintf('Finished row %d of %d\n', ii, size(upGrid, 1));
end

%% Nearest-neighbor demosiacking
redGrid = isRed .* rawList{1};
greenGrid = isGreen .* rawList{1};
blueGrid = isBlue .* rawList{1};

redFilt = [1 1; 1 1];
greenFilt = [0.5, 0.5; 0.5, 0.5];
blueFilt = redFilt;

redDm = imfilter(redGrid, redFilt);
greenDm = imfilter(greenGrid, greenFilt);
blueDm = imfilter(blueGrid, blueFilt);

stack = cat(3, redDm, greenDm, blueDm);
figure(1); imshow(stack);

%% Linear interpolation demosaicking
redGrid = isRed .* rawList{1};
greenGrid = isGreen .* rawList{1};
blueGrid = isBlue .* rawList{1};

redDm = imfilter(redGrid, [0.5, 1, 0.5]); % Interpolate horizontally
redDm = imfilter(redDm, [0.5; 1; 0.5]); % Interpolate vertically
greenDm = imfilter(greenGrid, [0, 0.25, 0; 0.25, 1, 0.25; 0, 0.25, 0]); % Or do it all at once
blueDm = imfilter(blueGrid, [0.25, 0.5, 0.25; 0.5, 1, 0.5; 0.25, 0.5, 0.25]);

stack = cat(3, redDm, greenDm, blueDm);
figure(2); imshow(stack);

%% MATLAB's built-in demosaicking
stack = demosaic(uint8(rawList{1} * 256), 'rggb');
figure(3); imshow(stack);

%% TODO: ISET demosaicking?

%%
figure(4); imshow(upGrid);
