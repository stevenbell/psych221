% simulate.m
% Top-level script for the demosiacking/super-resolution process
% Steven Bell, Donguk Moon, and Artem Vassiliev
% PSYCH 221, Spring 2013

clear;

%% Read in the source image

% Use a large JPEG image, but downsample to reduce any artifacts from the
% original image and to have something of reasonable size to work with
sourceFull = im2double(imread('../testimages/PC150634.JPG'));
source = imresize(sourceFull, 0.25);
[srcX, srcY, width, height] = deal(20, 20, 640, 480);

% Alternatively, create a synthetic image
% source = createSynthetic('zoneplate', 400, 300);
% [srcX, srcY, width, height] = deal(20, 20, 320, 240);

% Don't crop the source; the mosaicking method will do that
imshow(source);
hold on;
rectangle('Position', [srcX, srcY, width, height], 'EdgeColor', [1 1 0.5]);

%% Create artificial mosaicked images
% Seed the random number generator so we get repeatable results
RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', 0));

% trueOffsetX = [0, 1, 0, 1];
% trueOffsetY = [0, 0, 1, 1];
trueOffsetX = [0, rand(1, 3) * 5]; % Always use 0 for first offset
trueOffsetY = [0, rand(1, 3) * 5];

% Multiply offsets by 2 to work in input-image-space
rawList = createMosaicked(source, length(trueOffsetX), srcX, srcY, width, height,...
  trueOffsetX * 2, trueOffsetY * 2);

%% Align the images
[dx, dy] = getAlignment(rawList);

trueOffsetX = trueOffsetX - trueOffsetX(1);
trueOffsetY = trueOffsetY - trueOffsetY(1);

alignmentError = [trueOffsetX', trueOffsetY'] - [dx, dy]

% Use the ground truth offsets
offsetX = trueOffsetX;
offsetY = trueOffsetY;

% Use the calculated offsets
% offsetX = dx;
% offsetY = dy;

%% Perform demosiacking
% This expects to find rawList, offsetX, and offsetY in the workspace
combineImages

%%
% combineImagesPass2
