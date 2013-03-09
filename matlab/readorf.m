% Steven Bell
% 21 February 2013

% Use dcraw to convert the image from Olympus' ORF (or anything else)
% into a TIFF file.
% Don't convert colorspace (-o 0), don't demosiac (-D), save as 16-bit (-6)
% tiff file (-T)
% dcraw -o 0 -D -6 -T <image name>

imagePath = '../testimages/';
% im = imread([imagePath, 'PC150620.tiff']);
im = imread([imagePath, 'PC150634.tiff']);

bitDepth = 16; % Surely not really?
bitScale = 2^bitDepth;

% Cut the image down to a manageable size for now
offset = [0, 0]; % X, Y - Must be a multiple of bayer size (2)
dim = [2000 2000]; % X, Y - Must be a multiple of bayer size (2)
im = im(offset(2) + (1:dim(2)), offset(1) + (1:dim(1)));

%% Split the image by color channel

% Specify the bayer matrix.
% This assumes a zero (or multiple of 2) offset
bayer = ['r', 'g'; 'g', 'b'];

isRed   = repmat(bayer == 'r', size(im) / 2);
isGreen = repmat(bayer == 'g', size(im) / 2);
isBlue  = repmat(bayer == 'b', size(im) / 2);

red = zeros(size(im));
green = zeros(size(im));
blue = zeros(size(im));

red(logical(isRed)) = im(logical(isRed));
green(logical(isGreen)) = im(logical(isGreen));
blue(logical(isBlue)) = im(logical(isBlue));

% Normalization (maybe this isn't the right thing to do?)
red = red / max(red(:));
green = green / max(green(:));
blue = blue / max(blue(:));

% imshow(stack / bitScale);

% imshow(blue / bitScale);

%% Demosiacking
% Do the dumbest possible thing - nearest neighbor
redFilt = [1 1; 1 1];
greenFilt = [0.5, 0.5; 0.5, 0.5];
blueFilt = redFilt;

redDm = imfilter(red, redFilt);
greenDm = imfilter(green, greenFilt);
blueDm = imfilter(blue, blueFilt);

% imshow(greenDm);


%% Color correction matrix


%
stack = cat(3, redDm, greenDm, blueDm);
imshow(stack);


%% Save the result
imwrite(stack, 'out.tif', 'TIFF');

