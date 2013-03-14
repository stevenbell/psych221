% Generates test images with various frequencies and patterns


width = 128;
height = 128;

[x, y] = meshgrid(1:width, 1:height);

% Constant spatial frequency
% im = 0.5 * sin(x + y / 3) + 0.5;

% Frequency sweep
% im = 0.5 * sin(x.^2 / width / 3) + 0.5;

% Radial pattern
im = 0.5 * sin(10 * atan2(x-(width/2), y-(height/2))) + 0.5;

imshow(im);

