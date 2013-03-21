function im = createSynthetic(type, width, height)
% createSynthetic(type, width, height)
% Creates a synthetic RGB image which can be mosaicked and reprocessed
%
% type - one of 'xsweep', 'ysweep', 'zoneplate' or 'radial', producing a
%   frequency sweep in the X or Y directions, a  or a radial pattern.
% width - Width of the resulting image
% height - Height of the resulting image
%
% The output image is full-color RGB, but all 3 channels are identical.

[x, y] = meshgrid(1:width, 1:height);

% Constant spatial frequency
% im = 0.5 * sin(x + y / 3) + 0.5;

if(strcmp(type, 'radial'))
  im = sin(10 * atan2(x-(width/2), y-(height/2)));
elseif(strcmp(type, 'xsweep'))
  im = sin(x.^2 / width);
elseif(strcmp(type, 'ysweep'))
 im = sin(y.^2 / width);
elseif(strcmp(type, 'zoneplate'))
  im = sin(((x - width/2).^2 + (y - height/2).^2) / width * pi / 2);
else
  error('Unknown type');
end

im = 0.5 * im + 0.5; % Convert from -1:1 to 0:1 range
im = repmat(im, [1 1 3]); % Make full color by replicating channels
