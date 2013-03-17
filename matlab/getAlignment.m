function [x, y] = getAlignment(imSet)
% Calculate the relative X and Y offsets needed to align a set of images.
% The returned offsets have subpixel resolution.

nImages = length(imSet);
width = size(imSet{1}, 2);
height = size(imSet{1}, 1);

% For each image, calculate the relative alignment to all other images.
% Since half of these are duplicates, only calculate vs the remaining ones.

alignX = zeros([nImages, nImages]);
alignY = zeros([nImages, nImages]);

% Precalculate all the FFTs
imFreq = cell([1, nImages]);
for ii = 1:nImages
  imFreq{ii} = fft2(imSet{ii});
end

for baseIdx = 1:nImages-1
  for targetIdx = baseIdx+1:nImages
    % Use phase correlation
    % TODO: drop high frequencies introduced by the Bayer pattern
    g = imFreq{baseIdx} .* conj(imFreq{targetIdx});
    pcorr = fftshift(abs(ifft2(g ./ abs(abs(g)))));
    [yy, xx] = find(pcorr == max(pcorr(:)));
    
    % Do subpixel registration by fitting a 2D parabola
    x = repmat(-3:3, [1, 7])';
    y = reshape(repmat(-3:3, [7, 1]), [49, 1]);

    v = [x.^2, y.^2, x.*y, x, y, ones(49, 1)];
    z = reshape(pcorr(yy-3:yy+3, xx-3:xx+3)', [49, 1]);

    c = v'*v \ v' * z; % Least squares fit to find the parameters

    % Solve for the peak
    xInterp = ((2*c(2)*c(4)/c(3)) - c(5)) / (c(3) - (4*c(2)*c(1)/c(3)));
    yInterp = -(2*c(1)*xInterp + c(4)) / c(3);

%     figure(3); clf;
%     imagesc(pcorr); hold on;
%     plot(xInterp, yInterp, 'wx');
%     pause();

    % Calculate the actual offset from the peak position
    alignX(baseIdx, targetIdx) = xInterp + xx - width / 2 - 1;
    alignY(baseIdx, targetIdx) = yInterp + yy - height / 2 - 1;
  end
end

%% Perform a least-squares fit to find the best offset guesses

% Each image pair gives us an equation (in each dimension) relating the
% absolute positions of the two images to the measurement of their offset.
% We also define the constraint that the first image has an offset of zero.
% There will be sum(1:(nImages-1)) + 1 equations (e.g, 4 images gives 3
% pairs + 2 pairs + 1 pair + 1 "grounding" constraint = 7.

nEq = sum(1:(nImages-1)) + 1;
A = zeros(nEq, nImages);
bx = zeros(nEq, 1);
by = zeros(nEq, 1);

A(1, 1) = 1; % Write the first equation
count = 2; % Row to put the next equation into

for baseIdx = 1:nImages-1
  for targetIdx = baseIdx+1:nImages
    A(count, baseIdx) = 1;
    A(count, targetIdx) = -1;
    bx(count) = alignX(baseIdx, targetIdx);
    by(count) = alignY(baseIdx, targetIdx);
    count = count + 1;
  end
end

x = (A'*A) \ A' * bx;
y = (A'*A) \ A' * by;

