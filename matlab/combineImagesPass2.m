

%% Get the second derivative of the red and blue channels
red = upGrid(:, :, 1);
green = upGrid(:, :, 2);
blue = upGrid(:, :, 3);

bddy = diff(blue, 2, 1); % Blue second derivative in the vertical direction
bddx = diff(blue, 2, 2); % Horizontal direction

rddy = diff(red, 2, 1); % Red second derivatives
rddx = diff(red, 2, 2);

ddx = (abs(rddx) + abs(bddx));
ddx = [zeros(size(ddx, 1), 1), ddx, zeros(size(ddx, 1), 1)]; % Pad back to full size

ddy = (abs(rddy) + abs(bddy));
ddy = [zeros(1, size(ddy, 2)); ddy; zeros(1, size(ddy, 2))]; % Pad back to full size

% figure(1); imagesc(ddx);
% figure(2); imagesc(ddy);

k = 0.1; % Regularization factor
alpha = (ddx + k) ./ (ddx + ddy + 2*k);
imagesc(alpha);

%%
upGrid2 = zeros([hrHeight, hrWidth, 3]);
posX = {redPosX, greenPosX, bluePosX};
posY = {redPosY, greenPosY, bluePosY};
samples = {redSamples, greenSamples, blueSamples};

% Take the nearest sample
for ii = 1:hrHeight
  for jj = 1:hrWidth
    interpPointX = jj / gridScale; % Convert to LR coordinate space
    interpPointY = ii / gridScale;

    % Gaussian weighting
    % If alpha is high, then the X derivative is high, meaning we should
    % trust the Y values more, meaning a higher variance gaussian.
    gWeight = exp(-((greenPosY - interpPointY).^2 * (1-alpha(ii, jj)) + ...
                    (greenPosX - interpPointX).^2 * alpha(ii, jj)) * 4);

    normval = sum(gWeight);
    upGrid2(ii, jj, 2) = greenSamples' * gWeight / normval;
  end
  fprintf('Finished row %d of %d\n', ii, size(upGrid, 1));
end


%%
upGrid2(:, :, 1) = red; % Add original red and blue channels
upGrid2(:, :, 3) = blue;

figure(1); imshow(upGrid2);

%%
oldUpGrid = zeros([hrHeight, hrWidth, 3]);
oldUpGrid(:, :, 2) = green;
figure(2); imshow(oldUpGrid);
