function getPSNR(img1,img2,bit_depth)

err = (double(img1) - double(img2)).^2;
s = sum(err(:));

s
mse = double(s)/numel(img1)
max_i = 2^bit_depth
PSNR = 20*log10(max_i) - 10*log10(mse)

end
