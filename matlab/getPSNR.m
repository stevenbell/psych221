function psnr = getPSNR(img1,img2,bit_depth)

err = (double(img1) - double(img2)).^2;
s = sum(err(:));

mse = double(s)/numel(img1);
max_i = 2^bit_depth;
psnr = 20*log10(max_i) - 10*log10(mse);

end
