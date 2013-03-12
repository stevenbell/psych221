function getPSNR(img1,img2,bit_depth)
% 
rows = size(img1,1);
cols = size(img1,2);
num_channel = size(img1,3);
count =0;
sum =0;
for ii = 1:rows
    for jj = 1: cols
        for cc = 1: num_channel
             sum = sum + (img1(ii,jj,cc)-img2(ii,jj,cc))^2;
             count = count +1;
    end
    end
end


sum
count
mse = double(sum)/count
max_i = 2^bit_depth
PSNR = 20*log10(max_i) - 10*log10(mse)

end
