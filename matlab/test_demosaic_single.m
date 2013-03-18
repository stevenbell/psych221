clear;
file = {'../testimages/fruits.jpg' , '../testimages/tiger.jpg' ,'../testimages/tiger2.jpg' ,'../testimages/fish.jpg' ,'../testimages/PC150634.JPG'}

%% GtoR, GtoB are optimal values in simulation results
GtoR = 0.097;
GtoB = 0.1;

for i = 1: size(file,2)
        [psnr(i),psnr_ref(i),scielab(i),scielab_ref(i)] = demosaic_single(file{i},GtoR,GtoB,1)
end

    
