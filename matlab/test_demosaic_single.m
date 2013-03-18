clear;
file = {'../testimages/fruits.jpg' , '../testimages/PC150634.JPG'}

%% GtoR, GtoB are optimal values in simulation results
GtoR = 0.097;
GtoB = 0.1;

for i = 1: size(file,2)
    
        [psnr,psnr_ref,scielab,scielab_ref] = demosaic_single(file{i},GtoR,GtoB,1)
            
end

    
