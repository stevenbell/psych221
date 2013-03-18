function [psnr, psnr_ref, scielab, scielab_ref] = demosaic_single (file, GtoR, GtoB, display_on)
%% Tunable Paramters
GtoRed =GtoR; % effectiveness of Green for Red Demosaicing 
GtoBlue = GtoB;% effectiveness of Green for Blue Demosaicing 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing Input Test Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read a test image
source = im2double(imread(file));
img_size = size(source);
img_size = img_size(1:2);

%% Make image size multiples of 2
img_size(1) = floor(img_size(1)/2) *2
img_size(2) = floor(img_size(2)/2) *2
source = source(1:img_size(1),1:img_size(2),:);

%% Create the artificial Bayer matrix
bayerPattern = [1, 2; 2, 3]; % 1=red, 2=green, 3=blue
bayerFull = false([img_size,3]);

%% Full Bayer Pattern generation
for c = 1:3
  bayerFull(:, :, c) = repmat(bayerPattern == c, img_size / 2);
end
  bayerFull_v = repmat(bayerPattern, img_size / 2);

%% Generating mosaiced RGB images from the test image
%% Apply Bayer Pattern to the source image
for c = 1:3
image_mosaic(:,:,c) = source(:,:,c) .* double(bayerFull(:,:,c));
end
image_mosaic_ref = image_mosaic;
image_mosaic_RGBcombined = image_mosaic(:,:,1) + image_mosaic(:,:,2) + image_mosaic (:,:,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actual Interpolation Kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Interpolation sequence
%% 1. interpolate Green channel with second-derivative of Red or Blue 
%% 2. interpolate Red and Blue channels with intercolor correlation factor  

%% Computation for Blue grandiednt matrix
%% Gbx2 = Bx-2 + Bx+2 -2*Bx
%% Gby2 = By-2 + By+2 -2*By
padd = [0,2];
Bx_n2 = padarray(image_mosaic (:,:,3),padd,0,'pre');
Bx_p2 = image_mosaic(:,3:(img_size(2)),3);
Bx_p2 = padarray(Bx_p2,[0,4],0,'post');
Bx = padarray(image_mosaic (:,:,3),padd,0,'post');
Gbx2 = Bx_n2 + Bx_p2 -2 *Bx;

padd = [2,0];
By_n2 = padarray(image_mosaic (:,:,3),padd,0,'pre');
By_p2 = image_mosaic(3:(img_size(1)),:,3);
By_p2 = padarray(By_p2,[4,0],0,'post');
By = padarray(image_mosaic (:,:,3),padd,0,'post');
Gby2 = By_n2 + By_p2 -2 *By;

%% Computation for Red grandiednt matrix
%% Grx2 = Rx-2 + Rx+2 -2*Rx
%% Gry2 = Ry-2 + Ry+2 -2*Ry
padd = [0,2];
Rx_n2 = padarray(image_mosaic (:,:,1),padd,0,'pre');
Rx_p2 = image_mosaic(:,3:(img_size(2)),1);
Rx_p2 = padarray(Rx_p2,[0,4],0,'post');
Rx = padarray(image_mosaic (:,:,1),padd,0,'post');
Grx2 = Rx_n2 + Rx_p2 -2 *Rx;

padd = [2,0];
Ry_n2 = padarray(image_mosaic (:,:,1),padd,0,'pre');
Ry_p2 = image_mosaic(3:(img_size(1)),:,1);
Ry_p2 = padarray(Ry_p2,[4,0],0,'post');
Ry = padarray(image_mosaic (:,:,1),padd,0,'post');
Gry2 = Ry_n2 + Ry_p2 -2 * Ry;


%% Filling in Green value
%% Using 2nd gradients in x,y direction to consider edges
for r = 3: (img_size(1)-3)
  for c = 3: (img_size(2)-3)
    
    if(bayerFull_v(r,c) == 1) % interpolation at Red pixel points  
     %% interpolation terms in row and column direction
     blur_r = image_mosaic(r-1,c,2) + image_mosaic(r+1,c,2);
     blur_c = image_mosaic(r,c-1,2) + image_mosaic(r,c+1,2);
     
     %% if both x and y grandients are zero, we give them the same weight 
     if (Grx2(r,c) == 0 && Gry2(r,c) == 0)
         r_weight = 0.5;
         c_weight = 0.5;
     else
       r_weight = abs(Grx2(r,c))/(abs(Grx2(r,c)) + abs(Gry2(r,c)));
       c_weight = abs(Gry2(r,c))/(abs(Grx2(r,c)) + abs(Gry2(r,c)));
     end 

     image_mosaic(r,c,2) = 0.5*r_weight*blur_r + 0.5*c_weight*blur_c;

    elseif(bayerFull_v(r,c) == 3) % interpolation at Blue pixel points  
     %% interpolation terms in row and column direction
     blur_r = image_mosaic(r-1,c,2) + image_mosaic(r+1,c,2);
     blur_c = image_mosaic(r,c-1,2) + image_mosaic(r,c+1,2);
     
     %% if both x and y grandients are zero, we give them the same weight 
     if (Gbx2(r,c) == 0 && Gby2(r,c) == 0)
         r_weight = 0.5;
         c_weight = 0.5;
     else
     r_weight = abs(Gbx2(r,c))/(abs(Gbx2(r,c)) + abs(Gby2(r,c)));
     c_weight = abs(Gby2(r,c))/(abs(Gbx2(r,c)) + abs(Gby2(r,c)));
      end
      
      image_mosaic(r,c,2) = 0.5*r_weight*blur_r + 0.5*c_weight*blur_c;

    end
  end
end

% Evaluate Green prediction error
g_err = image_mosaic(3:img_size(1)-3, 3: img_size(2)-3,2) - source(3:img_size(1)-3, 3: img_size(2)-3,2);
g_err = abs(g_err);
gerr_avg= sum(sum(g_err,1),2)/(size(g_err,1))/ (size(g_err,2));


%% 2. Interpolation for the Red and Blue channels 
%% Filling in R,B holes in two stages
%% First Pass  - fill the unknown values at pixel where adjecent values of the same color  are available    
%% Second Pass - fill the rest of unkown values by bilinear interpolation

for r = 2: (img_size(1)-1)
  for c = 2: (img_size(2)-1)
    
    %% Filling in Red values
    %% g_factor implies green value contribution
    if(bayerFull_v(r,c) == 2 && bayerFull_v(r,c-1) == 1) %RGR  
     if(image_mosaic(r,c-1,2) ~= 0 && image_mosaic(r,c+1,2) ~=0)
       g_factor = 0.5*(image_mosaic(r,c-1,1)/image_mosaic(r,c-1,2) + image_mosaic(r,c+1,1)/image_mosaic(r,c+1,2));
     elseif (image_mosaic(r,c-1,2) ~= 0 && image_mosaic(r,c+1,2) ==0)
       g_factor = (image_mosaic(r,c-1,1)/image_mosaic(r,c-1,2));
     elseif (image_mosaic(r,c-1,2) == 0 && image_mosaic(r,c+1,2) ~=0)
       g_factor = (image_mosaic(r,c+1,1)/image_mosaic(r,c+1,2));
     else 
       g_factor =1.0;
     end
     image_mosaic(r,c,1) = (1.0-GtoRed)*0.5*(image_mosaic(r,c-1,1) + image_mosaic(r,c+1,1)) + GtoRed*min(image_mosaic(r,c,2)*g_factor,1.0);
    elseif(bayerFull_v(r,c) == 2 && bayerFull_v(r-1,c) == 1) %R;G;R  
     if(image_mosaic(r-1,c,2) ~= 0 && image_mosaic(r+1,c,2) ~=0)
       g_factor = 0.5*(image_mosaic(r-1,c,1)/image_mosaic(r-1,c,2) + image_mosaic(r+1,c,1)/image_mosaic(r+1,c,2));
     elseif(image_mosaic(r-1,c,2) ~= 0 && image_mosaic(r+1,c,2) ==0)
       g_factor = (image_mosaic(r-1,c,1)/image_mosaic(r-1,c,2));
     elseif(image_mosaic(r-1,c,2) == 0 && image_mosaic(r+1,c,2) ~=0)
       g_factor = (image_mosaic(r+1,c,1)/image_mosaic(r+1,c,2));
     else 
       g_factor =1.0;
     end
     image_mosaic(r,c,1) = (1.0-GtoRed)*0.5*(image_mosaic(r-1,c,1) + image_mosaic(r+1,c,1)) + GtoRed*min(image_mosaic(r,c,2)*g_factor,1.0);
    end

    %% Filling in Blue values
    if(bayerFull_v(r,c) == 2 && bayerFull_v(r,c-1) == 3) %BGB  
     if(image_mosaic(r,c-1,2) ~= 0 && image_mosaic(r,c+1,2) ~=0)
       g_factor = 0.5*(image_mosaic(r,c-1,3)/image_mosaic(r,c-1,2) + image_mosaic(r,c+1,3)/image_mosaic(r,c+1,2));
     elseif(image_mosaic(r,c-1,2) ~= 0 && image_mosaic(r,c+1,2) ==0)
       g_factor = (image_mosaic(r,c-1,3)/image_mosaic(r,c-1,2));
     elseif(image_mosaic(r,c-1,2) == 0 && image_mosaic(r,c+1,2) ~=0)
       g_factor = (image_mosaic(r,c+1,3)/image_mosaic(r,c+1,2));
     else 
       g_factor =1.0;
     end
     image_mosaic(r,c,3) = (1.0-GtoBlue)*0.5*(image_mosaic(r,c-1,3) + image_mosaic(r,c+1,3)) + GtoBlue*min(image_mosaic(r,c,2)*g_factor,1.0);

    elseif(bayerFull_v(r,c) == 2 && bayerFull_v(r-1,c) == 3) %B;G;B  
     if(image_mosaic(r-1,c,2) ~= 0 && image_mosaic(r+1,c,2) ~=0)
       g_factor = 0.5*(image_mosaic(r-1,c,3)/image_mosaic(r-1,c,2) + image_mosaic(r+1,c,3)/image_mosaic(r+1,c,2));
     elseif(image_mosaic(r-1,c,2) ~= 0 && image_mosaic(r+1,c,2) ==0)
       g_factor = (image_mosaic(r-1,c,3)/image_mosaic(r-1,c,2));
     elseif(image_mosaic(r-1,c,2) == 0 && image_mosaic(r+1,c,2) ~=0)
       g_factor = (image_mosaic(r+1,c,3)/image_mosaic(r+1,c,2));
     else 
       g_factor =1.0;
     end
     image_mosaic(r,c,3) = (1.0-GtoRed)*0.5*(image_mosaic(r-1,c,3) + image_mosaic(r+1,c,3)) + GtoBlue*min(image_mosaic(r,c,2)*g_factor,1.0);
    end

  end % inner for loop end
end % outer for loop end


%% Second Pass - fill the rest of holes (bilinear interpolation)
for r = 2: (img_size(1)-1)
  for c = 2: (img_size(2)-1)
    % filling in red
    if(bayerFull_v(r,c) == 3)   % Blue spots should be filled with Red
      image_mosaic(r,c,1) = 0.25* ( image_mosaic(r-1,c,1) + image_mosaic(r+1,c,1) + image_mosaic(r,c-1,1)+ image_mosaic(r,c+1,1));
    end
    if(bayerFull_v(r,c) == 1)   % Red spots should be filled with Blue
      image_mosaic(r,c,3) = 0.25* ( image_mosaic(r-1,c,3) + image_mosaic(r+1,c,3) + image_mosaic(r,c-1,3)+ image_mosaic(r,c+1,3));
    end
  end
end


%% Excluding boundary region
source_cropped = source (7: (img_size(1)-6), 7: (img_size(2)-6), :);
image_cropped = image_mosaic(7: (img_size(1)-6), 7: (img_size(2)-6), :);
%image_mosaic_ref_cropped = image_mosaic_ref(7: (img_size(1)-6), 7: (img_size(2)-6), :);
image_demosaic_ref = demosaic(uint8(image_mosaic_RGBcombined *2^8),'rggb');
image_demosaic_ref_cropped = image_demosaic_ref(7: (img_size(1)-6), 7: (img_size(2)-6), :);

%% Performance Metric
psnr = getPSNR(source_cropped*2^8,image_cropped*2^8,8);
psnr_ref = getPSNR(source_cropped*2^8,image_demosaic_ref_cropped*2^8,8);
scielab = getscielab(source_cropped*2^8,image_cropped*2^8,8,0);
scielab_ref = getscielab(source_cropped*2^8,image_demosaic_ref_cropped*2^8,8,0);

%% Display
if(display_on ==1)
figure;
imshow(source_cropped);
figure;
imshow(image_cropped);
figure;
imshow(image_demosaic_ref_cropped);
end
