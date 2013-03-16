% Test our demosaicing algorithm with a small image.
% input image soulbe be RGB image
clear;

%% Tunable Paramters
GtoRed = 0.1; %effectiveness of Green for Red Demosaicing 
GtoBlue = 0.1;

source = im2double(imread('../testimages/rock.jpg'));
%source = source(1:10,1:10,:);
img_size = size(source);
img_size = img_size(1:2);


%% Create the artificial Bayer matrix

bayerPattern = [1, 2; 2, 3]; % 1=red, 2=green, 3=blue
bayerFull = false([img_size,3]);


%% Full Bayer Pattern generation
for c = 1:3
  bayerFull(:, :, c) = repmat(bayerPattern == c, img_size / 2);
end
  bayerFull_v = repmat(bayerPattern, img_size / 2);

%% Apply Bayer Pattern to the source image
for c = 1:3
image_mosaic(:,:,c) = source(:,:,c) .* double(bayerFull(:,:,c));
end

image_mosaic_RGBcombined = image_mosaic(:,:,1) + image_mosaic(:,:,2) + image_mosaic (:,:,3);


%% Interpolation
%% 1. interpolate green channel with second-derivative weighting
%% Gx2 = Bx-2 + Bx+2 -2*Bx

% Blue grandiednt matrix
padd = [0,2]
Bx_n2 = padarray(image_mosaic (:,:,3),padd,0,'pre');
Bx_p2 = image_mosaic(:,3:(img_size(2)),3);
Bx_p2 = padarray(Bx_p2,[0,4],0,'post');
Bx = padarray(image_mosaic (:,:,3),padd,0,'post');
Gbx2 = Bx_n2 + Bx_p2 -2 *Bx;

padd = [2,0]
By_n2 = padarray(image_mosaic (:,:,3),padd,0,'pre');
By_p2 = image_mosaic(3:(img_size(1)),:,3);
By_p2 = padarray(By_p2,[4,0],0,'post');
By = padarray(image_mosaic (:,:,3),padd,0,'post');
Gby2 = By_n2 + By_p2 -2 *By;

% Red grandiednt matrix
padd = [0,2]
Rx_n2 = padarray(image_mosaic (:,:,1),padd,0,'pre');
Rx_p2 = image_mosaic(:,3:(img_size(2)),1);
Rx_p2 = padarray(Rx_p2,[0,4],0,'post');
Rx = padarray(image_mosaic (:,:,1),padd,0,'post');
Grx2 = Rx_n2 + Rx_p2 -2 *Rx;

padd = [2,0]
Ry_n2 = padarray(image_mosaic (:,:,1),padd,0,'pre');
Ry_p2 = image_mosaic(3:(img_size(1)),:,1);
Ry_p2 = padarray(Ry_p2,[4,0],0,'post');
Ry = padarray(image_mosaic (:,:,1),padd,0,'post');
Gry2 = Ry_n2 + Ry_p2 -2 * Ry;


%% Filling in Green value
for r = 3: (img_size(1)-3)
  for c = 3: (img_size(2)-3)
    
    if(bayerFull_v(r,c) == 1) %R  
     blur_r = image_mosaic(r-1,c,2) + image_mosaic(r+1,c,2);
     blur_c = image_mosaic(r,c-1,2) + image_mosaic(r,c+1,2);
     
     if (Grx2(r,c) == 0 && Gry2(r,c) == 0)
         r_weight = 0.5;
         c_weight = 0.5;
     else
       r_weight = abs(Grx2(r,c))/(abs(Grx2(r,c)) + abs(Gry2(r,c)));
       c_weight = abs(Gry2(r,c))/(abs(Grx2(r,c)) + abs(Gry2(r,c)));
     end 

     image_mosaic(r,c,2) = 0.5*r_weight*blur_r + 0.5*c_weight*blur_c;

    elseif(bayerFull_v(r,c) == 3) %B  
     blur_r = image_mosaic(r-1,c,2) + image_mosaic(r+1,c,2);
     blur_c = image_mosaic(r,c-1,2) + image_mosaic(r,c+1,2);
     
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
gerr_avg= sum(sum(g_err,1),2)/(size(g_err,1))/ (size(g_err,2))

%% Filling in R,B holes by two stages
%% First Pass - fill the holes between Reds
for r = 2: (img_size(1)-1)
  for c = 2: (img_size(2)-1)
    
    %% Filling in Red value
    if(bayerFull_v(r,c) == 2 && bayerFull_v(r,c-1) == 1) %RGR  
     if(image_mosaic(r,c-1,2) ~= 0 && image_mosaic(r,c+1,2) ~=0)
       g_factor = 0.5*min((image_mosaic(r,c-1,1)/image_mosaic(r,c-1,2) + image_mosaic(r,c+1,1)/image_mosaic(r,c+1,2)),1.0);
     else 
       g_factor =0.5;
     end
     image_mosaic(r,c,1) = (1.0-GtoRed)*0.5*(image_mosaic(r,c-1,1) + image_mosaic(r,c+1,1)) + GtoRed*image_mosaic(r,c,2)*g_factor;
    elseif(bayerFull_v(r,c) == 2 && bayerFull_v(r-1,c) == 1) %R;G;R  
     if(image_mosaic(r-1,c,2) ~= 0 && image_mosaic(r+1,c,2) ~=0)
       g_factor = 0.5*min((image_mosaic(r-1,c,1)/image_mosaic(r-1,c,2) + image_mosaic(r+1,c,1)/image_mosaic(r+1,c,2)),1.0);
     else 
       g_factor =0.5;
     end
     image_mosaic(r,c,1) = (1.0-GtoRed)*0.5*(image_mosaic(r-1,c,1) + image_mosaic(r+1,c,1)) + GtoRed*image_mosaic(r,c,2)*g_factor;
    end

    %% Filling in Blue value
    if(bayerFull_v(r,c) == 2 && bayerFull_v(r,c-1) == 3) %BGB  
     if(image_mosaic(r,c-1,2) ~= 0 && image_mosaic(r,c+1,2) ~=0)
       g_factor = 0.5*min((image_mosaic(r,c-1,3)/image_mosaic(r,c-1,2) + image_mosaic(r,c+1,3)/image_mosaic(r,c+1,2)),1.0);
     else 
       g_factor =0.5;
     end
     image_mosaic(r,c,3) = (1.0-GtoBlue)*0.5*(image_mosaic(r,c-1,3) + image_mosaic(r,c+1,3)) + GtoBlue*image_mosaic(r,c,2)*g_factor;
    elseif(bayerFull_v(r,c) == 2 && bayerFull_v(r-1,c) == 3) %B;G;B  
     if(image_mosaic(r-1,c,2) ~= 0 && image_mosaic(r+1,c,2) ~=0)
       g_factor = 0.5*min((image_mosaic(r-1,c,3)/image_mosaic(r-1,c,2) + image_mosaic(r+1,c,3)/image_mosaic(r+1,c,2)),1);
     else 
       g_factor =0.5;
     end
     image_mosaic(r,c,3) = (1.0-GtoRed)*0.5*(image_mosaic(r-1,c,3) + image_mosaic(r+1,c,3)) + GtoRed*image_mosaic(r,c,2)*g_factor;
    end

  end
end


%% Second Pass
%% First Pass - fill the rest of holes (bilinear interpolation)

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

figure;
imshow(source);
figure;
imshow(image_mosaic);

img_demosaic_ref = demosaic(uint8(image_mosaic_RGBcombined *2^8),'rggb');
figure;
imshow(img_demosaic_ref);