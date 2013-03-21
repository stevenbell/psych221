%% Sample use
% clear all;
% 
% width = 5;
% height = 5;
% 
% nSamples = 4;
% x = [0, 0.1, 0  , 0.3];
% y = [0, 0  , 0.2, 0.5];
% 
% NEW_SCALE = 0.6;
% 
% data_R = cell(nSamples,1);
% for ii=1:nSamples
%     %%Create array with 3 elements out of 4 being 0
%     data_R{ii}  = insertZeros(rand(height,width),2);     
% end
% 
% a = CombineLR( data_R, [1 0; 0 0], x, y , NEW_SCALE);
% 
% mesh(a);
%%


function [ res ] = CombineLR( IMG, PATTERN, x, y , NEW_SCALE)
    SHOW_IMG = true;%false;

    nSamples = size(IMG,1);

    [height width] = size(IMG{1});

    %width = 5;
    %height = 5;

    %nSamples = 4;
    %x = [0, 0.1, 0  , 0.3];
    %y = [0, 0  , 0.2, 0.5];

    data_R = cell(nSamples,1);
    X = cell(nSamples,1);
    Y = cell(nSamples,1);

    data_R_s = cell(nSamples,1);
    X_s = cell(nSamples,1);
    Y_s = cell(nSamples,1);

    if SHOW_IMG
        figure1 = figure;

        axes1 = axes('Parent',figure1,'XTick',1:width,'XMinorTick','on', ...
           'YTick',1:width,'YMinorTick','on');
        view(axes1,[0.5 90]);
        grid(axes1,'on');
        hold(axes1,'all');
    end

    SCALE=10;
    PITCH=1/SCALE;

    Mask = zeros(SCALE*height,SCALE*width);
    data_Sum  = zeros(SCALE*height,SCALE*width);

    %PATTERN = [1 0; 0 1];%Green

    Data_Pattern = repmat(PATTERN, ceil(height/2), ceil(width/2));
    Data_Pattern = Data_Pattern(1:height,1:width);

    
    for ii=1:nSamples
        data_R{ii}  = IMG{ii};%insertZeros(rand(height,width),2);  

        data_Sum = data_Sum + remapImage( data_R{ii}, SCALE, SCALE*x(ii), SCALE*y(ii) );

        Mask = Mask + remapImage( Data_Pattern, SCALE, SCALE*x(ii), SCALE*y(ii) );

        data_R_s{ii} = remapImage( data_R{ii}, SCALE, 0, 0);


        [X{ii} Y{ii}] = meshgrid(1+x(ii):width+x(ii), 1+y(ii):height+y(ii));
        [X_s{ii} Y_s{ii}] = meshgrid(1:PITCH:width+1-PITCH, 1:PITCH:height+1-PITCH);
        X_s{ii} = X_s{ii} + x(ii);
        Y_s{ii} = Y_s{ii} + y(ii);    
    end


    %remove empty data
    for ii=1:nSamples
        X{ii}(data_R{ii}==0)=[];
        Y{ii}(data_R{ii}==0)=[];
        data_R{ii}(data_R{ii}==0)=[];

        X_s{ii}(data_R_s{ii}==0)=[];
        Y_s{ii}(data_R_s{ii}==0)=[];
        data_R_s{ii}(data_R_s{ii}==0)=[];
        
        if SHOW_IMG
            scatter3(X{ii}(:), Y{ii}(:),data_R{ii}(:),'LineWidth',1.7); hold on
        end
    end

    %% Map to 10x scale and interpolate

    [X_f Y_f] = meshgrid(1:PITCH:width+1-PITCH, 1:PITCH:height+1-PITCH);

    [y_c x_c] = find(Mask==0);

    T = data_Sum ./ Mask;
    T(find(Mask==0)) = [];
    X_sample = X_f(find(Mask~=0));
    Y_sample = Y_f(find(Mask~=0));

    R_super = griddata( X_sample , Y_sample, T ,X_f,Y_f, 'cubic');

    if SHOW_IMG
        %mesh(X_f, Y_f,R_super);
        contour(X_f, Y_f,R_super);
        alpha(0.1);
        %title('LR samples and interpolated data');
        [X_sample Y_sample] = meshgrid(1:NEW_SCALE:width+1, 1:NEW_SCALE:height+1);
        mesh(X_sample, Y_sample,X_sample*0); hold on
        alpha(0.1);
        
        
        figure;
        mesh(X_f, Y_f,Mask);
        %title('Number of samples per point');
    end

    %% Downsize the d
    %NEW_SCALE = 1;%0.95;

    [X_res Y_res] = meshgrid(1:NEW_SCALE:width+1, 1:NEW_SCALE:height+1);
    d_res = R_super(1:SCALE*NEW_SCALE:SCALE*height,1:SCALE*NEW_SCALE:SCALE*width);

    [a0 b0] = size(d_res);
    [a1 b1] = size(X_res);
    
    a=min(a0,a1);
    b=min(b0,b1);
    
    if SHOW_IMG
        figure;
        mesh(X_res(1:a,1:b), Y_res(1:a,1:b),d_res(1:a,1:b)); hold on;
        alpha(0.1);
        
        x_1 = X_res(1:a,1:b);
        y_1 = Y_res(1:a,1:b);
        z_1 = d_res(1:a,1:b);
        
        
        for ii=1:nSamples
            scatter3(X{ii}(:), Y{ii}(:),data_R{ii}(:),'LineWidth',1.7); hold on
        
        end
        
        
        plot3(x_1(:), y_1(:),z_1(:),'o','LineWidth',1.5);
        hold off;
        %title('Result: subsampled interpolated data');
    
    
        %% Display alligned data from LR images

        figure2=figure;
        
        axes1 = axes('Parent',figure2);
        view(axes1,[0.5 90]);
        grid(axes1,'on');
        hold(axes1,'all');


        for ii=1:nSamples
            scatter3(X_s{ii}(:), Y_s{ii}(:),data_R_s{ii}(:),'LineWidth',1.5);hold on
        end


        %NEW_SCALE = 0.8;

        [X_sample Y_sample] = meshgrid(1:NEW_SCALE:width+2, 1:NEW_SCALE:height+2);
        mesh(X_sample, Y_sample,X_sample*0); hold on
        alpha(0.1);
        
        %title('Map LR samples to squares');
        
        hold off;

        %% Display samples and weight function on the new grid
        figure5 = figure;
        
        axes1 = axes('Parent',figure5);
        view(axes1,[0.5 90]);
        grid(axes1,'on');
        hold(axes1,'all');
        
        

        for ii=1:nSamples
            scatter3(X{ii}(:), Y{ii}(:),data_R{ii}(:),'LineWidth',1.5); hold on
            %mesh(X{ii}, Y{ii},data_R{ii}); hold on
        end


        %NEW_SCALE = 0.8;

        [X_sample Y_sample] = meshgrid(1:NEW_SCALE:width+1, 1:NEW_SCALE:height+1);
        mesh(X_sample, Y_sample,X_sample*0); hold on
        alpha(0.1);

        center_x =  X_sample(1,round((width+1)/NEW_SCALE/2));
        center_y =  Y_sample(round((height+1)/NEW_SCALE/2),1);


        [x y] = meshgrid(-2+center_x:0.1:2+center_x, -2+center_y:0.1:2+center_y);
        weight = exp(-((y-center_y).^2 + (x-center_x).^2));
        contour(x,y,weight); hold on;
        %mesh(x,y,weight); hold on;
        %surfc(x,y,weight); hold on;
        %shading interp;
        alpha(.1);
        
        %title('LR samples and weighting function');
        
        hold off;

        
        
        
        figure;
        for ii=1:nSamples
            scatter3(X{ii}(:), Y{ii}(:),data_R{ii}(:),'LineWidth',1.5); hold on
        end
        mesh(X_sample, Y_sample,X_sample*0); hold on
        alpha(0.1);
        mesh(x,y,weight); hold on;
        alpha(.1);
        hold off;
        
        
        
        
    end
    %%
    res = d_res;

end

