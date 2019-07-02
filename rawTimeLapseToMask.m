%rawTimelapseToMask
clear variables
close all
commandwindow;


%MASKING PARAMETER SETUP
numCh = str2double(input('How many channels are there? \n','s'));
numFrames = str2double(input('How many frames are there? \n','s'));
filePath = uigetdir('Select folder with images to mask');
addpath(filePath);
fileList = dir([filePath '/*.tif']);


%IMPORT IMAGES
numFiles = size(fileList,1);
images = cell(numFiles,numCh);
for i = 1:numFiles
    k = 0;
    h = 1;
    for j = 1:numCh*numFrames
        k = k+1;
        images{i,k,h} = imread(fileList(i).name,j);
        if rem(k,numCh) == 0
            k = 0;
            h = h+1;
        end
    end
end


%MANUAL ROI SELECTION
maskCh = str2double(input(['Of the ' num2str(numCh) ' channels, which should be used for masking the cell? \n'], 's'));
roi = cell(1,numFiles);
for i = 1:numFiles
    brightImg = imadjust(images{i,maskCh,1},[0 0.2]);
    imshow(brightImg);
    roi{i} = roipoly;
    if isempty(roi{i})
        try
            roi{i} = imbinarize(ones(size(images{1,3,1},1),size(images{1,3,1},2)));
        catch
            roi{i} = imbinarize(ones(size(images{1,1,1},1),size(images{1,1,1},2)));            
        end
        
    end
end


%MASK CHANNELS
close(figure(1));
for i = 1:numFiles
    mkdir(fullfile(filePath, strrep(fileList(i).name,'.tif','')));
    savePath = fullfile(filePath, strrep(fileList(i).name,'.tif',''));    
    %Apply ROI to image and blur
    roi{i} = immultiply(roi{i},images{i,maskCh,1});
    roi{i} = medfilt2(roi{i}, [15 15]);
    
    %Calculate threshold from blurred image
    t = double(multithresh(roi{i},3));
    t = t(1)/65535;
    
    %Make binary, fill holes, exclude small objects
    mask = bwareaopen(imfill(imbinarize(roi{i},t),'holes'),150);
    
    %Save mask
    if numFiles > 1
        imwrite(mask,fullfile(savePath, 'cellmask.tif','compression','none'));
    else
        imwrite(mask,fullfile(savePath,'cellmask.tif','compression','none'));
    end
    
    %Mask all other channels
    for j = 1:size(images,2)
        for k = 1:size(images,3)
            images{i,j,k} = images{i,j,k}.*uint16(mask);
        end
    end
end


%THRESHOLD PARAMETER SETUP
segIndex = zeros(1,numCh);
for i = 1:numCh
    %Ask how to segment other channels
    imshow(imadjust(images{1,i,1},[0 0.2]));
    response = questdlg(['What is the morphology of channel ' num2str(i) '? Close to use defaults.'],...
        'Channel segmentation workflow selection','Particle',...
        'Large organelle','Sheet','Large organelle');
    %Assign segmentation index
    switch response
        case 'Particle'
            segIndex(i) = 0;
        case 'Large organelle'
            segIndex(i) = 1;
        case 'Sheet'
            segIndex(i) = 2;
        otherwise
            if numCh == 6 
                segIndex = [0 1 2 0 1 0];
            else
                segIndex = zeros(1,numCh);
            end
            break
    end
end


%THRESHOLD AND SAVE
close(figure(1));
ball = strel('disk',1,0);
for i = 1:numFiles
    disp(['Working on cell ' num2str(i) '...']);
    savePath = fullfile(filePath, strrep(fileList(i).name,'.tif',''));
%     figure(i);              %Important: Only goes up to 6 channels!
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     hold on;
    for j = 1:size(images,2)
        for k = 1:size(images,3)
            %Particle segmentation
            if segIndex(j) == 0
                %Gaussian filter and opening
                images{i,j,k} = imgaussfilt(images{i,j,k});
                t = double(multithresh(images{i,j,k},2));
                t = t(2)/65535;
                images{i,j,k} = imbinarize(images{i,j,k},t);
                images{i,j,k} = imopen(images{i,j,k},ball);

                %Watershed
                D = -bwdist(~images{i,j,k});
                M = imextendedmin(D,2);
                D2 = imimposemin(D,M);
                L = watershed(D2,8);
                L(~images{i,j,k}) = 0;
                images{i,j,k} = bwareaopen((L>0),10);            

            %Large organelle segmentation
            elseif segIndex(j) == 1
                %Median filter and opening
                images{i,j,k} = medfilt2(images{i,j,k});
                t = double(multithresh(images{i,j,k},2));
                t = t(2)/65535;
                images{i,j,k} = imbinarize(images{i,j,k},t);
                images{i,j,k} = imopen(images{i,j,k},ball);

            %Sheet segmentation
            else
                images{i,j,k} = medfilt2(images{i,j,k});
                t = adaptthresh(images{i,j,k}, 'NeighborhoodSize', 51);
                images{i,j,k} = imbinarize(images{i,j,k},t);
            end


            %Save image
            if k == 1
            imwrite(images{i,j,k},...
                fullfile(savePath,['ch0' num2str(j) '.tif']),'Compression','none');
            else
                imwrite(images{i,j,k},fullfile(savePath,['ch0' num2str(j) '.tif']),...
                    'Compression','none','WriteMode','append');
            end
        end
    end
end

disp('Complete.');



