%INITIALIZATION
%rawToMask.m is intended to generate binary masks from single-slice
%multispectral images. To begin, all .tif stacks should be in the same
%folder. The images will be saved in directory 'fileName' as cellmask.tif,
%ch01.tif, ch02.tif, etc.
clear variables
close all
commandwindow;


%MASKING PARAMETER SETUP
numCh = str2double(input('How many channels are there? \n','s'));
filePath = uigetdir('Select folder with images to mask');
addpath(filePath);
fileList = dir([filePath '/*.tif']);


%IMPORT IMAGES
numFiles = size(fileList,1);
images = cell(numFiles,numCh);
for i = 1:numFiles
    for j = 1:numCh
        images{i,j} = imread(fileList(i).name,j);
    end
end


%MANUAL ROI SELECTION
maskCh = str2double(input(['Of the ' num2str(numCh) ' channels, which should be used for masking the cell? \n'], 's'));
roi = cell(1,numFiles);
for i = 1:numFiles
    brightImg = imadjust(images{i,maskCh},[0 0.2]);
    imshow(brightImg);
    roi{i} = roipoly;
    if isempty(roi{i})
        try
            roi{i} = imbinarize(ones(size(images{3},1),size(images{3},2)));
        catch
            roi{i} = imbinarize(ones(size(images{1},1),size(images{1},2)));
        end
    end
end


%MASK CHANNELS
close(figure(1));
for i = 1:numFiles
    mkdir(fullfile(filePath, strrep(fileList(i).name,'.tif','')));
    savePath = fullfile(filePath, strrep(fileList(i).name,'.tif',''));    
    %Apply ROI to image and blur
    roi{i} = immultiply(roi{i},images{i,maskCh});
    roi{i} = medfilt2(roi{i}, [15 15]);
    
    %Calculate threshold from blurred image
    t = double(multithresh(roi{i},3));
    t = t(1)/65535;
    
    %Make binary, fill holes, exclude small objects
    mask = bwareaopen(imfill(imbinarize(roi{i},t),'holes'),150);
    
    %Save mask
    if numFiles > 1
        imwrite(mask,fullfile(savePath, 'cellmask.tif'),'compression','none');
    else
        imwrite(mask,fullfile(savePath,'cellmask.tif'),'compression','none');
    end
    
    %Mask all other channels
    for j = 1:numCh
        images{i,j} = images{i,j}.*uint16(mask);
    end
end


%THRESHOLD PARAMETER SETUP
segIndex = zeros(1,numCh);
for i = 1:numCh
    %Ask how to segment other channels
    imshow(imadjust(images{1,i},[0 0.2]));
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
B = ones(50,50)/(50^2);
for i = 1:numFiles
    disp(['Working on cell ' num2str(i) '...']);
    savePath = fullfile(filePath, strrep(fileList(i).name,'.tif',''));
    
    for j = 1:numCh
        %Particle segmentation
        if segIndex(j) == 0
            %Gaussian filter and opening
            C = uint16(conv2(images{i,j},B,'same'));
            images{i,j} = imgaussfilt(imsubtract(images{i,j},C),1.5);
            t = double(multithresh(images{i,j},2));
            t = t(2)/65535;
            images{i,j} = imbinarize(images{i,j},t);
            images{i,j} = imopen(images{i,j},ball);
            
            %Watershed
            D = -bwdist(~images{i,j});
            M = imextendedmin(D,2);
            D2 = imimposemin(D,M);
            L = watershed(D2,8);
            L(~images{i,j}) = 0;
            images{i,j} = bwareaopen((L>0),10);
            
        %Large organelle segmentation
        elseif segIndex(j) == 1
            %Median filter and opening
            images{i,j} = medfilt2(images{i,j});
            t = double(multithresh(images{i,j},2));
            t = t(2)/65535;
            images{i,j} = imbinarize(images{i,j},t);
            images{i,j} = imopen(images{i,j},ball);
            
        %Sheet segmentation
        else
            images{i,j} = medfilt2(images{i,j});
            t = adaptthresh(images{i,j}, 'NeighborhoodSize', 51);
            images{i,j} = imbinarize(images{i,j},t);
        end
        
        imwrite(images{i,j},fullfile(savePath,['ch0' num2str(j) '.tif']),'Compression','none');
    end
end

disp('Complete.');



