%INITIALIZATION
%rawZStackToMask.m is intended to generate binary masks from single-slice
%multispectral images. To begin, all .tif stacks should be in the same
%folder. The images will be saved in directory 'fileName' as cellmask.tif,
%ch01.tif, ch02.tif, etc.
%It is HIGHLY recommended that you do not change the filenames resulting
%from maskMultispectral - doing so may compromise your ability to use
%countMultispectral and interactomeMultispectral. 
clear variables
close all
commandwindow;


%MASKING PARAMETER SETUP
numCh = str2double(input('How many channels are there? \n','s'));
filePath = uigetdir('Select folder with images to mask');
addpath(filePath);
fileList = dir([filePath '/*.tif']);
%temp = imread(fileList(i).name,1);


%IMPORT IMAGES
numFiles = size(fileList,1);
numSlices = zeros(1,numFiles);
%images = cell(numFiles,numCh);
for i = 1:numFiles
    z = 1;
    while z > 0
        try
            j = rem(z-1,numCh)+1;
 %           images{i,j,ceil(z/2)} = imread(fileList(i).name,z);
            images(:,:,i,j,ceil(z/numCh)) = imread(fileList(i).name,z);
            z = z+1;
        catch
            numSlices(i) = (z-1)/numCh;
            z = 0;
        end
    end
    disp(i);
end



%MANUAL ROI SELECTION
maskCh = str2double(input(['Of the ' num2str(numCh) ' channels, which should be used for masking the cell? \n'], 's'));
roi = cell(1,numFiles);
for i = 1:numFiles
    brightImg = imadjust(max(images(:,:,i,maskCh,:),[],5),[0 0.2]);
    imshow(brightImg);
    roi{i} = roipoly;
    if isempty(roi{i})
        try
            roi{i} = imbinarize(ones(size(images(:,:,3),1),size(images(:,:,3),2)));
        catch
            roi{i} = imbinarize(ones(size(images(:,:,1),1),size(images(:,:,1),2)));
        end
    end
end


%MASK CHANNELS
close(figure(1));
maskType = str2double(input('Use selected mask channel for tighter mask? 0 for Yes, 1 for No. Recommended only if using ER as mask channel. \n','s'));
for i = 1:numFiles
    mkdir(fullfile(filePath, strrep(fileList(i).name,'.tif','')));
    savePath = fullfile(filePath, strrep(fileList(i).name,'.tif',''));    

    %Make binary, fill holes, exclude small objects
%     mask = bwareaopen(imfill(imbinarize(roi{i},t),'holes'),150);
    if maskType == 1
        mask = roi{i};
    else
        mask = max(images(:,:,i,maskCh,:),[],5);
        try
            mask = medfilt2(mask,[15 15]).*uint16(roi{i});
            disp('uint16');
        catch
            mask = medfilt2(mask,[15 15]).*uint8(roi{i});
            disp('uint8');
        end
        t = adaptthresh(mask,.9, 'NeighborhoodSize', 51);
        mask = imfill(imclose(imbinarize(mask,t),strel('disk',1,0)),'holes');
    end
    
    %Save mask
    if numFiles > 1
        imwrite(mask,fullfile(savePath, 'cellmask.tif'),'Compression','none');
    else
        imwrite(mask,fullfile(savePath,'cellmask.tif'),'Compression','none');
    end
    
%     Mask all other channels
    for j = 1:numCh
        for z = 1:numSlices(i)
            try
                images(:,:,i,j,z) = images(:,:,i,j,z).*uint16(mask);
            catch
                images(:,:,i,j,z) = images(:,:,i,j,z).*uint8(mask);
            end
        end
    end
end


%THRESHOLD PARAMETER SETUP
segIndex = zeros(1,numCh);
for i = 1:numCh
    %Ask how to segment other channels
    imshow(images(:,:,1,i,1));
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
mipImages = cell(numFiles,numCh);
tf = isa(images,'uint8');
if tf == 1
    div = 255;
else
    div = 65535;
end
for i = 1:numFiles
    disp(['Working on cell ' num2str(i) '...']);
    savePath = fullfile(filePath, strrep(fileList(i).name,'.tif',''));
%     figure(i);              %Important: Only goes up to 6 channels!
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     hold on;
    for j = 1:numCh
            mipImages{i,j} = max(images(:,:,i,j,:),[],5);
            if tf == 1
                C = uint8(conv2(mipImages{i,j},B,'same'));
            else
                C = uint16(conv2(mipImages{i,j},B,'same'));
            end
            %Particle segmentation
            if segIndex(j) == 0
                %Gaussian filter and opening
                mipImages{i,j} = imgaussfilt(imsubtract(mipImages{i,j},C),1.5);
                t = double(multithresh(mipImages{i,j},2));
                t = t(2)/div;
                disp(t);
                mipImages{i,j} = imbinarize(mipImages{i,j},t);      
                mipImages{i,j} = imopen(mipImages{i,j},ball);

                %Watershed
                D = -bwdist(~mipImages{i,j});
                M = imextendedmin(D,2);
                D2 = imimposemin(D,M);
                L = watershed(D2,8);
                L(~mipImages{i,j}) = 0;
                mipImages{i,j} = bwareaopen(imbinarize(L,0),10);
%                 bw = bwdist(~mipImages{i,j});
%                 bw = -bw;
%                 bw(~mipImages{i,j}) = Inf;
%                 bw = watershed(bw,8);
%                 bw(~mipImages{i,j}) = 0;
%                 mipImages{i,j} = ~bw;            

            %Large organelle segmentation
            elseif segIndex(j) == 1
                %Median filter and opening
                mipImages{i,j} = medfilt2(imsubtract(mipImages{i,j},C));
                t = double(multithresh(mipImages{i,j},2));
                t = t(2)/div;
                mipImages{i,j} = imbinarize(mipImages{i,j},t);
                mipImages{i,j} = imopen(mipImages{i,j},ball);

            %Sheet segmentation
            else
                mipImages{i,j} = medfilt2(imsubtract(mipImages{i,j},C));
                t = adaptthresh(mipImages{i,j}, 'NeighborhoodSize', 51);
                mipImages{i,j} = imbinarize(mipImages{i,j},t);
            end

            
            
            
        %Save image
%         if numFiles > 1
%             imwrite(images{i,j},fullfile(savePath,...
%                 ['ch0' num2str(j) '.tif']));
%         else
        imwrite(mipImages{i,j},fullfile(savePath,['ch0' num2str(j) '.tif']),'Compression','none');
%         end
    end
end

disp('Complete.');



