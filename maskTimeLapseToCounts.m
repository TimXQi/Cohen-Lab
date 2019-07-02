%INITIALIZATION
%maskTimeLapseToCounts.m allows users to select folders created by
%maskMultispectral.m. These folders should contain only *ch0#.tif and
%*cellmask.tif images. The user is prompted to desginate the organelle type
%in each channel, which affects the computed parameters. 

clear variables
close all
commandwindow;

%FOLDER SELECTION
numFrames = str2double(input('How many frames are there? \n','s'));
folderList = uigetdir2(0,'Select folders to compute');
addpath(folderList{:});
numFolders = size(folderList,2);
numCh = size(dir([folderList{1} '/*ch*.tif']),1);
images = cell(numFolders,numCh+1);
pixSize = str2double(input('How many pixels per micron? \n','s'))^2;


%READ IN IMAGES
for i = 1:numFolders
    fileList = dir([folderList{i} '/*.tif']);
    for j = 1:numCh+1
        if j ==1
            images{i,j,1} = imread(fullfile(folderList{i},fileList(j).name));
        else
            for k = 1:numFrames
                images{i,j,k} = imread(fullfile(folderList{i},fileList(j).name),k);
            end
        end
    end
end
orgIndex = zeros(1,numCh+1); orgIndex(1) = -1;
orgNames = cell(1,numCh+1); orgNames{1} = 'mask';

%ORGANELLE DESIGNATION
for i = 2:numCh+1
    %Ask how to segment other channels
    %Skip mask in slot 1, assign -1
    imshow(images{1,i});
    if isempty(orgNames{i}) == 1
        try
            orgNames(i) = inputdlg('What will you name this channel? Close to use defaults.');
        catch
            if numCh == 6
                orgIndex = [-1 0 1 2 0 1 0];
                orgNames = {'mask','lysosome','mitochondria','ER','peroxisome','Golgi','lipid droplet'};
            else
                orgIndex = zeros(1,numCh);
            end
            break
        end
    end
    response = questdlg(['What is the morphology of channel ' num2str(i-1) '?'],...
        'Channel segmentation workflow selection','Particle',...
        'Large organelle','Sheet','Large organelle');
    %Assign segmentation index
    switch response
        case 'Particle'
            orgIndex(i) = 0;
        case 'Large organelle'
            orgIndex(i) = 1;
        case 'Sheet'
            orgIndex(i) = 2;
    end
end
close(figure(1));

%GENERATE RESULT MATRIX
numParam = zeros(1,4);
for i = -1:2
    numParam(i+2) = sum(orgIndex == i);
end
numParam = sum(numParam.*[1 4 4 2]);
results = cell(numFolders,numParam+1,numFrames+1);
for i = 1:size(results,1)
    index = 2;
    results{i,1,1} = fliplr(strtok(fliplr(folderList{i}),'\'));
    for j = 1:length(orgIndex)
        if orgIndex(j) == -1
            results{i,index,1} = 'Surface area';
            index = index+1;
        elseif orgIndex(j) == 2
            results{i,index,1} = [orgNames{j} ' area'];
            results{i,index+1,1} = [orgNames{j} ' area fraction'];
            index = index+2;
        else
            results{i,index,1} = [orgNames{j} ' number'];
            results{i,index+1,1} = [orgNames{j} ' mean area'];
            results{i,index+2,1} = [orgNames{j} ' median area'];
            results{i,index+3,1} = [orgNames{j} ' area fraction'];
            index = index+4;
        end
    end

    for j = 2:size(results,3)
        results{i,1,j} = ['frame ' num2str(j-1)];
    end
end
% COMPUTATION
for i = 1:numFolders       %new loop every cell
    index = 2;
    for j = 1:numCh+1             %new loop every channel
        for k = 2:numFrames+1
            if orgIndex(j) == -1 && k == 2   %mask
                results{i,index,k} = bwarea(images{i,j,k-1})/pixSize;

            elseif orgIndex(j) == 2 %sheet
                results{i,index,k} = bwarea(images{i,j,k-1})/pixSize;
                if sum(orgIndex == -1) == 1
                    maskIndx = find(orgIndex==-1);
                    results{i,index+1,k} = bwarea(images{i,j,k-1})/bwarea(images{i,maskIndx,1});
                else
                    results{i,index+1,k} = 'Cannot compute without mask';
                end

            else                    %particles and large organelles
                bw = bwconncomp (images{i,j,k-1});
                results{i,index,k} = bw.NumObjects;

                areas = zeros(1,bw.NumObjects);
                for l = 1:bw.NumObjects
                    areas(l) = length(bw.PixelIdxList{l});
                end
                results{i,index+1,k} = mean(areas)/pixSize;
                results{i,index+2,k} = median(areas)/pixSize;


                if sum(orgIndex == -1) == 1
                    maskIndx = find(orgIndex==-1);
                    results{i,index+3,k} = bwarea(images{i,j,k-1})/bwarea(images{i,maskIndx,1});
                else
                    results{i,index+3,k} = 'Cannot compute without mask';
                end
            end
        end
        if orgIndex(j) == -1
            index = index+1;
        elseif orgIndex(j) == 2
            index = index+2;
        else
            index = index+4;
        end
    end
end


%SAVE RESULT
for i = 1:numFolders
    writetable(cell2table(squeeze(results(i,:,:))),[fliplr(strrep(fliplr(...
        folderList{1}),strtok(fliplr(folderList{1}),'\'),'')) 'timelapse ' num2str(i) ' organelle counts.csv'],...
        'WriteVariableNames',0);
end

disp('Complete. Results salved in parent folder.');
%MULTIFOLDER SELECT FUNCTION
function [pathname] = uigetdir2(start_path, dialog_title)
    import javax.swing.JFileChooser;

    if nargin == 0 || start_path == 0 
        start_path = pwd;
    end

    jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);
    jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

    if nargin > 1
        jchooser.setDialogTitle(dialog_title);
    end

    jchooser.setMultiSelectionEnabled(true);
    status = jchooser.showOpenDialog([]);

    if status == JFileChooser.APPROVE_OPTION
        jFile = jchooser.getSelectedFiles();
        pathname{size(jFile, 1)}=[];
        for i=1:size(jFile, 1)
            pathname{i} = char(jFile(i).getAbsolutePath);
        end

    elseif status == JFileChooser.CANCEL_OPTION
        pathname = [];
    else
        error('Error occured while picking file.');
    end
end
