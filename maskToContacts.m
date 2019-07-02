%INITIALIZATION
%maskToContacts
%DO NOT CHANGE FILE NAMES AFTER USING maskMULTISPECTRAL
%interactomeMultispectral.m takes the resulting folders of
%maskMultispectral.m and analyzes contact sites between the channels. Any
%number of channels n can be selected, and contacts of any multiplicity up
%to n can be computed.
%Contacts are assumed to be within 1 pixel (see computeContact function).
%You can expand your inclusion threshold by increasing the radius of
%strel('disk',r) being fed into imdilate.
clear variables
close all
commandwindow;

%FOLDER SELECTION
folderList = uigetdir2(0,'Select folders to compute');
addpath(folderList{:});
numFolders = size(folderList,2);

%WARNING: uigetfile has issues with long file names. See comment in line 3.
chList = uigetfile(fullfile(folderList{1},'/*ch*.tif'),...
    'Select channels to compute contacts for.','Multiselect','on');
numCh = size(chList,2);
chListNum = zeros(1,numCh);
for i = 1:numCh           
    %searches for integers in channel i and returns their location(s)
    temp = regexp(chList{i},'\d');
    %extremely bad section heavily dependent on filename
    %literally: temp looks something like [3 4]
    %takes characters 3 and 4 from chList{i}, which looks like 'ch04.tif'
    %returns '04' --str2double--> 4
    chListNum(i) = str2double(chList{i}(temp(1):temp(2))); 
end
fileList = dir([folderList{1} '/*ch*.tif']);
images = cell(1,numCh);

%NAME ASSIGNMENT
defaultNames = cell(1,numCh);
for i = 1:numCh
    realCh = strrep(fliplr(strtok(fliplr(chList{i}),'_')),'.tif','');
    switch realCh
        case 'ch01'
            defaultNames{i} = ' lysosome ';
        case 'ch02'
            defaultNames{i} = ' mitochondria ';
        case 'ch03'
            defaultNames{i} = ' ER ';
        case 'ch04'
            defaultNames{i} = ' peroxisome ';
        case 'ch05'
            defaultNames{i} =  ' Golgi ';
        case 'ch06'
            defaultNames{i} =  ' lipid droplet ';
        otherwise
            defaultNames{i} = ' organelle ';
    end
end
orgNames = cell(1,numCh);
for i = 1:numCh
    imshow(imread(fullfile(folderList{1},chList{i})));
    response = inputdlg(...
        ['What will you name this channel? Close at any point to use default array of'...
        defaultNames{:}]);
    if isempty(response) == 0
        orgNames{i} = response{1};
    else
        orgNames = defaultNames;
        break;
    end
end
close(figure(1));

%DEGREE OF CONTACT SELECTION
degList = 2:numCh;
opts = cell(1,length(degList));
for i = 1:length(degList)
    opts{i} = num2str(degList(i));
end
opts = {cat(1,opts{:})};
[indx,tf] = listdlg('Name','Contact degree','PromptString','How many degrees of contact?',...
    'ListString',opts,'SelectionMode','single','ListSize',[175 75],...
    'OKString','Apply','CancelString','Binary only');
if isempty(indx) == 1 || tf == 0
    contactDegree = 2;
else
    contactDegree = degList(indx);
end

%CONTACT ANALYSIS
results = {''};
for i = 1:size(folderList,2)
    disp(['Working on cell ' num2str(i)]);
    index = 2;                  %Tracks current row of result matrix
    if i == 1
        results{1,1} = 'cell 1';
    else
        results{1,end+1} = ['cell ' num2str(i)];
    end
    results{1,end+1} =  'number';
    results{1,end+1} = 'number fraction';
    results{1,end+1} = 'area fraction';
    
    %Read images from selected channels
    for j = 1:numCh
        fileName = fullfile(folderList{i},fileList(chListNum(j)).name);
        images{j} = imread(fileName);
    end
   
    for j = 2:contactDegree
        %Generates list of unique indices with multiplicity j
        allPerms = npermutek(1:numCh,j);    %do not change
        for k = size(allPerms,1):-1:1   %count down
            %If the number of unique elements in a permutation differs from
            %the total number of elements in a permutation, delete the row.
            %Counts down from top to avoid indexing changes
            if length(unique(allPerms(k,:))) ~= size(allPerms,2) == 1
                allPerms(k,:) = [];
            end
        end
        
        %For multiplicities above 2, the order of indices 2:end do not
        %matter. LD touching [lysosome mitonchondria golgi] is the same as
        %LD touching [mitochondria golgi lysosome] and [lysosome golgi
        %mitochondria].
        if j>2
            %Sorts all permutations of indices 2:end sequentially and eliminates duplicates.
            uniquePerms = unique(sort(allPerms(:,2:end),2),'rows');
            c = 1;
            for k = 1:numCh
                for l = 1:size(uniquePerms,1)
                    %Scans each row of permutations 2:end and eliminates rows
                    %containing index 1. For example, whereas [LD lysosome]
                    %and [mitochondria lysosome] are both unique
                    %permutations, the former will be eliminated when
                    %examinig LD contacts with [organelle1 organelle2].
                    if ismember(k,uniquePerms(l,:)) == 0
                        if c == 1
                            usePerms = [k uniquePerms(l,:)];
                            c = c+1;
                        else
                            usePerms(c,:) = [k uniquePerms(l,:)];
                            c = c+1;
                        end
                    end
                end
            end
        else
            %For binary contacts, all unique permutations are relevant.
            usePerms = allPerms;
        end
        for k = 1:size(usePerms,1)
            %Feeds parameters into computeContacts function and stores
            %results in cell matrix.
            [results{index,(i-1)*4+1}, results{index,(i-1)*4+2}, results{index,(i-1)*4+3}, results{index,(i-1)*4+4}] = ...
                computeContacts(images,orgNames,usePerms(k,:));
            index = index+1;
        end
        
    end
   
end

%SAVE RESULTS
writetable(cell2table(results),[fliplr(strrep(fliplr(...
    folderList{1}),strtok(fliplr(folderList{1}),'\'),'')) 'organelle interactome.csv'],...
    'WriteVariableNames',0);
disp('Complete. Results salved in parent folder.');

%CONTACT COMPUTATION FUNCTION
function [name, number, nFrac, aFrac] = computeContacts(images, names, indices)
    name = [names{indices(1)} ':' names{indices(2:end)}];
    baseArea = bwarea(~images{indices(1)});
    %Establishes blank image template
    imageBottom = ones(size(images{1},1),size(images{1},2));
    for n = 2:length(indices)
        %Adds together all images 2:end
        imageBottom = imageBottom & images{indices(n)};
    end
    %generates stats from overlapping images
    imageContact = bwconncomp(imdilate(images{indices(1)},strel('disk',1))...
        & imageBottom,4);
    %number of dilated objects
    baseNumber = max(max(bwlabel(imdilate(images{indices(1)},strel('disk',1)),4)));
    %number of contact objects
    number = imageContact.NumObjects;
    
    %creates map of dilated objects
    cc = bwlabel(imdilate(images{indices(1)},strel('disk',1)),4);
    %list of pixels in each object in the overlap image
    aaPix = [imageContact.PixelIdxList];
    nFrac = 0;
    aFrac = 0;    
    if ~isempty(aaPix)
        aaIdxList = 0;
        for z = 1:length(aaPix)
            %takes one pixel from each object in overlap image
            aaIdxList(z) = aaPix{z}(1);
        end
        %cross-checks pixels to identify which object they correspond to
        bb = cc(aaIdxList);
        bb = unique(bb);
        %number fraction calculated based on number of objects
        %post-dilation.
        nFrac = double(length(bb)/baseNumber);
        %area fraction calculated from non-dilation overlap
        aFrac = bwarea(~images{indices(1)} & imageBottom)/baseArea;
    end
end
%PERMUTATION FUNCTION
function [Matrix,Index] = npermutek(N,K)
%NPERMUTEK Permutation of elements with replacement/repetition.
% MAT = NPERMUTEK(N,K) returns all possible samplings of length K from   
% vector N of type: ordered sample with replacement.  
% MAT has size (length(N)^K)-by-K, where K must be a scalar.
% [MAT, IDX] = NPERMUTEK(N,K) also returns IDX such that MAT = N(IDX).
% N may be of class: single, double, or char.  If N is single or double,
% both MAT and IDX will be of the same class.
%
% For N = 1:M, for some integer M>1, all(MAT(:)==IDX(:)), so there is no 
% benefit to calling NPERMUTEK with two output arguments.
%
% Examples:
%         MAT = npermutek([2 4 5],2)
%
%  MAT =
% 
%       2     2
%       2     4
%       2     5
%       4     2
%       4     4
%       4     5
%       5     2
%       5     4
%       5     5
%
% NPERMUTEK also works on characters.
%
%         MAT = npermutek(['a' 'b' 'c'],2)
%  MAT =
%        
%      aa
%      ab
%      ac
%      ba
%      bb
%      bc
%      ca
%      cb
%      cc
%
% See also perms, nchoosek
%
% Also on the web:
% http://mathworld.wolfram.com/BallPicking.html
% See the section on Enumerative combinatorics below: 
% http://en.wikipedia.org/wiki/Permutations_and_combinations
% Author:  Matt Fig
% Contact:  popkenai@yahoo.com
    if nargin ~= 2
        error('NPERMUTEK requires two arguments. See help.')
    end
    if isempty(N) || K == 0,
       Matrix = [];  
       Index = Matrix;
       return
    elseif floor(K) ~= K || K<0 || ~isreal(K) || numel(K)~=1 
        error('Second argument should be a real positive integer. See help.')
    end
    LN = numel(N);  % Used in calculating the Matrix and Index.
    if K==1
        Matrix = N(:); % This one is easy to calculate.
        Index = (1:LN).';
        return
    elseif LN==1
        Index = ones(K,1);
        Matrix = N(1,Index);
        return  
    end
    CLS = class(N);
    if ischar(N)
        CLS = 'double';  % We will deal with this at the end.
        flg = 1;
        N = double(N);
    else
        flg = 0;
    end
    L = LN^K;  % This is the number of rows the outputs will have.
    Matrix = zeros(L,K,CLS);  % Preallocation.
    D = diff(N(1:LN));  % Use this for cumsumming later.
    LD = length(D);  % See comment on LN. 
    VL = [-sum(D) D].';  % These values will be put into Matrix.
    % Now start building the matrix.
    TMP = VL(:,ones(L/LN,1,CLS));  % Instead of repmatting.
    Matrix(:,K) = TMP(:);  % We don't need to do two these in loop.
    Matrix(1:LN^(K-1):L,1) = VL;  % The first column is the simplest.
    if nargout==1
        % Here we only have to build Matrix the rest of the way.
        for ii = 2:K-1
            ROWS = 1:LN^(ii-1):L;  % Indices into the rows for this col.
            TMP = VL(:,ones(length(ROWS)/(LD+1),1,CLS));  % Match dimension.
            Matrix(ROWS,K-ii+1) = TMP(:);  % Build it up, insert values.
        end
    else
        % Here we have to finish Matrix and build Index.
        Index = zeros(L,K,CLS);  % Preallocation.
        VL2 = ones(size(VL),CLS);  % Follow the logic in VL above.
        VL2(1) = 1-LN;  % These are the drops for cumsum.
        TMP2 = VL2(:,ones(L/LN,1,CLS));  % Instead of repmatting.
        Index(:,K) = TMP2(:);  % We don't need to do two these in loop.
        Index(1:LN^(K-1):L,1) = 1;  
        for ii = 2:K-1
            ROWS = 1:LN^(ii-1):L;  % Indices into the rows for this col.
            F = ones(length(ROWS)/(LD+1),1,CLS);  % Don't do it twice!
            TMP = VL(:,F);  % Match dimensions.
            TMP2 = VL2(:,F);
            Matrix(ROWS,K-ii+1) = TMP(:); % Build them up, insert values.
            Index(ROWS,K-ii+1) = TMP2(:);  
        end

        Index(1,:) = 1;  % The first row must be 1 for proper cumsumming.
        Index = cumsum(Index);  % This is the time hog.
    end
    Matrix(1,:) = N(1);  % For proper cumsumming.
    Matrix = cumsum(Matrix);  % This is the time hog.
    if flg
        Matrix = char(Matrix);  % char was implicitly cast to double above.
    end
end
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





