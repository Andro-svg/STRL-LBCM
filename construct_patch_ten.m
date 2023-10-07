function [patchTen, patchNumber, patchPosition] = construct_patch_ten(img, patchSize, slideStep)
if ~exist('patchSize', 'var')
    patchSize = 50;
end

if ~exist('slideStep', 'var')
    slideStep = 10;
end

% img = reshape(1:9, [3 3])
% img = reshape(1:12, [3 4])
% patchSize = 2;
% slideStep = 1;
[imgHei, imgWid] = size(img);

rowPatchNum = ceil((imgHei - patchSize) / slideStep) + 1;
colPatchNum = ceil((imgWid - patchSize) / slideStep) + 1;
rowPosArr = [1 : slideStep : (rowPatchNum - 1) * slideStep, imgHei - patchSize + 1];
colPosArr = [1 : slideStep : (colPatchNum - 1) * slideStep, imgWid - patchSize + 1];

%% arrayfun version, identical to the following for-loop version
% [meshCols, meshRows] = meshgrid(colPosArr, rowPosArr);
% idx_fun = @(row,col) img(row : row + patchSize - 1, col : col + patchSize - 1);
% patchCell = arrayfun(idx_fun, meshRows, meshCols, 'UniformOutput', false);
% patchTen = cat(3, patchCell{:});

%% for-loop version
% 常规的排队法：
%  1  2 
% 3  4  5
%  6  7  8  9  10
% patchTen = zeros(patchSize, patchSize, rowPatchNum * colPatchNum);
% patchPosition = zeros(1,2,rowPatchNum*colPatchNum);
% k = 0;
% for col = colPosArr
%     for row = rowPosArr
%         k = k + 1;
%         tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1);
%         patchTen(:, :, k) = tmp_patch;
%         patchPosition(:,:,k) = [row , col];
%     end
% end
% patchNumber=k;

% s型
%  1  2  3  4  5
%  10  9  8  7  6
patchTen = zeros(patchSize, patchSize, rowPatchNum * colPatchNum);
patchPosition = zeros(1,2,rowPatchNum*colPatchNum);
k = 0;
for row = rowPosArr
    if mod(find(rowPosArr==row),2)==0     % 偶数行
        colPosArr1 = colPosArr;
        for col = colPosArr1
            k = k + 1;
            tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1);
            patchTen(:, :, k) = tmp_patch;
            patchPosition(:,:,k) = [row , col];
        end
    else
        colPosArr2 = fliplr(colPosArr);
        for col = colPosArr2
            k = k + 1;
            tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1);
            patchTen(:, :, k) = tmp_patch;
            patchPosition(:,:,k) = [row , col];
        end
    end
end
patchNumber=k;
