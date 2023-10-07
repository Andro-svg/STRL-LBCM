function [patchTen, patchPosition] = construct_local_patch_ten(img,rowPosArr,colPosArr, patchSize,N,M)

% patchTen = zeros(patchSize, patchSize, N * M);
% patchPosition = zeros(1,2,N * M);
% k = 0;
% for col = colPosArr
%     for row = rowPosArr
%         k = k + 1;
%         tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1);
%         patchTen(:, :, k) = tmp_patch;
%         patchPosition(:,:,k) = [row , col];
%     end
% end


[meshCols, meshRows] = meshgrid(colPosArr, rowPosArr);
idx_fun1 = @(row,col) img(row : row + patchSize - 1, col : col + patchSize - 1);

patchCell = arrayfun(idx_fun1, meshRows, meshCols, 'UniformOutput', false);
patchTen = cat(3, patchCell{:}); 

idx_fun2 = @(row,col) [row,col];
patchPosCell = arrayfun(idx_fun2, meshRows, meshCols, 'UniformOutput', false);
patchPosition = cat(3, patchPosCell{:});


end