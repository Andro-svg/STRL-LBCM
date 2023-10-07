function img = reconstruct_image(img,ten,tenPos,patchSize,patchNum,idx)

% img = zeros([imgHei, imgWid]);
tenImage = ten(:,:,(idx-1)*patchNum+1 : idx*patchNum, :);
tenPosition = tenPos(:,:,(idx-1)*patchNum+1 : idx*patchNum, :);
[~,~,~,n4]=size(tenImage);

for i = 1:n4
    for j = 1:patchNum
        position = tenPosition(:,:,j,i);
        row = position(1);
        col = position(2);
        img(row: row + patchSize-1, col:col + patchSize-1) = tenImage(:,:,j,i);
    end
end