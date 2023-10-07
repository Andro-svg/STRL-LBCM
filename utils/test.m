clc
clear
a(:,:,1)=[1,2;3,4];
a(:,:,2)=[3,9;7,6];
a(:,:,3)=[5,2;1,0];
a 

b = my_Unfold(a,size(a),2)

c = my_Fold(b,size(a),2)