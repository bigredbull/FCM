function [I2,I3]=defuzzy(U1,I,V2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入：隶属度矩阵U1，原始图像I
%输出：经过标记的矩阵I2和分割后的图像I3
%功能：计算标记矩阵和分割后的图像
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(I);
I2=zeros(m,n);
for i=1:m
    for j=1:n
        [C,index]=max(U1{i,j});
        I2(i,j)=index;
        I3(i,j)=V2(index);
    end
end