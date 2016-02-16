clear all%清除workspace中的所有变量
clc;
 I=imread('Dfour244true.bmp');%读入图像文件
%I=imread('corner.bmp');

% [counts,x]=imhist(I);%灰度统计
% figure(1);
% imhist(I);
[m,n,h]=size(I);%测图像的大小
if h~=1
    I=rgb2gray(I);
end
%figure(1);
%imshow(I);title('原图像');
 I=imnoise(I,'salt & pepper',0.03);
% I=imnoise(I,'gaussian',0,0.03);
 figure(2);           
 imshow(I);title('加噪图像'); %显示原来的图像
% imwrite(I,'0.03_gaussian_four.bmp');

U9=cell(m,n);alpha=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=4;          %聚类个数
e=0.01;
V1=zeros(1,c);     %V1是旧的聚类中心
V2=zeros(1,c);     %V2是新的聚类中心
U1=cell(m,n);
%U9=cell(m,n);
m1=2;              %确定加权指数m
I1=I;
[U1,V1]=Initialization(I1,c);
count=0;         %循环次数
flag=1;
w=2.0;
I1=double(I1);     %转换变量的类型
%%%%%%%%%%%%%%%%%%%%%%%%%初始化完成%%%%%%%%%%%%%%%%%%%%%%
[C]=calculateC(I1,w);%基于方差均值比
%%%%%%%%%%%%%%%%%%%%%%%%%%更新%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;%tic1 
t1=clock; 

while (flag==1&&count<500)
    tic ;%tic2 
    t2=clock; 
 %   pause(3*rand) 
    
    [D1]=Distance(V1,c,I1);
    [G1,U1]=calculateG(U1,D1,I1,m1,c,C);%计算G矩阵和隶属度矩阵（c*(m*n）)
    V2=center(U1,c,m1,I1,D1);
    count=count+1;
     for i=1:m
         for j=1:n
             for k=1:c
                 costfunction(count)=U1{i,j}(k,1)^m1*(1-D1{i,j}(k,1))+G1{i,j}(k,1)*U1{i,j}(k,1)^m1;
             end
         end
     end
     
%%%%%%%%%%%%%抑制式修改
%        for i=1:m
%          for j=1:n
%              for k=1:c    
%                 U9{i,j}(:,1)=U1{i,j}(:,1).*alpha;
%          
%                 [~,max_loca]=max(U1{i,j}(:,1));
%                 tempU=alpha.*U1{i,j}(max_loca,1)+(1-alpha);
%                 U9{i,j}(max_loca,1)=tempU;
%           
%                 
%              end
%          end
%        end
%        U1=U9;
%%%%%%%%%%%%%%%%%%%
     
     
    fprintf('iter.count = %d, obj. fcn = %f\n',count, costfunction(count));
    if max(abs(V2-V1))<e;%迭代停止条件
        flag=0;
    else
        V1=V2;
    end
   %disp(sprintf('Iterations = %d',count));
   %计算每次循环的时间 
    disp(['etime计算第',num2str(count),'次循环运行时间：',num2str(etime(clock,t2))]); 
end
%计算运行的总时间
disp(['etime程序总运行时间：',num2str(etime(clock,t1))]);

[I2,I3]=defuzzy(U1,I,V2);%去模糊化，对原始图像进行标记

%%%%%%%%%%%%%%%%%%%%%显示分割结果%%%%%%%%%%%%%%%%%%%%%
figure(3);
I3=uint8(I3);   % 转换为无符号型整数 8表示8位二进制整数  范围0~255
imshow(I3);title('分割后的图像');   %显示分割后的图像
% imwrite(I3,'result_four.bmp');