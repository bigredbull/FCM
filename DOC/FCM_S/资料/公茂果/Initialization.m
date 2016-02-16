function [U2,V2]=Initialization(I,c)

[counts,x]=imhist(I); %灰度频率
 figure(5);
 imhist(I);
[m,n]=size(I);%测矩阵的大小
U2=cell(m,n);
for i=1:m
    for j=1:n
        U2{i,j}=zeros();
    end
end
e=0.001;
V1=zeros(1,c);     %V1是旧的聚类中心
V2=zeros(1,c);     %V2是新的聚类中心
U1=zeros(c,256);
m1=2;              %确定加权指数m
num=0;             %初始迭代次数
V1(1)=50;          %初始化聚类中心V
V1(2)=150;
V1(3)=200;
V1(4)=255;
I=double(I);%转换变量的类型
flag=1;
%%%%%% 开始迭代 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (flag==1&&num<500) 
%%%%%%%%%%%%%%% 更新隶属度矩阵 %%%%%%%%%%%%%%%
 for j=1:256 
       for i=1:c
           sum=0;
             if(x(j)-V1(i)==0)
                  U1(i,j)=1;
                  for k=1:i-1
                      U1(k,j)=0;
                  end
                  for k=i+1:c
                      U1(k,j)=0; 
                  end
                  break;
              else
                  for k=1:c
                      sum=sum+(x(j)-V1(k))^(2/(m1-1));
                  end
                  U1(i,j)=sum/(x(j)-V1(i))^(2/(m1-1));
               end
       end
   end
%%%%%%%%%%%%% 归一化 %%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:256
    temp=0;
    for i=1:c
        temp=temp+U1(i,j);
    end
    for i=1:c
         U1(i,j)=U1(i,j)/temp;
    end
end
%%%%%%%%%%%%%%%更新聚类中心%%%%%%%%%%%%%%%%%%%%%%
for i=1:c
    sum=0;sum1=0;
    for j=1:256
        sum=sum+counts(j)*U1(i,j)^m1;
        sum1=sum1+x(j)*counts(j)*U1(i,j)^m1;  
    end
    V2(i)=sum1/sum;
end
%%%%%%%%%%%%%%% 终止判断 %%%%%%%%%%%%%%%%%%%%%%
num=num+1;
    if max(abs(V2-V1))<e;%迭代停止条件
        flag=0;
    else
        V1=V2;
    end
    disp(sprintf('Iterations = %d',num));
end
for i=1:m
    for j=1:n
        for k=1:256
            if x(k)==I(i,j)
                U2{i,j}=U1(:,k);
            end
        end         
    end
end