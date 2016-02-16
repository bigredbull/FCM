
clear all
clc;

a=imread('MRI.jpg');
I=imnoise(a,'salt & pepper',0.03);
figure(1);           
 imshow(I);title('加噪图像');

[height,width,c]=size(a);
if c~=1
    a=rgb2gray(a);
end

a=double(a);
[row,column]=size(a); 
data = a(:);
data_n = size(data,1);

beta = 80;
cluster_num = 4;%将图像分为四类
default_options = [2.0;	% exponent for the partition matrix U
		300;	% max. number of iteration
		0.01;	% min. amount of improvement
		1];	% info display during iteration 
    
options = default_options;

m = options(1);		% p，图像的质量参数
iter_num = options(2);		% Max. iteration 最大迭代次数
threshold = options(3);		% Min. improvement  最小进化步长
display = options(4);		% Display info or not 显示信息与否
costfunction = zeros(iter_num, 1);

tic
% 初始化隶属度并归一
membership = zeros(height,width,cluster_num);
for i=1:height 
    for j=1:width
        member_sum=0;
        for k=1:cluster_num
            membership(i,j,k)=rand();
            member_sum=member_sum + membership(i,j,k);
        end        
        for p =1:cluster_num
            membership(i,j,p)=membership(i,j,p)/member_sum;
        end
    end
end

for o = 1:iter_num %迭代次数控制
    
    %计算初始中心
    center = zeros(cluster_num,1);
    for u = 1:cluster_num
        sum = 0;
        sum1 = 0;
             for h=1:height 
                for t=1:width
                 sum = sum + (membership(h,t,u).^m) * a(h,t);
                 sum1 = sum1 + membership(h,t,u).^m;
                 end
             end
             center(u) = sum/sum1;
    end         
    
   for i=1:height 
                    for j =1:width
                        up = i-1;
                        down = i+1;
                        left = j-1;
                        right = j+1;
                        if( up < 1)
                            up = 1;
                        end 
                        if( down > height)
                            down = height;
                        end 
                        if( left< 1)
                            left = 1;
                        end
                        if( right > width)
                            right = width;
                        end                        
                        s = 0;
                         for x = up : down
                              for y = left : right
                                    for u = 1:cluster_num
                                        for r = 1:cluster_num
                                            s = s + membership(x,y,r).^m;
                                        end
                                        s = s - membership(x,y,u).^m;
                                    end
                              end
                         end
                    end
   end             
               
    for h = 1:height
        for t = 1:width
            for k = 1:cluster_num
                    costfunction(o) = costfunction(o) + membership(h,t,k).^m*(a(h,t) - center(k))^2 + (beta/2) * membership(h,t,k)^m*s;
                    tmp = ((a(h,t)-center(k))^2  + beta * s).^(-1/(m - 1));
                    tomp = 0;
                    for p = 1:cluster_num
                        tomp = tomp + (a(h,t)-center(p))^2;
                        tp = (tomp + beta * s) .^(-1/(m - 1));
                    end
                    membership(h,t,k) = tmp./tp;
            end
        end
    end
    
    %%%%%%归一化隶属度
    for i=1:height 
        for j = 1:width
            member_sum = 0;
            for k = 1:cluster_num
                member_sum = member_sum + membership(i,j,k);
            end
            for p = 1:cluster_num
                membership(i,j,p) = membership(i,j,p) / member_sum;
            end
        end
    end    
    
    if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', o, costfunction(o));
        %输出迭代次数和函数的结果
	end
	% check termination condition
	if o > 1,  %进化步长控制
		if abs(costfunction(o) - costfunction(o-1)) < threshold, break; end,
    end
end

 toc   

  A = ones(height,width,1);
for i = 1:height
    for j = 1:width
        if (fix(a(i,j) / 85) == 1)
            A(i,j,1) = 2;
        end
        if (fix(a(i,j) / 85) == 2)
            A(i,j,1) = 3;
        end
        if (fix(a(i,j,1) / 85) == 3)
            A(i,j,1) = 4;
        end
    end
end    
A = reshape(A,1,data_n);

  newing = zeros(row,column);
for i=1:row
    for j=1:column         
        maxmembership=membership(i,j,1);
        index=1;
        for k=2:cluster_num            
            if(membership(i,j,k)>maxmembership)
                maxmembership=membership(i,j,k);
                index=k;
            end
        end
        newing(i,j) = round(255 * (1-(index-1)/(cluster_num-1)));
    end 
end

B = reshape(newing,1,data_n);
b = fix((max(B) - B(1,1)) / cluster_num);
for i = 2:data_n
    if B(1,i) == B(1,1)
        B(1,i) = 1;
    elseif (fix(B(1,i) / b) == 2)
        B(1,i) = 2;
    elseif (fix(B(1,i) / b)  == 3)
        B(1,i) = 3;
    else
        B(1,i) = 4;
    end
end
B(1,1) = 1;

sum = 0;
for i = 1:data_n
    if ( A(1,i) ~= B(1,i))
        sum = sum + 1;
    end
end
MCR = sum / data_n;
fprintf('MCR = %d\n',MCR);

S = 0;
for i = 1:data_n
    S = S + (A(1,i) - B(1,i)).^2;
end

RMS = sqrt(S / (data_n * (data_n - 1)));
fprintf('RMS = %d\n',RMS);

figure(2);
imshow((uint8(newing)));
title('RFCM分割后的图像');
