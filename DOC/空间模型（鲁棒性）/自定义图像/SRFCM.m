
clear all
clc;

image = zeros(96,96); 

for i = 1:58
    for j = 49:96
        image(i,j) = 75;
    end
end

for i = 59:77
    for j = 1:96
        image(i,j) = 50;
    end
end

for i = 78:96
    for j = 1:48
        image(i,j) = 50;
    end
end

for i = 78:96
    for j = 49:96
        image(i,j) = 25;
    end
end

a = image;

noise = uint8(image);
noise = imnoise(noise,'gaussian',100 / (255 * 255));
figure(1);
imagesc(noise);
colormap(gray);
title('����ͼ��');

[height,width,c]=size(a);
if c~=1
    a=rgb2gray(a);
end

a=double(a);
[row,column]=size(a); 
data = a(:);
data_n = size(data,1);

beta = 435;
alpha=0.5;
cluster_num = 4;%��ͼ���Ϊ����
default_options = [2.0;	% exponent for the partition matrix U
		300;	% max. number of iteration
		0.01;	% min. amount of improvement
		1];	% info display during iteration 
    
options = default_options;

m = options(1);		% p��ͼ�����������
iter_num = options(2);		% Max. iteration ����������
threshold = options(3);		% Min. improvement  ��С��������
display = options(4);		% Display info or not ��ʾ��Ϣ���
costfunction = zeros(iter_num, 1);

tic
% ��ʼ�������Ȳ���һ
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

for o = 1:iter_num %������������
    
    %�����ʼ����
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
             center(u) = sum / sum1;
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
      
   %%%%%%%����������
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
    
    %%%%%%��һ��������
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
    
    %%%%%%%����ʽ�޸�
    ha = 0;
    for i = 1:height
        for j = 1:width
            for p = 1:cluster_num
            ha = ha + membership(i,j,p).^2;
            end
        end
    end 
  
    alpha = ( cluster_num/(cluster_num - 1)) * (1 -(1/data_n) * ha);
    
     membership1 = membership.*alpha; 
    
    for i = 1:height
        for j = 1:width
             [max_data,max_local] = max(membership(i,j,:));
             temembership = alpha.* membership(i,j,max_local) + (1 - alpha);
             membership1(i,j,max_local) = temembership;
        end
    end
    membership = membership1 * alpha;
    
     %%%%%%��һ��������
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
        %������������ͺ����Ľ��
	end
	% check termination condition
	if o > 1,  %������������
		if abs(costfunction(o) - costfunction(o-1)) < threshold, break; end,
    end
end

 toc   
 
   A = ones(height,width,1);
for i = 1:height
    for j = 1:width
        if (a(i,j) == 75)
            A(i,j,1) = 2;
        end
        if (a(i,j) == 50)
            A(i,j,1) = 3;
        end
        if (a(i,j,1) == 25)
            A(i,j,1) = 4;
        end
    end
end    
A = reshape(A,1,data_n);

  B = ones(row,column,1);
for i=1:row
    for j=1:column       
        maxmembership=membership(i,j,1);
        for k=2:cluster_num            
            if(membership(i,j,k)>maxmembership)
                maxmembership = membership(i,j,k);
                B(i,j,1)=k;
            end
        end
    end
end

B = reshape(B,1,data_n);

t = MCR(A,B);
  %���ͼ��
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
figure(2);
imagesc(uint8(newing));
colormap(gray);
title('SRFCM�ָ���ͼ��');       
           
    
  