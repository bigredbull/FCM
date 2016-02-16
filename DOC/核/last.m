%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc;

 a=imread('picture.bmp');
I=imnoise(a,'salt & pepper',0.03);
figure(2);           
 imshow(I);title('����ͼ��');
 
[height,width,c]=size(a);
if c~=1
    a=rgb2gray(a);
end

a=double(a);
[row,column]=size(a); 

% ��������
beta = 1.5;
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
%%%%%%%%%%%%%%%%%%%%%%%

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
K = zeros(height,width,cluster_num);
center = zeros(cluster_num,1);

for i = 1:cluster_num
    for j = 1:height
        for t = 1:width
            K(j,t,i) = exp((-1) * (a(j,t)-center(i))^2 / (beta.^2));
        end
    end
end

for o = 1:iter_num %������������
    
    %�����ʼ����
    
    for u = 1:cluster_num
        sum = 0;
        sum1 = 0;
             for h=1:height 
                for t=1:width
                 sum =  sum + membership(h,t,u).^m * K(h,t,u) * a(h,t);
                 sum1 = sum1 + membership(h,t,u).^m  * K(h,t,u) / (beta.^2);
                 end
             end
             center(u) = sum / sum1;
    end   
    
    for i = 1:cluster_num
        for j = 1:height
            for t = 1:width
                K(j,t,i) = exp((-1) * (a(j,t)-center(i))^2 / (beta.^2));
            end
        end
    end
               
    for h = 1:height
        for t = 1:width
            for k = 1:cluster_num 
                costfunction(o) = costfunction(o) + 2 * membership(h,t,k).^m * (1 - K(h,t,k));
                temp =(1 - K(h,t,k)).^(-1/(m - 1));
                top = 0;
                for p = 1:cluster_num
                    top = top + (1 - K(h,t,k)).^(-1/(m - 1));
                end
                membership(h,t,k) = temp./top;
            end
        end
    end  
       
    if display,
        fprintf('iter.count = %d, obj. fcn = %f\n',o, costfunction(o));
    end
    if o >1
        if abs(costfunction(o)-costfunction(o-1)) < threshold , break;end;
    end

end
 toc   
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
        newimg(i,j) = round(255 * (1-(index-1)/(cluster_num-1)));
    end
end
figure(5);
imshow(uint8(newimg));title('�ָ���ͼ��');       