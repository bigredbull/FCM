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

cluster_num = 4;
default_options = [10;	% �����Ⱥ������ݴη�
		300;	% ����������
		1e-5;	%����
		1];	% �ж�����

options = default_options;


expo = options(1);		% Exponent for U �����Ⱥ������ݴη�
max_iter = options(2);		% Max. iteration ����������
min_impro = options(3);		% Min. improvement  ��С��������
display = options(4);		% Display info or not ��ʾ��Ϣ���

obj_fcn = zeros(max_iter, 1);	% Array for objective function
membership = zeros(height,width,cluster_num);
center = zeros(cluster_num,1);

tic
% ��ʼ�������Ȳ���һ
for i=1:height 
    for j=1:width
        member_sum=0;
        for k=1:cluster_num
            membership(i,j,k)=rand();
            member_sum = member_sum + membership(i,j,k);
        end        
        for p =1:cluster_num
            membership(i,j,p) = membership(i,j,p) / member_sum;
        end
    end
end

tic
for i = 1:max_iter,%������������
    mf = membership.^expo;     
    %%%%%%%%������������
    for m = 1:cluster_num
        to = 0;
        tp =0;
        for j = 1:height
            for t = 1:width
                to = to + membership(j,t,m) * a(j,t);
                tp = tp + membership(j,t,m);
            end
        end
        center(m,1) = to / tp; 
    end
    
    %%%%%%%�õ�ŷʽ�����Լ�Ŀ�꺯��
    out = zeros(height,width,cluster_num);
    for m =1:height
        for j =1:width
            for t = 1:cluster_num
                out(m,j,t) = abs(a(m,j) - center(t,1));
                obj_fcn(i) = obj_fcn(i) + (membership(m,j,t).^expo) * (out(m,j,t).^2);
            end
        end
    end
    
    for m = 1:height
        for j = 1:width
            for r = 1:cluster_num
                top =0;
                for t = 1:cluster_num
                    top = top + (out(m,j,r) / out(m,j,t)).^(expo - 1);
                end
                membership(m,j,r) = 1 / top;
            end
        end
    end
    
    %%%%%%��һ��������
    for m=1:height 
        for j = 1:width
            member_sum = 0;
            for k = 1:cluster_num
                member_sum = member_sum + membership(m,j,k);
            end
            for p = 1:cluster_num
                membership(m,j,p) = membership(m,j,p) / member_sum;
            end
        end
    end
          
     if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
        %������������ͺ����Ľ��
	end
	% check termination condition
	if i > 1,  %������������
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
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
title('FCM�ָ���ͼ��');    