clear all 
clc;

a=imread('MRI.jpg');
I=imnoise(a,'salt & pepper',0.05);
figure(1);           
 imshow(I);title('����ͼ��');

[height,width,c]=size(a);
if c~=1
    a=rgb2gray(a);
end

a=double(a);
[row,column]=size(a); 
data = a(:);
data_n = size(data,1);
       
cluster_num = 4;
default_options = [2.0;	% �����Ⱥ������ݴη�
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

%%%%%%֤�����Զ���ͼ���е�MCR���ܼ��㣬���ڴ˼�������ֱ����newing��A��Ƚ�

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

RMS = sqrt(S / (data_n * (data_n -1)));
fprintf('RMS = %d\n',RMS);

figure(2);
imshow((uint8(newing)));
title('FCM�ָ���ͼ��');   
