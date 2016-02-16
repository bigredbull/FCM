clear all
clc;

load glass.txt;
cluster_n=6;

data = glass;
data(:,1) = [];
data(:,size(data,2)) = [];
data_n = size(data, 1); %���ݵĸ���    size(A,1)����A��������size(A,2)����A������
in_n = size(data, 2);% ����ά��

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
%��ʼ�������Ⱦ���
U = rand(cluster_n, data_n); %rand()����������� 
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);%��һ��

for i = 1:max_iter,%������������
	
    tic
    mf = U.^expo;      % MF matrix after exponential modification
    center = mf*data./((ones(size(data, 2), 1)*sum((mf')))'); % �����µľ�������
   
       out = zeros(size(center, 1), size(data, 1));  %ÿ���㵽ÿ�����ĵľ��룬����Ϊ������
     
        if size(center, 2) > 1,%������ά������һִ�����³���
            for k = 1:size(center, 1),%��K��ֵ
 	             out(k, :) = sqrt(sum((((data-ones(size(data, 1), 1)*center(k, :)).^2)')));%�õ�ŷ�Ͼ���
            end
       else	% 1-D data
            for k = 1:size(center, 1),
 	             out(k, :) = abs(center(k)-data)';
            end
        end
    % dist =out;
   %%%%%%%%%%%%%%%%%%%%
    obj_fcn(i) = sum(sum((out.^2).*mf));  % Ŀ�꺯��
    tmp = out.^(-2/(expo-1));      % �����µ������Ⱦ�����⹫ʽ���
    U= tmp./(ones(cluster_n, 1)*sum(tmp));  % �µ������Ⱦ���

	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
        %������������ͺ����Ľ��
	end
	% check termination condition
	if i > 1,  %������������
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
    end
    toc
end

iter_n = i;	% ʵ�ʵĵ������� 
obj_fcn(iter_n+1:max_iter) = [];

B = ones(data_n,1);
for i = 1:size(U,1)
    for j = 1:size(U,2)
        if U(B(j,1),j)< U(i,j)
            B(j,1) = i;
        end
    end
end

t = size(glass,2);
A = glass(:,t);

B = reshape(B,1,data_n);
A = reshape(A,1,data_n);
t = nmi(A,B);

plot(U(1,:),'-ro');
grid on
hold on
plot(U(2,:),'-g*');
plot(U(3,:),'-b+');
plot(U(4,:),'-m*');
plot(U(5,:),'-y+');
plot(U(6,:),'-k*');


ylabel('Membership degrees')
legend('FCM1','FCM2','FCM3','FCM4','FCM5','FCM6','location','northeast');

toc