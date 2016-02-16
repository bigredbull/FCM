clear all
clc;
%% ��������
load wine.txt;
cluster_n=6;

data = wine;
data(:,1) = [];
data(:,size(data,2)) = [];
data_n = size(data, 1); %���ݵĸ��� 
in_n = size(data, 2);% ����ά��

%% �������
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

tic
%% ��ʼ�������Ⱦ��󲢹�һ
U = rand(cluster_n, data_n); %rand()����������� 
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);%��һ��

%% ��ʼ����
for i = 1:max_iter,%������������
	tic
    mf = U.^expo;      % MF matrix after exponential modification
    center = mf*data./((ones(size(data, 2), 1)*sum((mf')))'); % �����µľ�������
   
       out = zeros(size(center, 1), size(data, 1));  %ÿ���㵽ÿ�����ĵľ��룬����Ϊ������
     
        if size(center, 2) > 1,%������ά������һִ�����³���
            for k = 1:size(center, 1),%��K��ֵ
                abc = ((data-ones(size(data,1),1) * center(k,:)).^2)';
 	             out(k, :) = sqrt(sum(abc));%�õ�ŷ�Ͼ���
            end
       else	% 1-D data
            for k = 1:size(center, 1),
 	             out(k, :) = abs(center(k)-data)';
            end
        end
   
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
toc

plot(U(1,:),'-ro');
grid on
hold on
plot(U(2,:),'-g*');
plot(U(3,:),'-b+');

ylabel('Membership degrees')
legend('FCM1','FCM2','FCM3','location','northeast');

toc