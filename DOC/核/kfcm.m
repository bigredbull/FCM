%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc;

 load data.txt;
 height = size(data,1);
 width = size(data,2);


% ��������
beta = 0.005;
cluster_num = 3;%��ͼ���Ϊ����
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
% membership1 = zeros(height,width,cluster_num);
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

center = zeros(cluster_num,1);
for o = 1:iter_num %������������
    
    %�����ʼ����
    
    for u = 1:cluster_num
        sum = 0;
        sum1 = 0;
             for h=1:height 
                for t=1:width
                 sum =  sum + membership(h,t,u).^m * exp((-1) * (data(h,t) - center(u).^2) / (beta.^2)) * data(h,t);
                 sum1 = sum1 + membership(h,t,u).^m  * exp((-1) * ((data(h,t) - center(u))^2) / (beta.^2));
                 end
             end
             center(u) = sum / sum1;
    end    
               
    for h = 1:height
        for t = 1:width
            for k = 1:cluster_num 
                costfunction(o) = costfunction(o) + 2 * membership(h,t,k).^m * (1 - exp( -1 * ((data(h,t) - center(k)).^2)/(beta.^2)));
                temp =(1 - exp( -1 * ((data(h,t) - center(k)).^2)/(beta.^2))).^(-1 / (m -1));
                top = 0;
                for p = 1:cluster_num
                    top = top + (1 - exp((-1 * (data(h,t) - center(p)).^2) /(beta.^2))).^(-1 / (m -1));
                end
                membership(h,t,k) = temp / top;
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
%  iter_n = o;	% ʵ�ʵĵ������� 
% %obj_fcn(iter_n+1:iter_num) = [];
% plot(membership(:,:,1),'-ro');
% grid on
% hold on
% plot(membership(:,:,2),'-g*');
% plot(membership(:,:,3),'-b+');
% 
% ylabel('Membership degrees')
% legend('FCM1','FCM2','FCM3','location','northeast');
% 
% toc  

plot(membership(1,:),'-ro');
grid on
hold on
plot(membership(2,:),'-g*');
plot(membership(3,:),'-b+');

ylabel('Membership degrees')
legend('FCM1','FCM2','FCM3','location','northeast');