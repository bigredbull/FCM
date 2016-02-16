clear all%���workspace�е����б���
clc;
 I=imread('Dfour244true.bmp');%����ͼ���ļ�
%I=imread('corner.bmp');

% [counts,x]=imhist(I);%�Ҷ�ͳ��
% figure(1);
% imhist(I);
[m,n,h]=size(I);%��ͼ��Ĵ�С
if h~=1
    I=rgb2gray(I);
end
%figure(1);
%imshow(I);title('ԭͼ��');
 I=imnoise(I,'salt & pepper',0.03);
% I=imnoise(I,'gaussian',0,0.03);
 figure(2);           
 imshow(I);title('����ͼ��'); %��ʾԭ����ͼ��
% imwrite(I,'0.03_gaussian_four.bmp');

U9=cell(m,n);alpha=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=4;          %�������
e=0.01;
V1=zeros(1,c);     %V1�Ǿɵľ�������
V2=zeros(1,c);     %V2���µľ�������
U1=cell(m,n);
%U9=cell(m,n);
m1=2;              %ȷ����Ȩָ��m
I1=I;
[U1,V1]=Initialization(I1,c);
count=0;         %ѭ������
flag=1;
w=2.0;
I1=double(I1);     %ת������������
%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ�����%%%%%%%%%%%%%%%%%%%%%%
[C]=calculateC(I1,w);%���ڷ����ֵ��
%%%%%%%%%%%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;%tic1 
t1=clock; 

while (flag==1&&count<500)
    tic ;%tic2 
    t2=clock; 
 %   pause(3*rand) 
    
    [D1]=Distance(V1,c,I1);
    [G1,U1]=calculateG(U1,D1,I1,m1,c,C);%����G����������Ⱦ���c*(m*n��)
    V2=center(U1,c,m1,I1,D1);
    count=count+1;
     for i=1:m
         for j=1:n
             for k=1:c
                 costfunction(count)=U1{i,j}(k,1)^m1*(1-D1{i,j}(k,1))+G1{i,j}(k,1)*U1{i,j}(k,1)^m1;
             end
         end
     end
     
%%%%%%%%%%%%%����ʽ�޸�
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
    if max(abs(V2-V1))<e;%����ֹͣ����
        flag=0;
    else
        V1=V2;
    end
   %disp(sprintf('Iterations = %d',count));
   %����ÿ��ѭ����ʱ�� 
    disp(['etime�����',num2str(count),'��ѭ������ʱ�䣺',num2str(etime(clock,t2))]); 
end
%�������е���ʱ��
disp(['etime����������ʱ�䣺',num2str(etime(clock,t1))]);

[I2,I3]=defuzzy(U1,I,V2);%ȥģ��������ԭʼͼ����б��

%%%%%%%%%%%%%%%%%%%%%��ʾ�ָ���%%%%%%%%%%%%%%%%%%%%%
figure(3);
I3=uint8(I3);   % ת��Ϊ�޷��������� 8��ʾ8λ����������  ��Χ0~255
imshow(I3);title('�ָ���ͼ��');   %��ʾ�ָ���ͼ��
% imwrite(I3,'result_four.bmp');