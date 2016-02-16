function [U2,V2]=Initialization(I,c)

[counts,x]=imhist(I); %�Ҷ�Ƶ��
 figure(5);
 imhist(I);
[m,n]=size(I);%�����Ĵ�С
U2=cell(m,n);
for i=1:m
    for j=1:n
        U2{i,j}=zeros();
    end
end
e=0.001;
V1=zeros(1,c);     %V1�Ǿɵľ�������
V2=zeros(1,c);     %V2���µľ�������
U1=zeros(c,256);
m1=2;              %ȷ����Ȩָ��m
num=0;             %��ʼ��������
V1(1)=50;          %��ʼ����������V
V1(2)=150;
V1(3)=200;
V1(4)=255;
I=double(I);%ת������������
flag=1;
%%%%%% ��ʼ���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (flag==1&&num<500) 
%%%%%%%%%%%%%%% ���������Ⱦ��� %%%%%%%%%%%%%%%
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
%%%%%%%%%%%%% ��һ�� %%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:256
    temp=0;
    for i=1:c
        temp=temp+U1(i,j);
    end
    for i=1:c
         U1(i,j)=U1(i,j)/temp;
    end
end
%%%%%%%%%%%%%%%���¾�������%%%%%%%%%%%%%%%%%%%%%%
for i=1:c
    sum=0;sum1=0;
    for j=1:256
        sum=sum+counts(j)*U1(i,j)^m1;
        sum1=sum1+x(j)*counts(j)*U1(i,j)^m1;  
    end
    V2(i)=sum1/sum;
end
%%%%%%%%%%%%%%% ��ֹ�ж� %%%%%%%%%%%%%%%%%%%%%%
num=num+1;
    if max(abs(V2-V1))<e;%����ֹͣ����
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