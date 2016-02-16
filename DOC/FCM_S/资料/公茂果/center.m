function [V2]=center(U1,c,m1,I1,D1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���룺�����Ⱦ���U1,������c,ģ����m1,I1
%�������������V1
%���ܣ������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V2=zeros(1,c);
[m,n]=size(U1);
U2=U1;
U3=U1;
for i=1:c
    temp1=0;
    temp2=0;
    for j=1:m
        for k=1:n
            temp1=temp1+U2{j,k}(i,1)^m1*I1(j,k)*D1{j,k}(i,1);
            temp2=temp2+U3{j,k}(i,1)^m1*D1{j,k}(i,1);
        end
    end
    V2(i)=temp1/temp2;
end        