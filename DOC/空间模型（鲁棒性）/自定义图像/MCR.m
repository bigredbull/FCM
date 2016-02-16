%%%%%%采用误分的像素点数/总像素点数的办法求得最终的误分率

function [mcr] = MCR(A,B)

total = size(A,2);
A_i = unique(A);
A_c = length(A_i);
B_i = unique(B);
B_c = length(B_i);


idA = double (repmat(A,A_c,1) == repmat(A_i',1,total));
idB = double (repmat(B,B_c,1) == repmat(B_i',1,total));

Sa = zeros(A_c,1);
Sb = zeros(B_c,1);
for i = 1:A_c
    for j = 1:total
        if idA(i,j) == 1
            Sa(i,1) = Sa(i,1) + 1;
        end
    end
end


for i = 1:B_c
    for j = 1:total
        if idB(i,j) ==1
            Sb(i,1) = Sb(i,1) + 1;
        end
    end
end

C = reshape(Sa,1,A_c);
D = reshape(Sb,1,B_c);


err = 0;
for i = 1:B_c
    err = err + abs(D(1,i) - C(1,i));
end
 
 mcr = (err / 2) / total;
 fprintf('MCR = %d\n',mcr);
 
 
