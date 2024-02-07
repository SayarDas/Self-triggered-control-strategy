function A_r = swaprow(A,i,j)
    v = A(i,:);
    A(i,:) = A(j,:);
    A(j,:) = v;
    A_r = A;
end