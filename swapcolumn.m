function A_c = swapcolumn(A,i,j)
    v = A(:,i);
    A(:,i) = A(:,j);
    A(:,j) = v;
    A_c = A;
end