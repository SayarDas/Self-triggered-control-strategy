function Y = gen_inv(X)
    r = rank(X);
    n = size(X,1);
    m = size(X,2);
    [C,col] = rref(X);
    [H,row] = rref(X');
    for i = 1:r
        X = swapcolumn(X,i,col(i));
        X = swaprow(X,i,row(i));
    end
    L = eye(n);
    R = eye(m);
    for i = 1:r
        R = swapcolumn(R,i,col(i));
    end
    for i = 1:r
        L = swaprow(L,i,row(i));
    end
    %%% R*Ginv(A)*L is the matrix we are looking for
    Y = R*[inv(X(1:r,1:r)) zeros(r,n-r);zeros(m-r,r) zeros(m-r,n-r)]*L;
end
                                               