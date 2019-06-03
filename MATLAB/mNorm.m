function y = mNorm(A, b)
    matrixL1Norm = norm(A, 1);
    matrixL2Norm = norm(A);
    matrixLinfNorm = norm(A, inf);
    vectorL1Norm = norm(b, 1);
    vectorL2Norm = norm(b);
    vectorLinfNorm = norm(b, inf);
    fprintf('The L1-Norm, L2-Norm, Linf-Norm of A: %f %f %f \n', matrixL1Norm, matrixL2Norm, matrixLinfNorm);
    fprintf('The L1-Norm, L2-Norm, Linf-Norm of b: %f %f %f \n', vectorL1Norm, vectorL2Norm, vectorLinfNorm);
end

