function omega = bestSOROmega(A)
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    lambda = max(abs(eig(-(D \ (L + U)))));
    omega = 2 / (1 + sqrt(1 - lambda * lambda));
end
