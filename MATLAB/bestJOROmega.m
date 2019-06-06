function omega = bestJOROmega(A)
    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);
    lambda_max = max(eig(-(D \ (L + U))));
    lambda_min = min(eig(-(D \ (L + U))));
    omega = 2 / (2 - lambda_max - lambda_min);
end

