def rachford_rice(z, K, tol=1e-8, max_iter=100):
    """Solve Rachford-Rice for vapor fraction β."""
    beta = 0.5  # initial guess
    for _ in range(max_iter):
        f = sum(z[i] * (K[i] - 1) / (1 + beta * (K[i] - 1)) for i in z)
        df = -sum(z[i] * (K[i] - 1) ** 2 / (1 + beta * (K[i] - 1)) ** 2 for i in z)
        if df == 0:
            break
        beta_new = beta - f / df
        if abs(beta_new - beta) < tol:
            beta = beta_new
            break
        beta = beta_new
    return max(0, min(1, beta))
