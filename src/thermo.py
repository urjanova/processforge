import CoolProp.CoolProp as CP

def get_enthalpy_molar(mixture, T, P):
    """
    Compute molar enthalpy [J/mol] for a pseudo-mixture.
    mixture: dict {component: mole fraction}
    """
    H = 0.0
    for comp, frac in mixture.items():
        # HMOLAR returns J/mol
        H += frac * CP.PropsSI("HMOLAR", "T", T, "P", P, comp)
    return H

def get_Cp_molar(mixture, T, P):
    """
    Compute molar heat capacity Cp,molar [J/mol/K] for a mixture.
    mixture: dict {component: mole fraction}
    """
    Cp = 0.0
    for comp, frac in mixture.items():
        # Cpmolar returns J/mol/K
        Cp += frac * CP.PropsSI("Cpmolar", "T", T, "P", P, comp)
    return Cp

def get_K_values(components, T, P):
    """Return K-values for each component at (T,P)."""
    Ks = {}
    for comp in components:
        try:
            fugL = CP.PropsSI("fugL", "T", T, "P", P, comp)
            fugV = CP.PropsSI("fugV", "T", T, "P", P, comp)
            Ks[comp] = fugL / fugV if fugV != 0 else 1.0
        except Exception:
            # Fallback: Wilson approximation or just set K=1
            Ks[comp] = 1.0
    return Ks

def rachford_rice(z, K, tol=1e-8, max_iter=100):
    """Solve Rachford-Rice for vapor fraction β."""
    beta = 0.5  # initial guess
    for _ in range(max_iter):
        f = sum(z[i]*(K[i]-1)/(1 + beta*(K[i]-1)) for i in z)
        df = -sum(z[i]*(K[i]-1)**2 / (1 + beta*(K[i]-1))**2 for i in z)
        if df == 0:
            break
        beta_new = beta - f/df
        if abs(beta_new - beta) < tol:
            beta = beta_new
            break
        beta = beta_new
    return max(0, min(1, beta))
