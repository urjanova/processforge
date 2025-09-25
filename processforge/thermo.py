import CoolProp.CoolProp as CP

def get_enthalpy(fluid, T, P):
    """Return molar enthalpy [J/mol]."""
    return CP.PropsSI("H", "T", T, "P", P, fluid)

def get_vapor_fraction(fluid, T, P):
    """Return vapor fraction (quality)."""
    return CP.PropsSI("Q", "T", T, "P", P, fluid)


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
        beta_new = beta - f/df if df != 0 else beta
        if abs(beta_new - beta) < tol:
            return max(0, min(1, beta_new))
        beta = beta_new
    return max(0, min(1, beta))
def get_enthalpy_molar(mixture, T, P):
    """
    Compute molar enthalpy [J/mol] for a pseudo-mixture.
    mixture: dict {component: mole fraction}
    """
    H = 0.0
    for comp, frac in mixture.items():
        H += frac * CP.PropsSI("HMOLAR", "T", T, "P", P, comp)
    return H