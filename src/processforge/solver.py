from scipy.integrate import solve_ivp

class Solver:
    def integrate(self, func, y0, t_span, t_eval):
        return solve_ivp(func, t_span, y0, t_eval=t_eval, method="BDF")
