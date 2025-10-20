from scipy.integrate import solve_ivp

class Solver:
    def integrate(self, func, y0, t_span, t_eval):
        """
        Integrates the given ordinary differential equation (ODE) function over the specified time span
        using the Backward Differentiation Formula (BDF) method from scipy's solve_ivp.
        Parameters
        ----------
        func : callable
            The right-hand side of the ODE system. It must be a function of the form func(t, y),
            where t is a scalar and y is an array-like object.
        y0 : array_like
            Initial conditions for the ODE system.
        t_span : tuple of float
            The interval of integration (t0, tf).
        t_eval : array_like, optional
            Times at which to evaluate the solution. If None, the solver will choose its own points.
        Returns
        -------
        Bunch
            The solution object returned by scipy.integrate.solve_ivp, containing the time points,
            solution values, and other information about the integration.
        """

        return solve_ivp(func, t_span, y0, t_eval=t_eval, method="BDF")