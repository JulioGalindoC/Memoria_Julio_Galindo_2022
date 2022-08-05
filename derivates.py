import sympy as sym

def partial_dr_dx(folds_num) :

    x_ , th_ = sym.symbols('x_ th_')
    r_min , r_med , r_max , c = sym.symbols('r_min r_med r_max c')


    gauss1 = (r_med - r_min) * sym.exp(-(x_**2/2/c**2))
    gauss2 = (r_max - r_med) * sym.exp(-(x_**2/2/c**2))

    r = gauss1 + gauss2 * sym.cos(folds_num * th_) + r_min

    derivative_fx = r.diff(x_)

    return derivative_fx

def partial_dr_dth(folds_num) :

    x_ , th_ = sym.symbols('x_ th_')
    r_min , r_med , r_max , c = sym.symbols('r_min r_med r_max c')

    gauss1 = (r_med - r_min) * sym.exp(-(x_**2/2/c**2))
    gauss2 = (r_max - r_med) * sym.exp(-(x_**2/2/c**2))

    r = gauss1 + gauss2 * sym.cos(folds_num * th_) + r_min

    derivative_fth = r.diff(th_)

    return derivative_fth

