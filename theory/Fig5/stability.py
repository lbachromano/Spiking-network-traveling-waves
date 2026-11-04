"""Translation of the uploaded MATLAB stability calculation.

Units and names follow the MATLAB sources: seconds, mV, rates in Hz;
``sigma_E/I`` store variances. See README.md for deliberate numerical changes.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace
import warnings
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import erfcx, zeta
from scipy.io import savemat



ALPHA = np.sqrt(2.0) * abs(zeta(0.5))

def _statistics(all_nu, p):
    e, i = all_nu
    mu_e = p.tau_me*(p.Je_ext*p.nu_ext+p.we*p.C_EE*e-p.we*p.C_EI*p.ge*i)
    mu_i = p.tau_mi*(p.Ji_ext*p.nu_ext+p.wi*p.C_IE*e-p.wi*p.C_II*p.gi*i)
    ea = p.tau_me*(p.Je_ext**2*p.nu_ext+p.we**2*p.C_EE*e)
    eg = p.tau_me*p.we**2*p.ge**2*p.C_EI*i
    ia = p.tau_mi*(p.Ji_ext**2*p.nu_ext+p.wi**2*p.C_IE*e)
    ig = p.tau_mi*p.wi**2*p.gi**2*p.C_II*i
    te = (ea+eg)/(ea/p.tilde_AMPA_e+eg/p.tilde_GABA_e)
    ti = (ia+ig)/(ia/p.tilde_AMPA_i+ig/p.tilde_GABA_i)
    return SimpleNamespace(mu_E=mu_e, mu_I=mu_i, sigma_E=ea+eg, sigma_I=ia+ig,
                           nu_E=e, nu_I=i, t_syn_e=te, t_syn_i=ti)


def _bounds(mu, variance, ts, tm, p):
    shift = ALPHA*np.sqrt(ts/tm)/2
    return (p.hvr-mu)/np.sqrt(variance)+shift, (p.theta-mu)/np.sqrt(variance)+shift


def _rate(lower, upper, tm, refractory):
    # erfcx(-y) == exp(y*y)*(1+erf(y)), without negative-y cancellation.
    integral = quad(lambda y: erfcx(-y), lower, upper, epsabs=1e-10, epsrel=1e-6)[0]
    return 1/(refractory+tm*np.sqrt(np.pi)*integral)


def _rate_residual(all_nu, param):
    r = _statistics(all_nu, param)
    e = _rate(*_bounds(r.mu_E,r.sigma_E,r.t_syn_e,param.tau_me,param),param.tau_me,param.taurpe)
    i = _rate(*_bounds(r.mu_I,r.sigma_I,r.t_syn_i,param.tau_mi,param),param.tau_mi,param.taurpi)
    return np.array([e-r.nu_E, i-r.nu_I])


def Compute_nu(all_nu, param):
    residual = _rate_residual(all_nu, param)
    return float(residual @ residual)


def Mean_field_rates_IC(param, initial_nu, *, return_result=False):
    """Same bounded objective as fmincon; SciPy SLSQP is a different optimizer.

    The MATLAB options variable is not passed to fmincon. Python settings below
    are explicit numerical choices, not a claim of identical MATLAB defaults.
    """
    result = minimize(Compute_nu, np.asarray(initial_nu, dtype=float), args=(param,),
                      method='SLSQP', bounds=[(1e-6,500)]*2,
                      options={'ftol':1e-12, 'maxiter':1000})
    if not result.success:
        raise RuntimeError(f'Mean-field optimization failed: {result.message}')
    if result.fun > 1e-8:
        warnings.warn(f'Mean-field optimizer stopped with squared residual {result.fun:g}', RuntimeWarning)
    r = _statistics(result.x, param)
    values = np.array([r.mu_E,r.mu_I,r.sigma_E,r.sigma_I,r.nu_E,r.nu_I,r.t_syn_e,r.t_syn_i])
    return (values, result) if return_result else values


def Set_rates_IC(param, initial_nu):
    values = Mean_field_rates_IC(param, initial_nu)
    return SimpleNamespace(**dict(zip(('mu_E','mu_I','sigma_E','sigma_I','nu_E','nu_I','t_syn_e','t_syn_i'), values)))


def _derivatives(param, rates, population):
    if population == 'E':
        mu,v,ts,tm,tr = rates.mu_E,rates.sigma_E,rates.t_syn_e,param.tau_me,param.taurpe
    else:
        mu,v,ts,tm,tr = rates.mu_I,rates.sigma_I,rates.t_syn_i,param.tau_mi,param.taurpi
    mu,v,ts = np.broadcast_arrays(mu,v,ts)
    dmu,dsig = np.empty(mu.shape),np.empty(mu.shape)
    for idx in np.ndindex(mu.shape):
        lo,hi = _bounds(mu[idx],v[idx],ts[idx],tm,param)
        nu = _rate(lo,hi,tm,tr)
        el,eu = erfcx(-lo),erfcx(-hi)
        dmu[idx] = tm*nu**2*np.sqrt(np.pi)*(eu-el)/np.sqrt(v[idx])
        # dtds = 0. This output is dnu/d(sqrt(variance)), as in MATLAB.
        der_s = eu*(-(param.theta-mu[idx])/v[idx])-el*(-(param.hvr-mu[idx])/v[idx])
        dsig[idx] = -tm*nu**2*np.sqrt(np.pi)*der_s
    return dmu,dsig


def Compute_derivatives_E(param,rates):
    return _derivatives(param,rates,'E')


def Compute_derivatives_I(param,rates):
    return _derivatives(param,rates,'I')


def _broadcast_scalar_function(fun, *args, two_outputs=False):
    arrays = np.broadcast_arrays(*args)
    if arrays[0].ndim == 0:
        return fun(*(a.item() for a in arrays))
    out = np.empty(arrays[0].shape,complex)
    other = np.empty(arrays[0].shape,int) if two_outputs else None
    for idx in np.ndindex(out.shape):
        value = fun(*(a[idx].item() for a in arrays))
        if two_outputs:
            out[idx],other[idx] = value
        else:
            out[idx] = value
    return (out,other) if two_outputs else out


@lru_cache(maxsize=32768)
def _gsurg(xx, yy):
    xx,yy = complex(xx),complex(yy)
    cof = (76.18009172947146,-86.50532032941677,24.01409824083091,
           -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5)
    def parts(z):
        if z.real <= 0: return z+1.5,z+6.5,z+2
        if z.real <= 1: return z+0.5,z+5.5,z+1
        return z-0.5,z+4.5,z
    xd,xg,x = parts(xx)
    yd,yg,y = parts(yy)
    ex = np.log(abs(xg))+1j*np.arcsin(xx.imag/abs(xg))
    ey = np.log(abs(yg))+1j*np.arcsin(yy.imag/abs(yg))
    nano = np.exp(ex*xd-xg-ey*yd+yg)
    sx=sy=1.000000000190015+0j
    for coefficient in cof:
        sx += coefficient/x
        sy += coefficient/y
        x += 1
        y += 1
    nano *= sx/sy
    if xx.real <= 0: nano /= xx*(1+xx)
    elif xx.real <= 1: nano /= xx
    if yy.real <= 0: nano *= yy*(1+yy)
    elif yy.real <= 1: nano *= yy
    return nano


def gsurg(xx,yy):
    return _broadcast_scalar_function(_gsurg,xx,yy)


def onefone(a,c,z):
    """MATLAB recurrence, absolute 1e-12 stop and 1000-iteration cap."""
    a,c,z = np.broadcast_arrays(np.asarray(a,complex),np.asarray(c,complex),np.asarray(z,complex))
    fac=np.ones(a.shape,complex)
    temp=fac.copy()
    deriv=np.zeros(a.shape,complex)
    aa,cc=a.copy(),c.copy()
    for n in range(1,1001):
        fac=fac*(aa/cc)
        deriv=deriv+fac
        fac=fac*(z/n)
        series=temp+fac
        if np.all(np.abs(series-temp)<1e-12):
            if n == 1000: raise RuntimeError('Convergence failure in onefone')
            return series,deriv
        temp=series
        aa=aa+1
        cc=cc+1
    raise RuntimeError('Convergence failure in onefone')


def correctionlargey(a,c,z):
    """C recurrence (temp=series), retaining MATLAB's absolute stop tolerance."""
    a,c,z = np.broadcast_arrays(np.asarray(a,complex),np.asarray(c,complex),np.asarray(z,complex))
    temp=np.ones(a.shape,complex)
    fac=temp.copy()
    aa,cc=a.copy(),a+1-c
    for n in range(1,11):
        fac=fac*(aa*cc)/(-n*z)
        series=temp+fac
        if np.all(np.abs(series-temp)<1e-12): break
        temp=series  # main.c; fixes MATLAB temp = temp .* series
        aa=aa+1
        cc=cc+1
    return series


@lru_cache(maxsize=65536)
def _u(om,y,derivative):
    om,y=complex(om),float(y)
    ah=0.5*om
    sqom=np.sqrt(2*om)
    if abs(om)<10.5 and y>-3.5:
        if derivative:
            series,_=onefone(1+ah,1.5,y*y)
            ut=y*om*series*_gsurg(ah,0.5+ah)
            series,_=onefone(0.5+ah,0.5,y*y)
            ut+=series
        else:
            series,_=onefone(ah,0.5,y*y)
            ut=0.5*series*_gsurg(ah,0.5+ah)
            series,_=onefone(0.5+ah,1.5,y*y)
            ut+=y*series
        return complex(2*np.sqrt(np.pi)*ut/_gsurg(ah,1)),1
    lower=-3.5-0.25*(abs(om)-10.5)
    upper=-3.5-0.15*(abs(om)-10.5)
    # C has strict inequalities that leave exact transition endpoints undefined.
    # Select the pure expansion at endpoints (the interpolation endpoint value).
    if abs(om)==10.5 and y==-3.5:
        raise ValueError('Undefined C asymptotic junction: abs(om)=10.5, y=-3.5')
    use_y = y <= upper
    use_om = y >= lower
    if use_y:
        ay=np.exp(-np.log(abs(y))*om)
        if derivative:
            ay=ay/abs(y)*correctionlargey(ah+1,1.5,y*y)*om
        else:
            ay=ay*correctionlargey(ah,0.5,y*y)
    if use_om:
        exponent=(y*sqom+(y**3/6-y/2)/sqom-0.24*y*y/sqom**2
                  +(-y**5/40+y**3/12+y/8)/sqom**3+(y**4/8-y*y/4)/sqom**4)
        ao=np.sqrt(np.pi)*np.exp(0.5*y*y)*np.exp(exponent)
        if derivative:
            prefac=(sqom+y+(0.5*y*y-0.5)/sqom-0.48*y/sqom**2
                    +(-y**4/8+y*y/4+1/8)/sqom**3+(y**3/2-y/2)/sqom**4)
            ao*=prefac
        ao/=_gsurg(0.5+ah,1)
    if y <= lower: return complex(ay),2
    if y >= upper: return complex(ao),3
    alpha=-((y+3.5)/(abs(om)-10.5)+0.15)/0.1
    return complex(alpha*ay+(1-alpha)*ao),4


def u01(om,y):
    return _broadcast_scalar_function(lambda o,v:_u(o,v,False),om,y,two_outputs=True)


def u23(om,y):
    return _broadcast_scalar_function(lambda o,v:_u(o,v,True),om,y,two_outputs=True)


def rmu(omt,thetatil,hvrtil):
    omt=np.asarray(omt,complex)
    denominator=u01(omt,thetatil)[0]-u01(omt,hvrtil)[0]
    numerator=u23(omt,thetatil)[0]-u23(omt,hvrtil)[0]
    return numerator/((1+omt)*denominator)


def rsi(omt,thetatil,hvrtil):
    omt=np.asarray(omt,complex)
    denominator=u01(omt,thetatil)[0]-u01(omt,hvrtil)[0]
    numerator=thetatil*u23(omt,thetatil)[0]-hvrtil*u23(omt,hvrtil)[0]
    return omt/(2+omt)+numerator/((2+omt)*denominator)


