 
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


def Set_parameters():
    p = SimpleNamespace()
    p.L_ee = p.L_ei = p.L_ie = p.L_ii = 1.0
    p.tau_me, p.tau_mi = 20e-3, 10e-3
    p.tau_L_eA = 2.57e-3
    p.tau_L_eG = p.tau_L_iA = p.tau_L_iG = 1.725e-3
    p.taurpe, p.taurpi = 2e-3, 1e-3
    p.tau_d_Ae = p.tau_d_Ai = 3.5e-3
    p.tau_r_Ae = p.tau_r_Ai = 0.7e-3
    p.tau_d_Ge = p.tau_d_Gi = 18e-3
    p.tau_r_Ge = p.tau_r_Gi = 1e-3
    p.tilde_AMPA_e = p.tau_d_Ae + p.tau_r_Ae
    p.tilde_AMPA_i = p.tau_d_Ai + p.tau_r_Ai
    p.tilde_GABA_e = p.tau_d_Ge + p.tau_r_Ge
    p.tilde_GABA_i = p.tau_d_Gi + p.tau_r_Gi
    p.C_EE = p.C_IE = 800
    p.C_EI = p.C_II = 200
    p.theta, p.hvr = 18.0, 11.0
    p.gi = 4.1
    p.ge = 1.07 * p.gi
    p.we, p.wi = 1.5 * 0.1050, 1.5 * 0.1750
    p.Ji_ext, p.Je_ext = 0.85, 0.455
    p.nutheta = p.C_EE*p.theta/(p.C_EE*p.tau_me*p.we)
    p.nu_ext = p.nutheta*0.1  # assigned before nutheta is overwritten, as in MATLAB
    def threshold_rate(tm, ts, j):
        a = (2.065/2)*j*np.sqrt(ts)  # source uses rounded alpha here only
        return ((a + np.sqrt(a*a + 4*tm*j*p.theta))/(2*tm*j))**2
    p.nt_alpha = (0.8*threshold_rate(p.tau_me, p.tilde_AMPA_e, p.Je_ext)
                  + 0.2*threshold_rate(p.tau_mi, p.tilde_AMPA_i, p.Ji_ext))
    p.nutheta = p.nt_alpha
    return p


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


def Compute_responsee(rhf,A0,x):
    x=np.abs(x)
    return A0+(-0.2989*x**(1.8070-1)+rhf*x**1.8070)/(2.1096+x**1.8070)


def Compute_responsei(rhf,A0,x):
    x=np.abs(x)
    return A0+(-0.5295*x**(2.0571-1)+rhf*x**2.0571)/(4.2311+x**2.0571)


def Compute_noise_responsee(omega):
    """Multiplicative recovered-noise amplitude correction; omega in rad/s."""
    omega=np.abs(omega)
    Tz, Tp = 4.3682e-3, 10.8458e-3  # seconds
    return np.sqrt((1+(omega*Tz)**2)/(1+(omega*Tp)**2))


def Compute_noise_responsei(omega):
    """Multiplicative recovered-noise amplitude correction; omega in rad/s."""
    omega=np.abs(omega)
    Tz, Tp = 3.1165e-3, 7.0674e-3  # seconds
    return np.sqrt((1+(omega*Tz)**2)/(1+(omega*Tp)**2))


def Compute_phase_shift_e(Phi0,x):
    x=np.abs(x)
    return np.abs(Phi0*(1-x**1.0970/(1.9558+x**1.0970)))


def Compute_phase_shift_i(Phi0,x):
    x=np.abs(x)
    return np.abs(Phi0*(1-x**1.1220/(1.5426+x**1.1220)))


def Compute_phase_shift_noise_e(Phi0,x):
    x=np.abs(x)
    b, xc, p = 1.2543, 0.7523, 1.0
    return (Phi0+b*x)/(1+(x/xc)**p)


def Compute_phase_shift_noise_i(Phi0,x):
    x=np.abs(x)
    b, xc, p = 0.75847414, 0.65089609, 1.0
    return (Phi0+b*x)/(1+(x/xc)**p)


def Function_determinant_abs_nov(cplx_om,param,rates):
    """log(abs(det)+1e-7); scalar or array lambda, scalar operating point.

    Corrections depend only on Im(lambda): this is the source's approximate,
    non-holomorphic continuation. No zero-frequency regularization is added.
    """
    if any(np.size(value) != 1 for value in vars(rates).values()):
        raise ValueError('All rates fields must be scalar')
    cplx_om = np.asarray(cplx_om,dtype=complex)
    cplx_om = np.where(cplx_om.imag < 0, np.conj(cplx_om), cplx_om)
    if np.any(cplx_om.imag == 0):
        raise ValueError('The supplied white-noise formulas are undefined at omega=0; use nonzero frequencies.')
    Vr=param.hvr
    theta=param.theta
    xnuE=rates.nu_E
    xnuI=rates.nu_I
    mu_i_ext=param.tau_mi *param.Ji_ext*param.nu_ext
    mu_i_A=param.tau_mi *param.wi*param.C_IE*xnuE
    mu_i_G=param.tau_mi *param.wi*param.C_II* param.gi*xnuI
    mu_e_ext=param.tau_me *param.Je_ext*param.nu_ext
    mu_e_A=param.tau_me *param.we*param.C_EE*xnuE
    mu_e_G=param.tau_me *param.we*param.C_EI* param.ge*xnuI
    sig_i_ext=param.tau_mi *param.Ji_ext**2*param.nu_ext
    sig_i_A=param.tau_mi *param.wi**2*param.C_IE*xnuE
    sig_i_G=param.tau_mi *param.wi**2*param.gi**2*param.C_II*xnuI
    sig_e_ext=param.tau_me *param.Je_ext**2*param.nu_ext
    sig_e_A=param.tau_me * param.we**2*param.C_EE*xnuE
    sig_e_G=param.tau_me *param.we**2*param.ge**2*param.C_EI*xnuI

    t_syn_i = (sig_i_ext + sig_i_A + sig_i_G) / ((sig_i_ext + sig_i_A) / param.tilde_AMPA_i+ sig_i_G / param.tilde_GABA_i)
    t_syn_e = (sig_e_ext + sig_e_A + sig_e_G) / ((sig_e_ext + sig_e_A) / param.tilde_AMPA_e+ sig_e_G / param.tilde_GABA_e)

    mu_i_l=mu_i_A-mu_i_G
    mu_e_l=mu_e_A-mu_e_G
    sigma_i_l=sig_i_A+sig_i_G
    sigma_e_l=sig_e_A+sig_e_G
    mu_i_tot=mu_i_l+mu_i_ext
    mu_e_tot=mu_e_l+mu_e_ext
    sigma_i_tot=sigma_i_l+sig_i_ext
    sigma_e_tot=sigma_e_l+sig_e_ext
    E_vr=(Vr-mu_e_tot)/np.sqrt(sigma_e_tot)
    E_th=(theta-mu_e_tot)/np.sqrt(sigma_e_tot)
    I_vr=(Vr-mu_i_tot)/np.sqrt(sigma_i_tot)
    I_th=(theta-mu_i_tot)/np.sqrt(sigma_i_tot)
    k_prime_E=np.sqrt(t_syn_e*mu_e_tot*mu_e_tot/(param.tau_me*sigma_e_tot))
    k_prime_I=np.sqrt(t_syn_i*mu_i_tot*mu_i_tot/(param.tau_mi*sigma_i_tot))
    DmuE, DsigE = Compute_derivatives_E(param,rates)
    DmuI, DsigI = Compute_derivatives_I(param,rates)


    #That factor 2 below converts standard deviation to drive variance.

    DsigE=DsigE/(2*np.sqrt(sigma_e_tot))
    DsigI=DsigI/(2*np.sqrt(sigma_i_tot))

    # Convert to derivatives with respect to filtered-current variance,
    # holding the effective synaptic time constant fixed.
    DsigE *= 2*t_syn_e / param.tau_me
    DsigI *= 2*t_syn_i / param.tau_mi

    
    S_eA=np.exp(-cplx_om*param.tau_L_eA)/((1+cplx_om*param.tau_r_Ae) * (1+cplx_om * param.tau_d_Ae))
    S_eG=np.exp(-cplx_om*param.tau_L_eG)/((1+cplx_om*param.tau_r_Ge) * (1+cplx_om * param.tau_d_Ge))
    S_iA=np.exp(-cplx_om*param.tau_L_iA)/((1+cplx_om*param.tau_r_Ai) * (1+cplx_om * param.tau_d_Ai))
    S_iG=np.exp(-cplx_om*param.tau_L_iG)/((1+cplx_om*param.tau_r_Gi) * (1+cplx_om * param.tau_d_Gi))
    Sn_eA=(param.tau_me/param.tilde_AMPA_e)*np.exp(-cplx_om*param.tau_L_eA)/(2+cplx_om*(param.tilde_AMPA_e))
    Sn_eG=(param.tau_me/param.tilde_GABA_e)*np.exp(-cplx_om*param.tau_L_eG)/(2+cplx_om*(param.tilde_GABA_e))
    Sn_iA=(param.tau_mi/param.tilde_AMPA_i)*np.exp(-cplx_om*param.tau_L_iA)/(2+cplx_om*(param.tilde_AMPA_i))
    Sn_iG=(param.tau_mi/param.tilde_GABA_i)*np.exp(-cplx_om*param.tau_L_iG)/(2+cplx_om*(param.tilde_GABA_i))
    lam = cplx_om
    omega_r = np.imag(lam)
    lambda_axis = 1j * omega_r
    om0 = 1j * 0.01
    Rmu0E  = np.abs(rmu(param.tau_me*om0,E_th,E_vr))
    Rmu0I  = np.abs(rmu(param.tau_mi*om0,I_th,I_vr))
    RSig0E = np.abs(rsi(param.tau_me*om0,E_th,E_vr))
    RSig0I = np.abs(rsi(param.tau_mi*om0,I_th,I_vr))
    WmuE_lambda = rmu(param.tau_me*lam,E_th,E_vr) / Rmu0E
    WmuI_lambda = rmu(param.tau_mi*lam,I_th,I_vr) / Rmu0I
    WsigE_lambda = rsi(param.tau_me*lam,E_th,E_vr) / RSig0E
    WsigI_lambda = rsi(param.tau_mi*lam,I_th,I_vr) / RSig0I
    WmuE_axis = rmu(param.tau_me*lambda_axis,E_th,E_vr) / Rmu0E
    WmuI_axis = rmu(param.tau_mi*lambda_axis,I_th,I_vr) / Rmu0I
    WsigE_axis =  rsi(param.tau_me*lambda_axis,E_th,E_vr) / RSig0E
    WsigI_axis = rsi(param.tau_mi*lambda_axis,I_th,I_vr) / RSig0I
    A_muE_0 = np.abs(WmuE_axis)
    Phi_muE_0 = -np.angle(WmuE_axis)
    A_muI_0 = np.abs(WmuI_axis)
    Phi_muI_0 = -np.angle(WmuI_axis)
    A_sigE_0 = np.abs(WsigE_axis)
    Phi_sigE_0 = -np.angle(WsigE_axis)
    A_sigI_0 = np.abs(WsigI_axis)
    Phi_sigI_0 = -np.angle(WsigI_axis)
    xE = omega_r * t_syn_e
    xI = omega_r * t_syn_i
    rhfE = 1.3238 * rates.nu_E / np.sqrt(sigma_e_tot) * np.sqrt(t_syn_e / param.tau_me) / DmuE
    rhfI = 1.3238 * rates.nu_I / np.sqrt(sigma_i_tot) * np.sqrt(t_syn_i / param.tau_mi) / DmuI
    A_muE_col = np.abs(Compute_responsee(rhfE,A_muE_0,xE))
    A_muI_col = np.abs(Compute_responsei(rhfI,A_muI_0,xI))
    Phi_muE_col = Compute_phase_shift_e(Phi_muE_0,xE)
    Phi_muI_col = Compute_phase_shift_i(Phi_muI_0,xI)
    amp_corr_muE = A_muE_col / A_muE_0
    phase_corr_muE = Phi_muE_col - Phi_muE_0
    amp_corr_muI = A_muI_col / A_muI_0
    phase_corr_muI = Phi_muI_col - Phi_muI_0
    R_E = DmuE * WmuE_lambda *        amp_corr_muE * np.exp(-1j * phase_corr_muE)
    R_I = DmuI * WmuI_lambda *        amp_corr_muI * np.exp(-1j * phase_corr_muI)

 ## EFFECTIVE APPROXIMATION   
    A_sigE_col = A_sigE_0 * Compute_noise_responsee(omega_r)
    A_sigI_col = A_sigI_0 * Compute_noise_responsei(omega_r)
    Phi_sigE_col =      Compute_phase_shift_noise_e(Phi_sigE_0,xE)
    Phi_sigI_col =      Compute_phase_shift_noise_i(Phi_sigI_0,xI)
    amp_corr_sigE = A_sigE_col / A_sigE_0
    phase_corr_sigE = Phi_sigE_col - Phi_sigE_0
    amp_corr_sigI = A_sigI_col / A_sigI_0
    phase_corr_sigI = Phi_sigI_col - Phi_sigI_0
    Rn_E = DsigE * WsigE_lambda *         amp_corr_sigE * np.exp(-1j * phase_corr_sigE)
    Rn_I = DsigI * WsigI_lambda *         amp_corr_sigI * np.exp(-1j * phase_corr_sigI)


## WHITE NOISE APPROXIMATION
 #   Rsig0_E = DsigE  
 #   Rsig0_I = DsigI  
 #   Rn_E = Rsig0_E * WsigE_lambda
 #   Rn_I = Rsig0_I * WsigI_lambda
                        
 ##STATIC APPROXIMATION
#    Rn_E = DsigE 
#    Rn_I = DsigI  

    RrpE=0
    RrpI=0
    AEE= RrpE + R_E*(mu_e_A/rates.nu_E)*S_eA*param.L_ee
    nAEE= Rn_E*(sig_e_A/rates.nu_E)*Sn_eA*param.L_ee
    AEI= R_E *(mu_e_G/rates.nu_I)*S_eG*param.L_ei
    nAEI=Rn_E*(sig_e_G/rates.nu_I)*Sn_eG*param.L_ei
    AIE= R_I*(mu_i_A/rates.nu_E)*S_iA*param.L_ie
    nAIE= Rn_I*(sig_i_A/rates.nu_E)*Sn_iA*param.L_ie
    AII= RrpI + R_I*(mu_i_G/rates.nu_I)*S_iG*param.L_ii 
    nAII=Rn_I*(sig_i_G/rates.nu_I)*Sn_iG*param.L_ii
    det_complex=  (1 -  AEE -  nAEE) * (1 +  AII - nAII) - (AEI - nAEI) * (- AIE - nAIE)
    funz = np.log(np.abs(det_complex) + 0.0000001)
    return funz


def _uniquetol_rows_indices(xy,tol=0.1):
    """Column-scaled, two-coordinate relative tolerance, sorted representatives.

    Nontransitive tolerance clusters can admit several valid representative
    sets. Exact MATLAB representative choices require a MATLAB comparison.
    """
    xy=np.asarray(xy,float).reshape(-1,2)
    if not len(xy): return np.empty(0,int)
    scale=np.max(np.abs(xy),axis=0)
    order=np.lexsort((xy[:,1],xy[:,0]))
    kept=[]
    for idx in order:
        if not kept or not np.any(np.all(np.abs(xy[kept]-xy[idx]) <= tol*scale,axis=1)):
            kept.append(idx)
    return np.asarray(kept,int)


def _round_matlab(x,decimals=3):
    scale=10.0**decimals
    return np.copysign(np.floor(np.abs(x)*scale+0.5),x)/scale


def _initial_prune(minima):
    if not len(minima): return np.empty((0,3))
    minima=minima[np.argsort(minima[:,2],kind='stable')]
    # Source divides by .05 THEN rounds to 3 decimals: not a .05 distance test.
    xy_scaled=_round_matlab(minima[:,:2]/np.array([0.05,0.05]),3)
    _,idx=np.unique(xy_scaled,axis=0,return_index=True)
    minima=minima[idx]
    return minima[minima[:,0]<=100]


def _find_minima(param,rates,starts,threshold=-10.0):
    def cost(x):
        # The MATLAB expression is undefined on the real axis. Mark invalid
        # trial points as unusable rather than inserting a new limit formula.
        if x[1] == 0: return np.inf
        try:
            with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
                value=float(Function_determinant_abs_nov(complex(*x),param,rates))
        except (ZeroDivisionError,OverflowError,FloatingPointError):
            return np.inf
        return value if np.isfinite(value) else np.inf
    found=[]
    for start in starts:
        result=minimize(cost,np.asarray(start,float),method='Nelder-Mead',
                        options={'maxiter':150,'maxfev':400,'xatol':1e-3,'fatol':1e-4,
                                 'adaptive':False,'disp':False})
        if result.success and result.fun<threshold:
            found.append([result.x[0],result.x[1]/(2*np.pi),result.fun])
    return np.asarray(found,float).reshape(-1,3)


def Compute_starting_points_fwd(param,initial_nu,input_ratio):
    """Return initial minima, updated rate guess and fallback start ranges."""
    param.nu_ext=param.nutheta*input_ratio
    rates=Set_rates_IC(param,initial_nu)
    initial_nu=np.array([rates.nu_E,rates.nu_I])
    starts=[(re,im) for re in np.linspace(-10,4,3)
            for im in 2*np.pi*np.linspace(3,20,3)]
    minima=_initial_prune(_find_minima(param,rates,starts,-10))
    return minima,initial_nu,np.linspace(-50,50,5),np.linspace(0.1*2*np.pi,80*2*np.pi,5)

def distinct_continuation_seeds(minima, re_tol=1e-3, freq_tol=1e-3):
    """Group numerical duplicates using absolute tolerances.

    Columns: Re(lambda) [1/s], frequency [Hz], log residual.
    Keep the candidate with the smallest residual in each group.
    """
    candidates = np.asarray(minima, dtype=float).reshape(-1, 3).copy()
    if not len(candidates):
        return candidates

    # The determinant maps negative frequencies to positive ones.
    candidates[:, 1] = np.abs(candidates[:, 1])

    candidates = candidates[
        np.argsort(candidates[:, 2], kind="stable")
    ]

    kept = []
    for row in candidates:
        duplicate = any(
            abs(row[0] - other[0]) <= re_tol
            and abs(row[1] - other[1]) <= freq_tol
            for other in kept
        )
        if not duplicate:
            kept.append(row)

    return np.asarray(kept, dtype=float).reshape(-1, 3)

def Wrap_minimization_forw_dec_OK(*,allg=(3.625,),N_pt=100,
                                 list_varying_param=None,output_dir='.',progress=True):
    """Run the original forward sweep and save one MATLAB-compatible .mat per g.

    Optional sweep arguments allow short checks; defaults match active MATLAB.
    Returns a list of dictionaries, one for each g.
    """
    param=Set_parameters()
    ratios=(np.linspace(0.7,1.5,N_pt) if list_varying_param is None
            else np.asarray(list_varying_param,float))
    if ratios.ndim != 1 or ratios.size == 0:
        raise ValueError('list_varying_param must be a nonempty one-dimensional sequence')
    output_dir=Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)
    results=[]
    for g in allg:
        param.gi=float(g)
        param.ge=1.07*param.gi
        initial_nu=np.array([0.001,0.02])
        previous,initial_nu,start_re_range,start_im_range=Compute_starting_points_fwd(
            param,initial_nu,ratios[0])
        minima_all=[]
        rate_history=[]
        for i,ratio in enumerate(ratios):
            param.nu_ext=param.nutheta*ratio
            rates=Set_rates_IC(param,initial_nu)
            initial_nu=np.array([rates.nu_E,rates.nu_I])
            rate_history.append(initial_nu.copy())


         #   if len(previous):
         #       previous = previous[[np.argmax(previous[:, 0])]] #this to keep the one with largest real part

    
            if len(previous):
                starts=[(row[0],row[1]*2*np.pi) for row in previous]
            else:
                starts=[(re,im) for re in start_re_range for im in start_im_range]

                
            minima=_find_minima(param,rates,starts)

            
          #  if len(minima): minima=minima[_uniquetol_rows_indices(minima[:,:2],0.1)]


            minima = distinct_continuation_seeds(minima)
            minima_all.append(minima)
            previous = minima
            
            if progress:
                print(
                    f"g={g:.3f}, {i+1}/{len(ratios)}, input={ratio:.6g}, "
                    f"rates=({rates.nu_E:.6g}, {rates.nu_I:.6g}) Hz, "
                    f"candidates={len(minima)}, "
                    f"distinct_seeds={len(previous)}",
                    flush=True,
                )
        cells=np.empty((len(ratios),1),object)
        for i,m in enumerate(minima_all): cells[i,0]=m
        path=output_dir/f'Theoretical_prediction_{g:.3f}.mat'
        savemat(path,{'list_varying_param':ratios.reshape(1,-1),'local_minima_all':cells})
        results.append({'g':g,'list_varying_param':ratios,'local_minima_all':minima_all,
                        'rate_history':np.asarray(rate_history),'output_path':path})
    return results
