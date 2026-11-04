"""Spatial determinant translated from Function_determinant_abs_nov_EE.m.

All model expressions are preserved except the two noise-amplitude corrections,
which multiply the white-noise amplitudes as in the supplied stability.py.
Vectorized Gauss–Legendre quadrature replaces MATLAB integral2..
"""
import numpy as np


from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace
import warnings
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import erfcx, zeta
from functools import lru_cache
from numpy.polynomial.legendre import leggauss
from stability import Compute_derivatives_E, Compute_derivatives_I, rmu, rsi

# Quadrature nodes per spatial coordinate.
# This is independent of n_u and n_v.
SPATIAL_QUAD_ORDER = 64

ALPHA = np.sqrt(2.0) * abs(zeta(0.5))


@lru_cache(maxsize=8)
def _quadrature_grid(Lbox, order):
    """Nodes, product weights, and distances on [0, Lbox/2]^2."""
    nodes, weights = leggauss(order)

    half = Lbox / 2.0

    # Map nodes and weights from [-1, 1] to [0, half].
    positions = (nodes + 1.0) * half / 2.0
    weights = weights * half / 2.0

    ru, rv = np.meshgrid(positions, positions, indexing="ij")
    product_weights = weights[:, None] * weights[None, :]
    distance = np.hypot(ru, rv)

    return ru, rv, product_weights, distance


@lru_cache(maxsize=16)
def _weighted_spatial_kernel(Lbox, length, rho, u, v, order):
    """Precompute everything that does not depend on lambda."""
    ru, rv, weights, distance = _quadrature_grid(Lbox, order)

    scale = length * np.sqrt(1.0 - rho**2)
    normalization = 1.0 / (
        2.0 * np.pi * length**2 * np.sqrt(1.0 - rho**2)
    )

    density = normalization * np.exp(
        -np.sqrt(
            (1.0 - rho) * ru**2
            + (1.0 + rho) * rv**2
        ) / scale
    )

    # Factor 4 accounts for all four quadrants.
    weighted_kernel = (
        4.0
        * weights
        * density
        * np.cos(u * ru)
        * np.cos(v * rv)
    )

    return weighted_kernel, distance


def _spatial_kernel(lam, length, rho, param):
    """Finite-square spatial integral, evaluated as one complex sum."""
    weighted_kernel, distance = _weighted_spatial_kernel(
        float(param.Lbox),
        float(length),
        float(rho),
        float(param.u),
        float(param.v),
        int(SPATIAL_QUAD_ORDER),
    )

    delay = np.exp(
        -complex(lam) * distance / float(param.velocity)
    )

    return np.sum(weighted_kernel * delay)


# The hard coded parameters in response functions come from the fit, see Table 2 in the paper.

def Compute_responsee(rhf,A0,x):
    x=np.abs(x)
    return A0+(-0.2989*x**(1.8070-1)+rhf*x**1.8070)/(2.1096+x**1.8070)


def Compute_responsei(rhf,A0,x):
    x=np.abs(x)
    return A0+(-0.5295*x**(2.0571-1)+rhf*x**2.0571)/(4.2311+x**2.0571)


def Compute_noise_responsee(omega):
    """Multiplicative recovered-noise amplitude correction; omega in rad/s."""
    omega=np.abs(omega)
    Tz, Tp = 4.3682e-3, 10.8458e-3  # seconds AMPLITUDE_WEIGHT_BY_STD = True
    #Tz, Tp = 0.00387834, 0.01026986  # seconds AMPLITUDE_WEIGHT_BY_STD = False
    return np.sqrt((1+(omega*Tz)**2)/(1+(omega*Tp)**2))


def Compute_noise_responsei(omega):
    """Multiplicative recovered-noise amplitude correction; omega in rad/s."""
    omega=np.abs(omega)
    Tz, Tp = 3.1165e-3, 7.0674e-3  # seconds AMPLITUDE_WEIGHT_BY_STD = True
    #Tz, Tp = 0.00303727, 0.00691497  # seconds AMPLITUDE_WEIGHT_BY_STD = False
    return np.sqrt((1+(omega*Tz)**2)/(1+(omega*Tp)**2))


def Compute_phase_shift_e(Phi0,x):
    x=np.abs(x)
    return np.abs(Phi0*(1-x**1.0970/(1.9558+x**1.0970)))


def Compute_phase_shift_i(Phi0,x):
    x=np.abs(x)
    return np.abs(Phi0*(1-x**1.1220/(1.5426+x**1.1220)))

def Compute_phase_shift_noise_e(Phi0, x):
    x = np.abs(x)
    a, xc, r = 0.96769, 1.23200366, 0.88385182
    return (Phi0 + a * x**r) / (1 + x / xc)


def Compute_phase_shift_noise_i(Phi0, x):
    x = np.abs(x)
    a, xc, r = 0.56166482, 1.46741053, 0.77472655
    return (Phi0 + a * x**r) / (1 + x / xc)

 
def Function_determinant_abs_nov_EE(cplx_om, param, rates):
    """Return log(abs(det) + 1e-7) for one complex eigenvalue.

    Spatial lengths, Lbox and velocity must use consistent length units;
    all times are seconds. The supplied response formulas require nonzero
    imaginary frequency. The wrapper searches only positive frequencies.
    """
    cplx_om = complex(cplx_om)
    if cplx_om.imag < 0:
        cplx_om = cplx_om.conjugate()
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
    [DmuE, DsigE] = Compute_derivatives_E(param,rates)
    [DmuI, DsigI] = Compute_derivatives_I(param,rates)

    DsigE = DsigE/(2*np.sqrt(sigma_e_tot))   # now dν/dQ
    DsigI = DsigI/(2*np.sqrt(sigma_i_tot))
    DqE, DqI = DsigE, DsigI
    
    S_eA=1/((1+cplx_om*param.tau_r_Ae) * (1+cplx_om * param.tau_d_Ae))
    S_eG=1/((1+cplx_om*param.tau_r_Ge) * (1+cplx_om * param.tau_d_Ge))
    S_iA=1/((1+cplx_om*param.tau_r_Ai) * (1+cplx_om * param.tau_d_Ai))
    S_iG=1/((1+cplx_om*param.tau_r_Gi) * (1+cplx_om * param.tau_d_Gi))
    Sn_eA=(param.tau_me/param.tilde_AMPA_e)/(2+cplx_om*(param.tilde_AMPA_e))
    Sn_eG=(param.tau_me/param.tilde_GABA_e)/(2+cplx_om*(param.tilde_GABA_e))
    Sn_iA=(param.tau_mi/param.tilde_AMPA_i)/(2+cplx_om*(param.tilde_AMPA_i))
    Sn_iG=(param.tau_mi/param.tilde_GABA_i)/(2+cplx_om*(param.tilde_GABA_i))
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
    
    Rsig0_E = DqE * (2.0 * t_syn_e / param.tau_me)
    Rsig0_I = DqI * (2.0 * t_syn_i / param.tau_mi)

    Rn_E = (Rsig0_E * WsigE_lambda * amp_corr_sigE * np.exp(-1j * phase_corr_sigE))
    Rn_I = (Rsig0_I * WsigI_lambda* amp_corr_sigI * np.exp(-1j * phase_corr_sigI))


    # Local values preserve MATLAB's pass-by-value treatment of param.
    L_ee = _spatial_kernel(cplx_om, param.l_ee, param.rho_EE, param)
    L_ii = _spatial_kernel(cplx_om, param.l_ii, param.rho_II, param)
    L_ei = _spatial_kernel(cplx_om, param.l_ei, param.rho_EI, param)
    L_ie = _spatial_kernel(cplx_om, param.l_ie, param.rho_IE, param)

    
    L_eff_EE = np.exp(-cplx_om*param.tau0)*L_ee
    L_eff_II = np.exp(-cplx_om*param.tau0)*L_ii

    L_eff_EI = np.exp(-cplx_om*param.tau0)*L_ei
    L_eff_IE = np.exp(-cplx_om*param.tau0)*L_ie

    RrpE=0
    RrpI=0

    
    AEE= RrpE + R_E*(mu_e_A/rates.nu_E)*S_eA*L_eff_EE
    nAEE= Rn_E*(sig_e_A/rates.nu_E)*Sn_eA*L_eff_EE
    AEI= R_E *(mu_e_G/rates.nu_I)*S_eG*L_eff_EI
    nAEI=Rn_E*(sig_e_G/rates.nu_I)*Sn_eG*L_eff_EI
    AIE= R_I*(mu_i_A/rates.nu_E)*S_iA*L_eff_IE
    nAIE= Rn_I*(sig_i_A/rates.nu_E)*Sn_iA*L_eff_IE
    AII= RrpI + R_I*(mu_i_G/rates.nu_I)*S_iG*L_eff_II 
    nAII=Rn_I*(sig_i_G/rates.nu_I)*Sn_iG*L_eff_II

    det_complex=  (1 -  AEE -  nAEE) * (1 +  AII - nAII) - (AEI - nAEI) * (- AIE - nAIE)
    funz = np.log(np.abs(det_complex) + 0.0000001)
    return funz
