# calculations/tension_superficial.py
import numpy as np
from scipy import optimize

def regla_lineal_molar(x, sigmas):
    """sigma_mix = sum(xi * sigma_i)"""
    return sum(xi * s for xi, s in zip(x, sigmas))

def regla_lineal_volumetrica(phi, sigmas):
    """sigma_mix = sum(phi_i * sigma_i)"""
    return sum(fi * s for fi, s in zip(phi, sigmas))

def macleod_sugden(x, paracor_list, mw_list, rho_liquid_mix):
    """
    Estima la tensión superficial de la mezcla usando la correlación de Macleod-Sugden.
    x: fracciones molares
    paracor_list: lista de Paracores de componentes puros
    mw_list: pesos moleculares
    rho_liquid_mix: densidad de la mezcla en kg/m3 (se convierte a g/cm3)
    """
    try:
        mw_mix = sum(xi * mw for xi, mw in zip(x, mw_list))
        p_mix = sum(xi * p for xi, p in zip(x, paracor_list))
        
        # rho_L en g/cm3 para la fórmula estándar
        rho_l_gcm3 = rho_liquid_mix / 1000.0
        
        # sigma = [ P_mix * rho_L / MW_mix ]^4
        sigma = (p_mix * rho_l_gcm3 / mw_mix)**4
        return sigma
    except Exception as e:
        return np.nan

def sprow_prausnitz(x, sigmas, mws, rhos, T):
    """
    Modelo de Sprow-Prausnitz para tensión superficial de una mezcla ideal.
    Resuelve sum(xi * exp((sigma - sigma_i)*Ai / RT)) = 1
    """
    try:
        R = 8.314 # J/molK
        T_val = T # Kelvin
        
        # Pre-calculamos constantes
        N_A = 6.022e23
        
        def objective_func(sigma_guess):
            sum_x_surf = 0
            for xi, sigma_i, mw, rho in zip(x, sigmas, mws, rhos):
                # Conversión a SI (rho en kg/m3, mw en kg/mol)
                mw_kg = mw / 1000.0
                vm = mw_kg / rho # m3/mol
                Ai_SI = (vm)**(2/3) * (N_A)**(1/3)
                
                sig_diff_SI = (sigma_guess - sigma_i) * 1e-3 # mN/m -> J/m2
                
                term = (sig_diff_SI * Ai_SI) / (R * T_val)
                sum_x_surf += xi * np.exp(term)
                
            return sum_x_surf - 1.0

        # Rango de búsqueda razonable
        sig_min = min(sigmas)
        sig_max = max(sigmas)
        
        sigma_sol = optimize.brentq(objective_func, sig_min - 10, sig_max + 10)
        return sigma_sol
    except:
        return np.nan
