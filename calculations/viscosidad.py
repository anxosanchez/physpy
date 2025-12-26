# calculations/viscosidad.py
import numpy as np

def metodo_lineal(fracciones_molares, viscosidades_puras):
    """Promedio lineal (molar)."""
    return sum(xi * mu for xi, mu in zip(fracciones_molares, viscosidades_puras))

def metodo_arrhenius(fracciones_molares, viscosidades_puras):
    """ln(mu_mix) = sum(xi * ln(mu_i))"""
    try:
        ln_mu = sum(xi * np.log(mu) for xi, mu in zip(fracciones_molares, viscosidades_puras))
        return np.exp(ln_mu)
    except:
        return np.nan

def metodo_kendall_monroe(fracciones_molares, viscosidades_puras):
    """(mu_mix)^(1/3) = sum(xi * mu_i^(1/3))"""
    try:
        mu_13 = sum(xi * (mu**(1/3)) for xi, mu in zip(fracciones_molares, viscosidades_puras))
        return mu_13**3
    except:
        return np.nan

def metodo_grunberg_nissan(fracciones_molares, viscosidades_puras, names):
    """
    Grunberg-Nissan con parámetros de interacción binaria (Gij).
    Implementa la lógica del script viscosidade.py
    """
    try:
        # Parte ideal (Arrhenius)
        ln_ideal = sum(xi * np.log(mu) for xi, mu in zip(fracciones_molares, viscosidades_puras))
        
        # Parámetros Gij estimados basados en viscosidade.py
        # Esto es una simplificación; en una DB real Gij estarían tabulados.
        excess_term = 0
        n = len(names)
        for i in range(n):
            for j in range(n):
                if i != j:
                    g_val = 0.0
                    name_i = names[i]
                    name_j = names[j]
                    
                    # Simulación de parámetros según viscosidade.py
                    # Butanol-Tolueno / Butanol-Acetona etc.
                    interact_pair = {name_i, name_j}
                    
                    if {"Butanol", "Tolueno"}.issubset(interact_pair) or {"n-Butanol", "Tolueno"}.issubset(interact_pair):
                        g_val = -0.4
                    elif {"Acetona", "Butanol"}.issubset(interact_pair) or {"Acetona", "n-Butanol"}.issubset(interact_pair):
                        g_val = -0.2
                    elif {"Agua", "Etanol"}.issubset(interact_pair):
                        g_val = 0.1 # Interacción positiva (puentes H)
                        
                    excess_term += g_val * fracciones_molares[i] * fracciones_molares[j]
        
        ln_mu_mix = ln_ideal + excess_term
        return np.exp(ln_mu_mix)
    except:
        return np.nan
