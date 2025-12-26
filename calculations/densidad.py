# calculations/densidad.py
import numpy as np

def densidad_ideal(fracciones_peso, densidades_puras):
    """
    Calcula la densidad ideal asumiendo volúmenes aditivos.
    Utiliza fracciones en peso: rho_mix = 1 / sum(wi / rhoi)
    """
    try:
        inv_rho = sum(w / rho for w, rho in zip(fracciones_peso, densidades_puras))
        return 1.0 / inv_rho
    except ZeroDivisionError:
        return np.nan

def modified_rackett(T, x, props):
    """
    Implementa la ecuación de Rackett siguiendo el modelo de densidade.py.
    T: Kelvin, x: fracciones molares
    """
    try:
        R = 8.314  # J/(mol*K)
        
        # Reglas de mezcla (Kay's Rule)
        Tc_mix = sum(xi * p['Tc'] for xi, p in zip(x, props))
        Pc_mix = sum(xi * p['Pc'] * 1e5 for xi, p in zip(x, props)) # bar -> Pa
        Zra_mix = sum(xi * p['Z_RA'] for xi, p in zip(x, props))
        MW_mix = sum(xi * p['MW'] for xi, p in zip(x, props))

        Tr = T / Tc_mix
        if Tr >= 1.0: return np.nan
        
        # Ecuación de Rackett
        term_exp = 1 + (1 - Tr)**(2/7)
        V_molar = (R * Tc_mix / Pc_mix) * (Zra_mix ** term_exp) # m3/mol
        
        # m3/mol -> cm3/mol
        V_molar_cm3 = V_molar * 1e6
        
        # g/mol / cm3/mol * 1000 = kg/m3
        density_kgm3 = (MW_mix / V_molar_cm3) * 1000
        return density_kgm3
    except:
        return np.nan

def costald_density(T, x, props):
    """
    Implementación de COSTALD (HBT) siguiendo la lógica de densidade.py.
    """
    try:
        MW_mix = sum(xi * p['MW'] for xi, p in zip(x, props))
        Tc_mix = sum(xi * p['Tc'] for xi, p in zip(x, props))
        w_mix = sum(xi * p['Omega'] for xi, p in zip(x, props))
        
        # Estimación de V* mix (según lógica en densidade.py)
        V_star_mix = 0
        for xi, p in zip(x, props):
            # V* ~ Vc * 0.29
            Pc_pa = p['Pc'] * 1e5
            Vc_est = (p['Tc'] * 8.314 / Pc_pa) * 0.29
            V_star_mix += xi * Vc_est

        Tr = T / Tc_mix
        if Tr >= 1.0: return np.nan

        # Funciones HBT para líquido saturado
        V_R0 = 1 - 1.52816*(1-Tr)**(1/3) + 1.43907*(1-Tr)**(2/3) - 0.81446*(1-Tr) + 0.190454*(1-Tr)**(4/3)
        V_Rdelta = (-0.296123 + 0.386914*Tr - 0.0427258*Tr**2 - 0.0480645*Tr**3) / (Tr - 1.00001)
        
        V_s = V_star_mix * V_R0 * (1 - w_mix * V_Rdelta)
        
        # Corrección empírica sugerida en densidade.py
        correction = 3.6
        V_molar_cm3 = (V_s / correction) * 1e6
        
        density_kgm3 = (MW_mix / V_molar_cm3) * 1000
        return density_kgm3
    except:
        return np.nan

def pr_peneloux_density(T, x, props, P=1.01325):
    """
    Peng-Robinson EOS con Volume Shift (Peneloux) según densidade.py.
    P: Presión en bar (default 1 atm)
    """
    try:
        R = 8.314
        P_pa = P * 1e5
        MW_mix = sum(xi * p['MW'] for xi, p in zip(x, props))
        
        a_mix = 0
        b_mix = 0
        c_shift_mix = 0
        
        for xi, p in zip(x, props):
            Tc = p['Tc']
            Pc = p['Pc'] * 1e5
            w = p['Omega']
            
            kappa = 0.37464 + 1.54226*w - 0.26992*w**2
            Tr = T / Tc
            alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
            
            ai = 0.45724 * (R*Tc)**2 / Pc * alpha
            bi = 0.07780 * R * Tc / Pc
            ci = 0.05 * bi # Peneloux shift
            
            a_mix += xi * np.sqrt(ai)
            b_mix += xi * bi
            c_shift_mix += xi * ci
            
        a_mix = a_mix**2
        
        A = a_mix * P_pa / (R*T)**2
        B = b_mix * P_pa / (R*T)
        
        coeff = [1, -(1-B), (A - 3*B**2 - 2*B), -(A*B - B**2 - B**3)]
        roots = np.roots(coeff)
        
        # Raíz menor real para la fase líquida
        real_roots = [r.real for r in roots if np.isreal(r) and r.real > 0]
        if not real_roots: return np.nan
        Z_liq = min(real_roots)
        
        V_eos = Z_liq * R * T / P_pa
        V_corrected = V_eos - c_shift_mix
        
        V_molar_cm3 = V_corrected * 1e6
        return (MW_mix / V_molar_cm3) * 1000
    except:
        return np.nan
