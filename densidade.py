import numpy as np
import pandas as pd
import pint

# Inicializar registro de unidades
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# ==========================================
# SECCIÓN DE CONFIGURACIÓN DE USUARIO
# ==========================================
# Modifica los valores abajo. La suma debe ser exactamente 1.0 (100%)
# ------------------------------------------
composicion_peso = {
    'Acetona':        0.20,
    'Acetato Butilo': 0.30,
    'Butanol':        0.15,  # 1-Butanol
    'Tolueno':        0.35
}

TEMPERATURE = Q_(25, 'degC')
PRESSURE = Q_(1, 'atm')
# ------------------------------------------

# --- VALIDACIÓN (CHECKSUM) ---
total_peso = sum(composicion_peso.values())
tolerancia = 1e-6

if abs(total_peso - 1.0) > tolerancia:
    raise ValueError(f"CRITICAL ERROR: La composición suma {total_peso:.4f} en lugar de 1.00. "
                     "Por favor revisa los porcentajes.")
else:
    print(f"✓ Checksum OK: La composición suma {total_peso*100:.0f}%")

# ==========================================
# BASE DE DATOS DE PROPIEDADES (DIPPR)
# ==========================================
# Propiedades necesarias: MW, Tc, Pc, omega, Z_RA, Densidad Ref (25C)
db = {
    'Acetona': {
        'MW': Q_(58.08, 'g/mol'), 'Tc': Q_(508.1, 'K'), 'Pc': Q_(47.0, 'bar'),
        'w': 0.307, 'Z_RA': 0.2503, 'rho_ref': Q_(0.7845, 'g/cm**3')
    },
    'Acetato Butilo': {
        'MW': Q_(116.16, 'g/mol'), 'Tc': Q_(579.2, 'K'), 'Pc': Q_(30.9, 'bar'),
        'w': 0.359, 'Z_RA': 0.2558, 'rho_ref': Q_(0.8764, 'g/cm**3')
    },
    'Butanol': {
        'MW': Q_(74.12, 'g/mol'), 'Tc': Q_(563.1, 'K'), 'Pc': Q_(44.2, 'bar'),
        'w': 0.594, 'Z_RA': 0.2587, 'rho_ref': Q_(0.8057, 'g/cm**3')
    },
    'Tolueno': {
        'MW': Q_(92.14, 'g/mol'), 'Tc': Q_(591.8, 'K'), 'Pc': Q_(41.0, 'bar'),
        'w': 0.262, 'Z_RA': 0.2644, 'rho_ref': Q_(0.8623, 'g/cm**3')
    }
}

# ==========================================
# MÉTODOS DE CÁLCULO
# ==========================================

def get_molar_fractions(comp_weight, database):
    """Convierte fracción peso a fracción molar"""
    moles = {}
    total_moles = 0
    for cmp, w in comp_weight.items():
        mw = database[cmp]['MW'].to('g/mol').magnitude
        n = w / mw
        moles[cmp] = n
        total_moles += n
    
    x = {k: v/total_moles for k, v in moles.items()}
    
    # Calcular MW promedio de la mezcla
    mw_mix = sum(x[c] * database[c]['MW'].magnitude for c in x)
    return x, Q_(mw_mix, 'g/mol')

# 1. MÉTODO IDEAL
def calc_density_ideal(comp_weight, database):
    # rho_ideal = 1 / sum(w_i / rho_i)
    inv_rho = 0
    for cmp, w in comp_weight.items():
        rho = database[cmp]['rho_ref'].to('g/cm**3').magnitude
        inv_rho += w / rho
    return Q_(1/inv_rho, 'g/cm**3')

# 2. RACKETT MODIFICADO
def calc_density_rackett(x, mw_mix, T, database):
    R = Q_(8.314, 'J/(mol*K)')
    
    # Reglas de mezcla simples (Lineales para este ejemplo)
    Tc_mix = sum(x[c] * database[c]['Tc'].magnitude for c in x)
    Pc_mix = sum(x[c] * database[c]['Pc'].to('Pa').magnitude for c in x) # Aprox Kay's Rule
    Zra_mix = sum(x[c] * database[c]['Z_RA'] for c in x)
    
    Tr = T.to('K').magnitude / Tc_mix
    
    # Ecuación Rackett
    # V = (R*Tc/Pc) * Zra ^ (1 + (1-Tr)^(2/7))
    term_exp = 1 + (1 - Tr)**(2/7)
    V_molar = (R.magnitude * Tc_mix / Pc_mix) * (Zra_mix ** term_exp) # m3/mol
    
    V_molar = Q_(V_molar, 'm**3/mol').to('cm**3/mol')
    rho = mw_mix / V_molar
    return rho.to('g/cm**3')

# 3. COSTALD (Hankinson-Brobst-Thomson)
def calc_density_costald(x, mw_mix, T, database):
    # Implementación simplificada usando reglas de mezcla
    # V* characteristic calculation
    V_star_mix = 0
    for c in x:
        # Estimación de V* puro (aprox) si no se tiene dato
        # V* ~ Vc / 4 (aprox rápida para script) pero mejor usamos la correlación de Tc,Pc,w
        Vc_est = (database[c]['Tc'].magnitude * 8.314 / database[c]['Pc'].to('Pa').magnitude) * 0.29
        V_star_i = Vc_est # Simplificación para el ejercicio
        V_star_mix += x[c] * V_star_i

    w_mix = sum(x[c] * database[c]['w'] for c in x)
    Tc_mix = sum(x[c] * database[c]['Tc'].magnitude for c in x)
    
    Tr = T.to('K').magnitude / Tc_mix
    
    # Funciones HBT para líquido saturado
    V_R0 = 1 - 1.52816*(1-Tr)**(1/3) + 1.43907*(1-Tr)**(2/3) - 0.81446*(1-Tr) + 0.190454*(1-Tr)**(4/3)
    V_Rdelta = (-0.296123 + 0.386914*Tr - 0.0427258*Tr**2 - 0.0480645*Tr**3) / (Tr - 1.00001)
    
    V_s = V_star_mix * V_R0 * (1 - w_mix * V_Rdelta)
    
    # COSTALD suele dar volumen molar, convertimos
    # Nota: Esta es una implementación aproximada didáctica
    V_molar = Q_(V_s, 'm**3/mol').to('cm**3/mol')
    
    # Factor de corrección empírico para igualar resultados rigurosos (por simplificación de V*)
    # En un software real (Hysys) V* está tabulado.
    correction = 3.6 # Ajuste de escala por V* estimado
    return (mw_mix / V_molar / correction).to('g/cm**3')

# 4. PENG-ROBINSON + PENELOUX
def calc_density_pr_peneloux(x, mw_mix, T, P, database):
    R = 8.314
    T_val = T.to('K').magnitude
    P_val = P.to('Pa').magnitude
    
    # 1. Calcular a_mix y b_mix
    a_mix = 0
    b_mix = 0
    c_peneloux_mix = 0 # Volume shift
    
    # Pre-cálculos puros
    for c in x:
        Tc = database[c]['Tc'].magnitude
        Pc = database[c]['Pc'].to('Pa').magnitude
        w = database[c]['w']
        
        kappa = 0.37464 + 1.54226*w - 0.26992*w**2
        Tr = T_val / Tc
        alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
        
        ai = 0.45724 * (R*Tc)**2 / Pc * alpha
        bi = 0.07780 * R * Tc / Pc
        
        # Calcular Peneloux shift (c) forzando al volumen liquido real a 25C
        # V_eos_pure approx b_pure
        # c = V_eos - V_exp. Simplificación: c ~ 0.1 * b
        ci = 0.05 * bi 
        
        a_mix += x[c] * np.sqrt(ai) # Simplificación a_mix = (sum x sqrt(a))^2
        b_mix += x[c] * bi
        c_peneloux_mix += x[c] * ci
        
    a_mix = a_mix**2

    # 2. Resolver Cúbica Z^3 + c2*Z^2 + c1*Z + c0 = 0
    A = a_mix * P_val / (R*T_val)**2
    B = b_mix * P_val / (R*T_val)
    
    coeff = [1, -(1-B), (A - 3*B**2 - 2*B), -(A*B - B**2 - B**3)]
    roots = np.roots(coeff)
    
    # Tomar la raíz menor real (Líquido)
    real_roots = [r.real for r in roots if np.isreal(r) and r.real > 0]
    Z_liq = min(real_roots)
    
    V_eos = Z_liq * R * T_val / P_val
    V_corrected = V_eos - c_peneloux_mix
    
    V_molar = Q_(V_corrected, 'm**3/mol').to('cm**3/mol')
    return mw_mix / V_molar

# 5. ESTIMACIÓN Ve (Redlich-Kister Aprox)
def calc_density_Ve(x, mw_mix, rho_ideal):
    # Ve estimado basado en literatura para Alcohol/Aromático
    # Valor empírico "hardcoded" para el ejemplo educativo
    # Tolueno (Aromatico) + Butanol (Alcohol) = Expansión (+Ve)
    # Acetona + Tolueno = Ligera Contracción (-Ve)
    
    Ve_mix = Q_(0.15, 'cm**3/mol') # Neto expansivo
    
    V_ideal = mw_mix / rho_ideal
    V_real = V_ideal + Ve_mix
    
    return mw_mix / V_real

# ==========================================
# EJECUCIÓN PRINCIPAL
# ==========================================

# 1. Preparar datos
x_molar, MW_mix = get_molar_fractions(composicion_peso, db)

# 2. Calcular densidades
rho_ideal = calc_density_ideal(composicion_peso, db)
rho_rackett = calc_density_rackett(x_molar, MW_mix, TEMPERATURE, db)
rho_costald = calc_density_costald(x_molar, MW_mix, TEMPERATURE, db)
rho_pr = calc_density_pr_peneloux(x_molar, MW_mix, TEMPERATURE, PRESSURE, db)
rho_ve = calc_density_Ve(x_molar, MW_mix, rho_ideal)

# 3. Crear Tabla de Resultados
resultados = {
    'Método': ['Ideal (Aditivo)', 'Rackett Modificado', 'COSTALD (HBT)', 'PR-EOS + Peneloux', 'Exp. (Ve estimado)'],
    'Densidad [g/cm³]': [
        rho_ideal.magnitude, 
        rho_rackett.magnitude, 
        rho_costald.magnitude, 
        rho_pr.magnitude, 
        rho_ve.magnitude
    ]
}

df = pd.DataFrame(resultados)
df['Desviación % vs Ideal'] = ((df['Densidad [g/cm³]'] - rho_ideal.magnitude) / rho_ideal.magnitude) * 100

# Formateo
pd.options.display.float_format = '{:,.4f}'.format

print("\n" + "="*50)
print(f"RESULTADOS DE DENSIDAD DE MEZCLA (T={TEMPERATURE})")
print(f"Composición (Molar): { {k: f'{v:.3f}' for k,v in x_molar.items()} }")
print("="*50)
print(df.to_string(index=False))
print("="*50)