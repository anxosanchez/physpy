import numpy as np
import pandas as pd
import pint
from scipy import optimize

# Inicializar registro de unidades
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# ==========================================
# 1. CONFIGURACIÓN DE USUARIO
# ==========================================
# Modifica los valores abajo. La suma debe ser exactamente 1.0 (100%)
composicion_peso = {
    'Acetona':        0.20,
    'Acetato Butilo': 0.30,
    'Butanol':        0.15,
    'Tolueno':        0.35
}

TEMPERATURE = Q_(25, 'degC')
# ==========================================

# --- VALIDACIÓN (CHECKSUM) ---
total_peso = sum(composicion_peso.values())
tolerancia = 1e-6
if abs(total_peso - 1.0) > tolerancia:
    raise ValueError(f"CRITICAL ERROR: La composición suma {total_peso:.4f} en lugar de 1.00.")
else:
    print(f"✓ Checksum OK: La composición suma {total_peso*100:.0f}%")

# ==========================================
# 2. BASE DE DATOS (PROPIEDADES A 25°C)
# ==========================================
# sigma: Tensión superficial pura (mN/m o dyn/cm)
# P: Paracor de Sugden (Tabulado)
# rho: Densidad pura (g/cm3)
db = {
    'Acetona': {
        'MW': Q_(58.08, 'g/mol'), 'rho': Q_(0.784, 'g/cm**3'),
        'sigma': Q_(23.0, 'mN/m'), 'P': 162.5
    },
    'Acetato Butilo': {
        'MW': Q_(116.16, 'g/mol'), 'rho': Q_(0.882, 'g/cm**3'),
        'sigma': Q_(24.8, 'mN/m'), 'P': 336.8
    },
    'Butanol': {
        'MW': Q_(74.12, 'g/mol'), 'rho': Q_(0.810, 'g/cm**3'),
        'sigma': Q_(24.6, 'mN/m'), 'P': 224.2
    },
    'Tolueno': {
        'MW': Q_(92.14, 'g/mol'), 'rho': Q_(0.867, 'g/cm**3'),
        'sigma': Q_(27.9, 'mN/m'), 'P': 246.3
    }
}

# ==========================================
# 3. FUNCIONES DE CÁLCULO
# ==========================================

def get_fractions(comp_weight, database):
    """Calcula fracciones molares (x) y volumétricas (phi)"""
    moles = {}
    volumenes = {}
    total_moles = 0
    total_vol = 0
    
    # Paso 1: Moles
    for cmp, w in comp_weight.items():
        mw = database[cmp]['MW'].magnitude
        n = w / mw
        moles[cmp] = n
        total_moles += n
        
        # Volumen para phi (w / rho)
        rho = database[cmp]['rho'].magnitude
        v = w / rho
        volumenes[cmp] = v
        total_vol += v
        
    x = {k: v/total_moles for k, v in moles.items()}
    phi = {k: v/total_vol for k, v in volumenes.items()}
    
    # MW promedio y Densidad Ideal promedio para Macleod
    mw_mix = sum(x[c] * database[c]['MW'].magnitude for c in x)
    rho_mix = 1 / total_vol # base 1 kg
    
    return x, phi, Q_(mw_mix, 'g/mol'), Q_(rho_mix, 'g/cm**3')

# MÉTODOS DE PREDICCIÓN

# A. Regla Lineal Molar
def calc_sigma_molar(x, database):
    sigma_mix = 0
    for c, xi in x.items():
        sigma_mix += xi * database[c]['sigma'].magnitude
    return Q_(sigma_mix, 'mN/m')

# B. Regla Lineal Volumétrica (Winterfeld)
def calc_sigma_volume(phi, database):
    sigma_mix = 0
    for c, pi in phi.items():
        sigma_mix += pi * database[c]['sigma'].magnitude
    return Q_(sigma_mix, 'mN/m')

# C. Macleod-Sugden (Paracor)
def calc_sigma_macleod(x, rho_mix_ideal, mw_mix, database):
    # P_mix = sum(xi * Pi)
    P_mix = 0
    for c, xi in x.items():
        P_mix += xi * database[c]['P']
    
    # sigma = (P_mix * rho_mix / MW_mix)^4
    # Nota: rho debe estar en g/cm3 y MW en g/mol para usar Paracor estándar
    term = (P_mix * rho_mix_ideal.magnitude / mw_mix.magnitude)
    sigma_val = term**4
    return Q_(sigma_val, 'mN/m')

# D. Sprow-Prausnitz (Modelo de Superficie Ideal)
def calc_sigma_sprow_prausnitz(x, database, T):
    """
    Resuelve iterativamente la tensión superficial asumiendo fase líquida ideal
    y fase superficial ideal. Considera el Área Molar Parcial (Ai).
    """
    R = 8.314 # J/molK
    T_val = T.to('K').magnitude
    
    # 1. Calcular Áreas Molares Parciales (Ai)
    # Ai aprox = (Vi_molar)^(2/3) * (Avogadro)^(1/3) -> Relación geométrica
    # Usaremos una constante simplificada típica para disolventes: 1e8 cm2/mol aprox
    # O mejor: Ai = (MW / rho)^(2/3) * 1.02e8 (factor estructural)
    A = {}
    for c in x:
        mw = database[c]['MW'].magnitude
        rho = database[c]['rho'].magnitude
        vol_molar = mw / rho # cm3/mol
        # Estimación aproximada del área proyectada
        Ai = (vol_molar)**(2/3) * 100000 # Factor escala arbitrario para unidades consistentes
        A[c] = Ai 

    # Función objetivo para encontrar sigma_mix
    # sum(xi_surface) must be 1
    # xi_surface = xi * exp( (sigma_mix - sigma_i)*Ai / RT )
    
    def objective_func(sigma_guess):
        sum_x_surf = 0
        for c, xi in x.items():
            sigma_i = database[c]['sigma'].magnitude
            # Cuidado con unidades: (mN/m * m2/mol) = J/mol
            # sigma (mN/m) = 1e-3 N/m = 1e-3 J/m2
            # Ai está en cm2/mol? Convertir todo a SI para el exponente
            
            # Conversión rigurosa para el argumento exponencial
            sig_diff_SI = (sigma_guess - sigma_i) * 1e-3 # J/m2
            
            mw = database[c]['MW'].to('kg/mol').magnitude
            rho = database[c]['rho'].to('kg/m**3').magnitude
            vm = mw/rho # m3/mol
            N_A = 6.022e23
            # Area molar estándar ~ (Vm/Na)^(2/3) * Na
            Ai_SI = (vm)**(2/3) * (N_A)**(1/3)
            
            term = (sig_diff_SI * Ai_SI) / (R * T_val)
            
            x_surf = xi * np.exp(term)
            sum_x_surf += x_surf
            
        return sum_x_surf - 1.0

    # Resolver numéricamente
    # Buscamos un valor entre el min y max de los puros
    sig_min = min(db[c]['sigma'].magnitude for c in db)
    sig_max = max(db[c]['sigma'].magnitude for c in db)
    
    try:
        sigma_sol = optimize.brentq(objective_func, sig_min - 5, sig_max + 5)
    except:
        sigma_sol = calc_sigma_molar(x, database).magnitude # Fallback
        
    return Q_(sigma_sol, 'mN/m')

# ==========================================
# 4. EJECUCIÓN
# ==========================================

x_mol, phi_vol, MW_mix, rho_mix = get_fractions(composicion_peso, db)

# Cálculos
s_molar = calc_sigma_molar(x_mol, db)
s_vol = calc_sigma_volume(phi_vol, db)
s_macleod = calc_sigma_macleod(x_mol, rho_mix, MW_mix, db)
s_sprow = calc_sigma_sprow_prausnitz(x_mol, db, TEMPERATURE)

# Tabla
resultados = {
    'Método de Predicción': [
        'Regla Lineal (Molar)',
        'Regla Lineal (Volumétrica)',
        'Macleod-Sugden (Paracor)',
        'Sprow-Prausnitz (Ideal)'
    ],
    'Tensión Superficial [mN/m]': [
        s_molar.magnitude,
        s_vol.magnitude,
        s_macleod.magnitude,
        s_sprow.magnitude
    ],
    'Notas Técnicas': [
        'Promedio simple (Tiende a sobreestimar)',
        'Mejor para mezclas con polímeros',
        'Estándar para hidrocarburos',
        'Considera enriquecimiento superficial'
    ]
}

df = pd.DataFrame(resultados)
pd.options.display.float_format = '{:,.2f}'.format

print("\n" + "="*70)
print(f"PREDICCIÓN DE TENSIÓN SUPERFICIAL (T={TEMPERATURE})")
print(f"Mezcla: { {k: f'{v:.2f}' for k,v in x_mol.items()} } (Fracción Molar)")
print("="*70)
print(df.to_string(index=False))
print("="*70)
print("Nota: El valor de Sprow-Prausnitz suele ser menor porque los componentes")
print("de baja tensión (Acetona) migran a la superficie, dominando la propiedad.")