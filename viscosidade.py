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
composicion_peso = {
    'Acetona':        0.20,
    'Acetato Butilo': 0.30,
    'Butanol':        0.15,  # 1-Butanol
    'Tolueno':        0.35
}

TEMPERATURE = Q_(25, 'degC') # La viscosidad es MUY sensible a T
# ==========================================

# --- VALIDACIÓN (CHECKSUM) ---
total_peso = sum(composicion_peso.values())
tolerancia = 1e-6

if abs(total_peso - 1.0) > tolerancia:
    raise ValueError(f"CRITICAL ERROR: La composición suma {total_peso:.4f} en lugar de 1.00.")
else:
    print(f"✓ Checksum OK: La composición suma {total_peso*100:.0f}%")

# ==========================================
# BASE DE DATOS Y PARÁMETROS
# ==========================================
# Datos de Viscosidad a 25°C (en cP o mPa.s)
db = {
    'Acetona':        {'MW': 58.08,  'mu': Q_(0.306, 'cP')},
    'Acetato Butilo': {'MW': 116.16, 'mu': Q_(0.734, 'cP')},
    'Butanol':        {'MW': 74.12,  'mu': Q_(2.544, 'cP')}, # Alto por puentes de H
    'Tolueno':        {'MW': 92.14,  'mu': Q_(0.560, 'cP')}
}

# MATRIZ DE INTERACCIÓN BINARIA (G_ij) PARA GRUNBERG-NISSAN
# Valor 0 = Ideal
# Valor Negativo = Ruptura de estructura (Viscosidad cae más de lo esperado)
# Valor Positivo = Formación de complejos (Viscosidad sube)
# Simulamos que Tolueno y Acetona "cortan" la viscosidad del alcohol.
G_ij = pd.DataFrame(0.0, index=db.keys(), columns=db.keys())

# Asignación de parámetros de interacción estimados (Simulación didáctica)
# Interacción Butanol-Tolueno (Fuerte ruptura de puentes H)
G_ij.loc['Butanol', 'Tolueno'] = -0.4 
G_ij.loc['Tolueno', 'Butanol'] = -0.4

# Interacción Acetona-Butanol
G_ij.loc['Acetona', 'Butanol'] = -0.2
G_ij.loc['Butanol', 'Acetona'] = -0.2

# ==========================================
# FUNCIONES AUXILIARES
# ==========================================

def get_molar_fractions(comp_weight, database):
    """Convierte fracción peso a fracción molar"""
    moles = {}
    total_moles = 0
    for cmp, w in comp_weight.items():
        mw = database[cmp]['MW']
        n = w / mw
        moles[cmp] = n
        total_moles += n
    return {k: v/total_moles for k, v in moles.items()}

# ==========================================
# MÉTODOS DE CÁLCULO DE VISCOSIDAD
# ==========================================

# 1. LINEAL (Incorrecto, solo referencia)
def calc_mu_linear(x, database):
    mu_mix = 0
    for c, xi in x.items():
        mu_mix += xi * database[c]['mu'].magnitude
    return Q_(mu_mix, 'cP')

# 2. ARRHENIUS (Regla Logarítmica - Estándar Ideal)
def calc_mu_arrhenius(x, database):
    # ln(mu_mix) = sum(xi * ln(mu_i))
    ln_mu_mix = 0
    for c, xi in x.items():
        ln_mu_mix += xi * np.log(database[c]['mu'].magnitude)
    
    mu_mix = np.exp(ln_mu_mix)
    return Q_(mu_mix, 'cP')

# 3. KENDALL-MONROE (Raíz Cúbica)
def calc_mu_kendall(x, database):
    # mu_mix^(1/3) = sum(xi * mu_i^(1/3))
    sum_root = 0
    for c, xi in x.items():
        val = database[c]['mu'].magnitude
        sum_root += xi * (val**(1/3))
    
    mu_mix = sum_root**3
    return Q_(mu_mix, 'cP')

# 4. GRUNBERG-NISSAN (Semi-empírico con interacciones)
def calc_mu_grunberg(x, database, interactions):
    # ln(mu) = sum(xi*ln(mui)) + sum(sum(xi*xj*Gij))
    
    # Parte ideal (Arrhenius)
    ln_ideal = 0
    for c, xi in x.items():
        ln_ideal += xi * np.log(database[c]['mu'].magnitude)
    
    # Parte de exceso (Interacciones)
    excess_term = 0
    comps = list(x.keys())
    for i in comps:
        for j in comps:
            if i != j:
                # G_ij * xi * xj
                g_val = interactions.loc[i, j]
                excess_term += g_val * x[i] * x[j]
    
    ln_mu_mix = ln_ideal + excess_term # El exceso suele ser negativo aquí
    return Q_(np.exp(ln_mu_mix), 'cP')

# ==========================================
# EJECUCIÓN
# ==========================================

# 1. Calcular molares
x_molar = get_molar_fractions(composicion_peso, db)

# 2. Calcular viscosidades
mu_linear = calc_mu_linear(x_molar, db)
mu_arrh = calc_mu_arrhenius(x_molar, db)
mu_kendall = calc_mu_kendall(x_molar, db)
mu_grunberg = calc_mu_grunberg(x_molar, db, G_ij)

# 3. Generar Tabla
resultados = {
    'Método de Predicción': [
        'Promedio Lineal (Referencia)', 
        'Arrhenius (Ideal Log)', 
        'Kendall-Monroe (Raíz Cúbica)', 
        'Grunberg-Nissan (Interacción)'
    ],
    'Viscosidad [cP]': [
        mu_linear.magnitude,
        mu_arrh.magnitude,
        mu_kendall.magnitude,
        mu_grunberg.magnitude
    ],
    'Comentario': [
        'Sobrestima (Inexacto)',
        'Estándar básico',
        'Mejor para hidrocarburos',
        'Considera ruptura de puentes H'
    ]
}

df = pd.DataFrame(resultados)
# Calcular diferencia respecto a Arrhenius (referencia base)
df['Diferencia %'] = ((df['Viscosidad [cP]'] - mu_arrh.magnitude) / mu_arrh.magnitude) * 100

pd.options.display.float_format = '{:,.3f}'.format

print("\n" + "="*60)
print(f"PREDICCIÓN DE VISCOSIDAD DE MEZCLA (T={TEMPERATURE})")
print(f"Composición Molar: { {k: f'{v:.2f}' for k,v in x_molar.items()} }")
print("="*60)
print(df.to_string(index=False))
print("="*60)
print("\nNota: La alta viscosidad del Butanol puro (2.54 cP) es reducida")
print("drásticamente por la Acetona y el Tolueno.")
