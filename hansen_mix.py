import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ==========================================
# 1. CONFIGURACIÓN DE USUARIO (INPUT MOLAR)
# ==========================================
# ¡ATENCIÓN! Estos valores ahora se interpretan como FRACCIÓN MOLAR
composicion_molar = {
    'Acetona':        0.20,
    'Acetato Butilo': 0.30,
    'Butanol':        0.15,
    'Tolueno':        0.35
}

# ==========================================
# 2. VALIDACIÓN (CHECKSUM)
# ==========================================
total_moles = sum(composicion_molar.values())
if abs(total_moles - 1.0) > 1e-4:
    raise ValueError(f"ERROR: La composición molar suma {total_moles:.4f}. Debe sumar 1.00")

# ==========================================
# 3. BASE DE DATOS
# ==========================================
# MW (g/mol), Rho (g/cm3), HSP (MPa^0.5)
db = {
    'Acetona': {
        'MW': 58.08, 'rho': 0.784, 
        'dD': 15.5, 'dP': 10.4, 'dH': 7.0, 'color': 'red'
    },
    'Acetato Butilo': {
        'MW': 116.16, 'rho': 0.882, 
        'dD': 15.8, 'dP': 3.7,  'dH': 6.3, 'color': 'blue'
    },
    'Butanol': {
        'MW': 74.12, 'rho': 0.810, 
        'dD': 16.0, 'dP': 5.7,  'dH': 15.8, 'color': 'green'
    },
    'Tolueno': {
        'MW': 92.14, 'rho': 0.867, 
        'dD': 18.0, 'dP': 1.4,  'dH': 2.0, 'color': 'orange'
    }
}

# ==========================================
# 4. MOTOR DE CÁLCULO (Molar -> Peso -> Volumen)
# ==========================================
def calcular_propiedades_mezcla(comp_molar, database):
    df = pd.DataFrame(index=comp_molar.keys())
    
    # 1. Entrada Molar (x)
    df['x_molar'] = pd.Series(comp_molar)
    
    # 2. Conversión a Peso (w)
    # w_i = (x_i * MW_i) / sum(x_j * MW_j)
    moles_x_mw = []
    for cmp in df.index:
        moles_x_mw.append(df.loc[cmp, 'x_molar'] * database[cmp]['MW'])
    
    df['masa_relativa'] = moles_x_mw
    peso_total = df['masa_relativa'].sum()
    df['w_peso'] = df['masa_relativa'] / peso_total
    
    # 3. Conversión a Volumen (phi)
    # phi_i = (w_i / rho_i) / sum(w_j / rho_j)
    vol_x_rho = []
    for cmp in df.index:
        vol_x_rho.append(df.loc[cmp, 'w_peso'] / database[cmp]['rho'])
        
    df['vol_relativo'] = vol_x_rho
    vol_total = df['vol_relativo'].sum()
    df['phi_vol'] = df['vol_relativo'] / vol_total
    
    # 4. Calcular Vector HSP de Mezcla (Promedio Volumétrico)
    dD_mix = 0
    dP_mix = 0
    dH_mix = 0
    
    for cmp in df.index:
        phi = df.loc[cmp, 'phi_vol']
        dD_mix += phi * database[cmp]['dD']
        dP_mix += phi * database[cmp]['dP']
        dH_mix += phi * database[cmp]['dH']
        
    return df, (dD_mix, dP_mix, dH_mix)

# Ejecutar
df_resultados, hsp_mix = calcular_propiedades_mezcla(composicion_molar, db)

# ==========================================
# 5. REPORTE NUMÉRICO
# ==========================================
print("\n" + "="*60)
print("REPORTE DE CONVERSIÓN Y PROPIEDADES DE MEZCLA")
print("="*60)
# Formato porcentual para mejor lectura
display_df = df_resultados[['x_molar', 'w_peso', 'phi_vol']].copy() * 100
display_df.columns = ['% Molar', '% Peso', '% Volumen']
pd.options.display.float_format = '{:,.2f}'.format
print(display_df)
print("-" * 60)
print(f"HSP MEZCLA FINAL [MPa^0.5]:")
print(f"  Dispersión (D): {hsp_mix[0]:.2f}")
print(f"  Polaridad  (P): {hsp_mix[1]:.2f}")
print(f"  Puentes H  (H): {hsp_mix[2]:.2f}")
print("="*60)

# ==========================================
# 6. VISUALIZACIÓN 3D (ESPACIO HANSEN)
# ==========================================
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# A. Plotear Componentes Puros
for cmp in composicion_molar:
    data = db[cmp]
    # El tamaño del punto será proporcional a su % en Volumen en la mezcla
    size_point = df_resultados.loc[cmp, 'phi_vol'] * 500 + 50 
    
    ax.scatter(data['dD'], data['dP'], data['dH'], 
               c=data['color'], s=size_point, label=cmp, alpha=0.7, edgecolors='k')
    
    # Etiqueta
    ax.text(data['dD'], data['dP'], data['dH']+0.5, cmp, fontsize=9)
    
    # Línea al plano cero (ayuda visual)
    ax.plot([data['dD'], data['dD']], [data['dP'], data['dP']], [0, data['dH']], 
            '--', color='gray', linewidth=0.5, alpha=0.3)

# B. Plotear La MEZCLA (Centro de Gravedad)
dD_m, dP_m, dH_m = hsp_mix
ax.scatter(dD_m, dP_m, dH_m, c='black', s=250, marker='X', label='MEZCLA RESULTANTE', zorder=10)
ax.text(dD_m, dP_m, dH_m+0.5, "  MEZCLA", fontweight='bold')

# C. Dibujar Esfera de Solubilidad (Referencia visual)
# Supongamos un radio de interacción R0 = 5.0 centrado en la mezcla
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
R_sphere = 5.0
x_sph = dD_m + R_sphere * np.cos(u) * np.sin(v)
y_sph = dP_m + R_sphere * np.sin(u) * np.sin(v)
z_sph = dH_m + R_sphere * np.cos(v)
ax.plot_wireframe(x_sph, y_sph, z_sph, color="black", alpha=0.1, label='Esfera Interacción (R=5)')

# D. Configuración de Ejes
ax.set_xlabel('$\delta_D$ (Dispersión)')
ax.set_ylabel('$\delta_P$ (Polaridad)')
ax.set_zlabel('$\delta_H$ (Puentes H)')
ax.set_title(f'Espacio de Hansen (HSP) - Composición Molar Convertida')

# Límites fijos para mejor perspectiva
ax.set_xlim(14, 19)
ax.set_ylim(0, 12)
ax.set_zlim(0, 18)

plt.legend(loc='upper left')
plt.tight_layout()
plt.show()
