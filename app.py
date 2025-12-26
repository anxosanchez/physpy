# app.py
import streamlit as st
import pandas as pd
import numpy as np
from data.database import SOLVENTS
from calculations import densidad, viscosidad, tension_superficial, hansen_mix

# Configuración de página con estado de sidebar
st.set_page_config(
    page_title="SolventMix AI | Next-Gen Simulator",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- SISTEMA DE DISEÑO: NEXT-GEN GLASSMORPHISM ---
st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&display=swap');

    /* Variables Globales */
    :root {
        --primary-gradient: linear-gradient(135deg, #00d2ff 0%, #3a7bd5 100%);
        --accent-gradient: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        --glass-bg: rgba(255, 255, 255, 0.03);
        --glass-border: rgba(255, 255, 255, 0.1);
        --text-main: #e0e5ec;
        --text-dim: #a0aec0;
    }

    /* Reset global y Tipografía */
    .stApp {
        background: radial-gradient(circle at top right, #1a1c2c, #0d0e12);
        font-family: 'Outfit', sans-serif;
        color: var(--text-main);
    }

    /* Sidebar Estilizado */
    [data-testid="stSidebar"] {
        background-color: rgba(15, 17, 25, 0.95);
        border-right: 1px solid var(--glass-border);
    }
    
    /* Títulos Impactantes */
    h1 {
        background: var(--primary-gradient);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-weight: 700 !important;
        font-size: 3rem !important;
        letter-spacing: -1px;
        margin-bottom: 0px !important;
    }
    h3 {
        color: var(--text-dim) !important;
        font-weight: 400 !important;
        margin-top: -10px !important;
    }

    /* Tarjetas de Métricas Custom (CSS puro para evitar el estilo "cutre") */
    .metric-card {
        background: var(--glass-bg);
        border: 1px solid var(--glass-border);
        backdrop-filter: blur(12px);
        padding: 1.5rem;
        border-radius: 20px;
        text-align: center;
        transition: all 0.3s ease;
        box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.37);
    }
    .metric-card:hover {
        transform: translateY(-5px);
        border-color: rgba(0, 210, 255, 0.4);
        box-shadow: 0 8px 32px 0 rgba(0, 210, 255, 0.1);
    }
    .metric-label {
        font-size: 0.9rem;
        color: var(--text-dim);
        text-transform: uppercase;
        letter-spacing: 2px;
        margin-bottom: 0.5rem;
    }
    .metric-value {
        font-size: 2.2rem;
        font-weight: 700;
        background: var(--primary-gradient);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    .metric-unit {
        font-size: 0.8rem;
        color: var(--text-dim);
        margin-left: 5px;
    }

    /* Botón de Acción Principal */
    .stButton>button {
        background: var(--primary-gradient);
        color: white;
        border: none;
        padding: 15px 30px;
        border-radius: 12px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 1px;
        box-shadow: 0 4px 15px rgba(0, 210, 255, 0.3);
        transition: all 0.3s ease;
        width: 100%;
    }
    .stButton>button:hover {
        box-shadow: 0 6px 20px rgba(0, 210, 255, 0.5);
        transform: scale(1.02);
    }

    /* Tablas Modernas */
    .stTable {
        background: var(--glass-bg);
        border-radius: 15px !important;
        overflow: hidden;
    }

    /* Personalización de Inputs */
    .stNumberInput, .stMultiSelect {
        background: var(--glass-bg) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 8px !important;
    }
    
    /* Expander */
    .streamlit-expanderHeader {
        background-color: transparent !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 12px !important;
        color: #00d2ff !important;
    }
    </style>
""", unsafe_allow_html=True)

# --- CABECERA ---
st.title("SOLVENT_MIX AI")
st.subheader("Molecular Engine for Industrial Formulations")
st.markdown("<br>", unsafe_allow_html=True)

# --- SIDEBAR: COMPOSITOR ---
with st.sidebar:
    st.markdown("### 🧬 COMPOSITION")
    base_calculo = st.radio("Unit Base", ["Volume Fraction", "Mass Fraction", "Molar Fraction"])
    
    selected_solvents = st.multiselect(
        "Industrial Components",
        options=list(SOLVENTS.keys()),
        default=["Agua", "Etanol"]
    )

    comps = {}
    if selected_solvents:
        st.markdown("---")
        for s in selected_solvents:
            comps[s] = st.number_input(f"{s} (%)", min_value=0.0, max_value=100.0, value=100.0/len(selected_solvents))

        total = sum(comps.values())
        progress_color = "#00d2ff" if abs(total - 100.0) < 0.01 else "#ff4b4b"
        st.markdown(f"**Total: <span style='color:{progress_color}; font-size:1.2rem'>{total:.2f}%</span>**", unsafe_allow_html=True)
        
        if abs(total - 100.0) > 0.01:
            st.warning("⚠️ Sum must be exactly 100%")
            st.stop()
        else:
            st.success("✅ Ready to simulate")

# --- ÁREA PRINCIPAL ---
col_config, col_spacer, col_res = st.columns([1, 0.1, 2.2])

with col_config:
    st.markdown("### 🛠️ ENGINE CONFIG")
    with st.expander("Physics Models", expanded=True):
        m_dens = st.selectbox("Density Model", ["Modified Rackett", "COSTALD (HBT)", "PR-Peneloux", "Ideal Mixing"], index=0)
        m_visc = st.selectbox("Viscosity Model", ["Arrhenius", "Grunberg-Nissan (Interaction)", "Kendall-Monroe", "Linear Average"], index=0)
        m_tens = st.selectbox("Surface Tension", ["Macleod-Sugden", "Sprow-Prausnitz (Ideal Surface)", "Linear Volumetric", "Linear Molar"], index=0)
        temp = st.slider("Process Temp (°C)", 10.0, 80.0, 25.0)
    
    st.markdown("<br>", unsafe_allow_html=True)
    calculate = st.button("🚀 RUN SIMULATION")

# --- LÓGICA DE CÁLCULO (Se activa al presionar el botón) ---
if calculate:
    T_kelvin = temp + 273.15
    f_input = np.array([v/100 for v in comps.values()])
    names = list(comps.keys())
    props = [SOLVENTS[n] for n in names]
    
    # Conversiones
    mws = np.array([p['MW'] for p in props])
    rhos_ref = np.array([p['rho_ref'] for p in props])
    
    if base_calculo == "Mass Fraction":
        f_mass = f_input
        moles_rel = f_mass / mws
        f_mol = moles_rel / np.sum(moles_rel)
        vols_rel = f_mass / rhos_ref
        f_vol = vols_rel / np.sum(vols_rel)
    elif base_calculo == "Molar Fraction":
        f_mol = f_input
        mass_rel = f_mol * mws
        f_mass = mass_rel / np.sum(mass_rel)
        vols_rel = mass_rel / rhos_ref
        f_vol = vols_rel / np.sum(vols_rel)
    else: # Volume
        f_vol = f_input
        mass_rel = f_vol * rhos_ref
        f_mass = mass_rel / np.sum(mass_rel)
        moles_rel = f_mass / mws
        f_mol = moles_rel / np.sum(moles_rel)

    try:
        # Cálculo de Propiedades
        if "Rackett" in m_dens:
            rho_mix = densidad.modified_rackett(T_kelvin, f_mol, props)
        elif "COSTALD" in m_dens:
            rho_mix = densidad.costald_density(T_kelvin, f_mol, props)
        elif "PR-Peneloux" in m_dens:
            rho_mix = densidad.pr_peneloux_density(T_kelvin, f_mol, props)
        else: # Ideal Mixing (usando fracciones en peso como en densidade.py)
            rho_mix = densidad.densidad_ideal(f_mass, [p['rho_ref']/1000 for p in props]) * 1000
        
        mu_puras = [p['visc_ref'] for p in props]
        if "Arrhenius" in m_visc: 
            mu_mix = viscosidad.metodo_arrhenius(f_mol, mu_puras)
        elif "Nissan" in m_visc: 
            mu_mix = viscosidad.metodo_grunberg_nissan(f_mol, mu_puras, names)
            st.info("💡 Grunberg-Nissan: Se aplican parámetros de interacción (Gij) para simular la ruptura de puentes de hidrógeno.")
        elif "Kendall" in m_visc:
            mu_mix = viscosidad.metodo_kendall_monroe(f_mol, mu_puras)
        else:
            mu_mix = viscosidad.metodo_lineal(f_mol, mu_puras)
            
        sig_puras = [p['sigma_ref'] for p in props]
        if "Macleod" in m_tens:
            paracors = [p['Paracor'] for p in props]
            tens_mix = tension_superficial.macleod_sugden(f_mol, paracors, mws, rho_mix)
        elif "Sprow" in m_tens:
            rhos_kgm3 = [p['rho_ref'] for p in props]
            tens_mix = tension_superficial.sprow_prausnitz(f_mol, sig_puras, mws, rhos_kgm3, T_kelvin)
            st.info("💡 Sprow-Prausnitz: El modelo considera que los componentes con menor tensión migran a la superficie.")
        elif "Volumetric" in m_tens: 
            tens_mix = tension_superficial.regla_lineal_volumetrica(f_vol, sig_puras)
        else: 
            tens_mix = tension_superficial.regla_lineal_molar(f_mol, sig_puras)
        
        h_d = [p['dD'] for p in props]; h_p = [p['dP'] for p in props]; h_h = [p['dH'] for p in props]
        h_mix = hansen_mix.calcular_hansen_mezcla(f_vol, h_d, h_p, h_h)

        # --- RESULTADOS ---
        with col_res:
            st.markdown("### 📊 MIXTURE PROFILE")
            
            # Cards de Métricas con HTML/CSS
            m1, m2, m3 = st.columns(3)
            with m1:
                val = f"{rho_mix:.2f}" if not np.isnan(rho_mix) else "N/A"
                st.markdown(f"""<div class='metric-card'><div class='metric-label'>Density</div><div class='metric-value'>{val}<span class='metric-unit'>kg/m³</span></div></div>""", unsafe_allow_html=True)
            with m2:
                val = f"{mu_mix:.3f}" if not np.isnan(mu_mix) else "N/A"
                st.markdown(f"""<div class='metric-card'><div class='metric-label'>Viscosity</div><div class='metric-value'>{val}<span class='metric-unit'>cP</span></div></div>""", unsafe_allow_html=True)
            with m3:
                val = f"{tens_mix:.2f}" if not np.isnan(tens_mix) else "N/A"
                st.markdown(f"""<div class='metric-card'><div class='metric-label'>Tension</div><div class='metric-value'>{val}<span class='metric-unit'>mN/m</span></div></div>""", unsafe_allow_html=True)

            st.markdown("<br>", unsafe_allow_html=True)
            
            with st.expander("🔍 Detailed Composition Data", expanded=False):
                df_comp = pd.DataFrame({
                    "Component": names,
                    "Molar (x)": f_mol,
                    "Mass (w)": f_mass,
                    "Vol (φ)": f_vol
                })
                st.table(df_comp.style.format({c: "{:.4f}" for c in df_comp.columns if c != "Component"}))

            st.markdown("---")
            fig = hansen_mix.plot_hansen_3d(names, h_d, h_p, h_h, f_vol, h_mix)
            st.plotly_chart(fig, use_container_width=True)

            # --- EXPORTAR INFORME ---
            st.markdown("<br>", unsafe_allow_html=True)
            report_data = {
                "Parameter": ["Temperature", "Density Model", "Viscosity Model", "Surface Tension Model", 
                              "Result: Density (kg/m3)", "Result: Viscosity (cP)", "Result: Surface Tension (mN/m)",
                              "HSP dD", "HSP dP", "HSP dH"],
                "Value": [f"{temp} °C", m_dens, m_visc, m_tens, 
                          f"{rho_mix:.4f}", f"{mu_mix:.4f}", f"{tens_mix:.4f}",
                          f"{h_mix[0]:.2f}", f"{h_mix[1]:.2f}", f"{h_mix[2]:.2f}"]
            }
            df_report = pd.DataFrame(report_data)
            
            # Combinar con datos de composición
            df_comp_report = df_comp.copy()
            df_comp_report.columns = ["Parameter", "Molar (x)", "Mass (w)", "Vol (phi)"]
            
            # Generar CSV string
            csv = df_report.to_csv(index=False).encode('utf-8')
            
            st.download_button(
                label="📥 DOWNLOAD TECHNICAL REPORT (CSV)",
                data=csv,
                file_name=f"solvent_mix_report_{temp}C.csv",
                mime="text/csv",
            )
            
    except Exception as e:
        st.error(f"Engine Stalled: {e}")
        import traceback
        st.code(traceback.format_exc())
