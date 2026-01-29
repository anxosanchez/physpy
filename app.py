import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- 1. GLOBAL DATABASE ---
DATA = {
    "Solvents (Mandatory)": {
        "Methyl Sunflowerate": {"hsp": [16.2, 3.2, 3.8], "rho": 880, "eacn": 1.5, "fp": 170, "price": 1.65, "ghs": []},
        "Methyl Soyate": {"hsp": [16.1, 3.1, 3.7], "rho": 885, "eacn": 1.4, "fp": 175, "price": 1.60, "ghs": []},
        "Methyl Palmitate": {"hsp": [16.3, 3.3, 3.9], "rho": 870, "eacn": 1.8, "fp": 180, "price": 1.55, "ghs": []},
        "DBE (Dibasic Esters)": {"hsp": [16.5, 7.5, 7.0], "rho": 1060, "eacn": -5.4, "fp": 108, "price": 2.85, "ghs": ["H319"]},
    },
    "Cosolvents (Optional)": {
        "None": {"hsp": [0, 0, 0], "rho": 1.0, "f_hld": 0, "fp": 200, "price": 0, "ghs": []},
        "Propylene Carbonate": {"hsp": [20.0, 18.0, 4.1], "rho": 1200, "f_hld": 0.1, "fp": 132, "price": 2.10, "ghs": ["H319"]},
        "Glycerin": {"hsp": [17.4, 12.1, 29.3], "rho": 1260, "f_hld": -0.8, "fp": 160, "price": 0.95, "ghs": []},
        "Butyl Diglicol": {"hsp": [16.0, 7.0, 10.6], "rho": 953, "f_hld": -0.2, "fp": 105, "price": 3.20, "ghs": ["H319"]},
    },
    "Surfactants": {
        "APG (Non-Ionic)": {"cc": 1.5, "hlb": 13.5, "rho": 1100, "hsp": [18.0, 12.0, 15.0], "price": 1.90, "ghs": ["H318"]},
        "SLES (Anionic)": {"cc": -2.0, "hlb": 40.0, "rho": 1050, "hsp": [17.5, 11.0, 9.5], "price": 1.85, "ghs": ["H315", "H318"]},
    },
    "Resins": {
        "Alquidic": {"hsp": [18.5, 4.5, 5.1], "r0": 8.0},
        "Nitrocellulose": {"hsp": [15.4, 10.1, 8.8], "r0": 11.5},
        "Polyurethane": {"hsp": [17.8, 10.5, 11.2], "r0": 9.0},
        "PVC": {"hsp": [18.8, 9.2, 3.5], "r0": 7.5},
    }
}

# --- 2. LOGIC ---
st.set_page_config(page_title="MicroSaaS Pro", layout="wide")
st.title("🧪 Industrial Microemulsion Manufacturing Suite")

with st.sidebar:
    st.header("Formulation Parameters")
    res_k = st.selectbox("Target Resin", list(DATA["Resins"].keys()))
    sol_k = st.selectbox("Solvent", list(DATA["Solvents (Mandatory)"].keys()))
    cos_k = st.selectbox("Cosolvent", list(DATA["Cosolvents (Optional)"].keys()))
    sur_k = st.selectbox("Surfactant", list(DATA["Surfactants"].keys()))
    
    p_sol = st.slider("% Solvent", 5, 50, 30)
    p_cos = st.slider("% Cosolvent", 0, 20, 5)
    p_sur = st.slider("% Surfactant", 5, 25, 15)
    p_alc = st.slider("% Cosurfactant (Alcohol)", 2, 15, 10)
    p_wat = 100 - p_sol - p_cos - p_sur - p_alc
    st.caption(f"Water: {p_wat}%")
    salinity = st.number_input("Salinity (% NaCl in mix)", 0.0, 5.0, 0.5)

# Calculations
S, C, Sf, R = DATA["Solvents (Mandatory)"][sol_k], DATA["Cosolvents (Optional)"][cos_k], DATA["Surfactants"][sur_k], DATA["Resins"][res_k]
ws = np.array([p_sol, p_cos, p_sur, p_alc, p_wat]) / 100
rhos = np.array([S['rho'], C['rho'], Sf['rho'], 810, 997])
rho_mix = 1 / np.sum(ws / rhos)

# --- 3. DASHBOARD TABS ---
t1, t2, t3 = st.tabs(["📊 Lab Analysis", "🛡️ Regulatory", "🏭 Batch Manufacturing"])

with t1:
    # (Existing RED/HLD/Hansen 3D Logic goes here)
    st.subheader("Performance Indicators")
    col1, col2 = st.columns(2)
    # Calculation of RED and HLD
    h_mix = np.dot((ws/rhos)/np.sum(ws/rhos), np.array([S['hsp'], C['hsp'], Sf['hsp'], [16,6,16], [15.5,16,42]]))
    red = np.sqrt(4*(h_mix[0]-R['hsp'][0])**2 + (h_mix[1]-R['hsp'][1])**2 + (h_mix[2]-R['hsp'][2])**2) / R['r0']
    hld = np.log(salinity + 0.001) - 0.17*S['eacn'] + Sf['cc']
    
    col1.metric("Solubility RED", f"{red:.2f}")
    col2.metric("Stability HLD", f"{hld:.2f}")

with t2:
    st.subheader("GHS Safety Phrases")
    st.info("Based on component concentrations > 10%")
    # Classification logic...

with t3:
    st.header("Batch Weight Calculator")
    batch_liters = st.number_input("Target Batch Size (Liters)", 10, 10000, 1000)
    brine_conc = st.slider("Brine Stock Concentration (% NaCl)", 5.0, 20.0, 10.0)
    
    total_kg = batch_liters * (rho_mix / 1000)
    
    # Calculate weights
    w_sol = total_kg * (p_sol/100)
    w_cos = total_kg * (p_cos/100)
    w_sur = total_kg * (p_sur/100)
    w_alc = total_kg * (p_alc/100)
    
    # Brine calculation
    # target_salt = total_kg * (salinity/100)
    # w_brine = target_salt / (brine_conc/100)
    w_brine = (total_kg * (salinity/100)) / (brine_conc/100)
    w_water = (total_kg * (p_wat/100)) - w_brine
    
    st.success(f"Total Mass to weigh: {total_kg:.2f} kg")
    
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("### 🛢️ Organic & Surfactant Phase")
        st.write(f"- **{S['nome'] if 'nome' in S else sol_k}:** {w_sol:.2f} kg")
        st.write(f"- **{C['nome'] if 'nome' in C else cos_k}:** {w_cos:.2f} kg")
        st.write(f"- **{Sf['nome'] if 'nome' in Sf else sur_k}:** {w_sur:.2f} kg")
        st.write(f"- **Cosurfactant:** {w_alc:.2f} kg")
    
    with col_b:
        st.markdown("### 💧 Aqueous Phase")
        st.write(f"- **NaCl Brine ({brine_conc}%):** {w_brine:.2f} kg")
        st.write(f"- **Pure Water:** {w_water:.2f} kg")
