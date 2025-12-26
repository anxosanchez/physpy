import plotly.graph_objects as go
import numpy as np

def calcular_hansen_mezcla(fracciones_vol, dDs, dPs, dHs):
    """Calcula el vector HSP de la mezcla por promedio volumétrico."""
    dD_mix = sum(f * d for f, d in zip(fracciones_vol, dDs))
    dP_mix = sum(f * p for f, p in zip(fracciones_vol, dPs))
    dH_mix = sum(f * h for f, h in zip(fracciones_vol, dHs))
    return dD_mix, dP_mix, dH_mix

def plot_hansen_3d(solv_names, dDs, dPs, dHs, f_vols, mixture_hsp):
    """Genera una figura 3D interactiva con Plotly."""
    fig = go.Figure()

    # 1. Puntos de componentes puros
    for name, d, p, h, f in zip(solv_names, dDs, dPs, dHs, f_vols):
        fig.add_trace(go.Scatter3d(
            x=[d], y=[p], z=[h],
            mode='markers+text',
            marker=dict(size=6 + f * 20, opacity=0.8),
            text=[name],
            textposition="top center",
            name=name
        ))

    # 2. Punto de la Mezcla (Estrella grande)
    dm, pm, hm = mixture_hsp
    fig.add_trace(go.Scatter3d(
        x=[dm], y=[pm], z=[hm],
        mode='markers',
        marker=dict(size=10, color='red', symbol='diamond', line=dict(width=2, color='white')),
        name='MEZCLA'
    ))

    # 3. Esfera de influencia (Wireframe / Superficie transparente)
    r = 5.0
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 20)
    x = dm + r * np.outer(np.cos(u), np.sin(v))
    y = pm + r * np.outer(np.sin(u), np.sin(v))
    z = hm + r * np.outer(np.ones(np.size(u)), np.cos(v))

    fig.add_trace(go.Surface(
        x=x, y=y, z=z,
        opacity=0.15,
        colorscale='Greys',
        showscale=False,
        name='Esfera de Influencia'
    ))

    # Configuración del Layout (Tema Oscuro Azulado)
    fig.update_layout(
        template="plotly_dark",
        scene=dict(
            xaxis_title='dD (Dispersión)',
            yaxis_title='dP (Polaridad)',
            zaxis_title='dH (Puentes H)',
            bgcolor="rgba(10, 25, 50, 1)"  # Azul oscuro profundo
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        title="Espacio de Solubilidad de Hansen (3D Interactivo)",
        paper_bgcolor="rgba(10, 25, 50, 1)",
        plot_bgcolor="rgba(0,0,0,0)"
    )

    return fig

