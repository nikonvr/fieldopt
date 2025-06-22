import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import numba
import io

# Essayer d'importer scipy pour l'optimisation avancée et l'intégration
try:
    from scipy.optimize import differential_evolution, minimize
    from scipy.integrate import trapezoid
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    # Si scipy n'est pas disponible, on peut définir une fonction trapézoïdale de remplacement
    def trapezoid(y, x):
        return np.trapz(y, x)

# ==============================================================================
# --- FONCTIONS DE CALCUL (NUMBA OPTIMIZED) ---
# ==============================================================================

@numba.njit(cache=True)
def get_layer_properties_numba(n1_r, n2_r, emp_factors_arr, l0):
    """Calcule les propriétés physiques des couches avec des indices purement réels."""
    n1 = n1_r + 0j
    n2 = n2_r + 0j
    n_layers = len(emp_factors_arr)
    indices_complex = np.empty(n_layers, dtype=np.complex128)
    ep_physical_nm = np.empty(n_layers, dtype=np.float64)

    for i in range(n_layers):
        is_odd_layer = (i % 2 == 0)
        n_layer_complex = n1 if is_odd_layer else n2
        n_layer_real = np.real(n_layer_complex)
        
        if n_layer_real <= 0:
            return np.empty(0, dtype=np.complex128), np.empty(0, dtype=np.float64)
        
        qwot_factor = emp_factors_arr[i]
        ep_physical_nm[i] = (qwot_factor * l0) / (4 * n_layer_real)
        indices_complex[i] = n_layer_complex
        
    return indices_complex, ep_physical_nm

def get_layer_properties_from_list(n1_r, n2_r, emp_factors, l0):
    """Prépare les tableaux de propriétés à partir d'une liste de facteurs."""
    if np.size(emp_factors) == 0:
        return np.empty(0, dtype=np.complex128), np.empty(0, dtype=np.float64)

    emp_factors_arr = np.array(emp_factors, dtype=np.float64)
    return get_layer_properties_numba(n1_r, n2_r, emp_factors_arr, l0)

@numba.njit(cache=True)
def _calcul_admittance_recursive_numba(indices_complex, ep_physical_nm, lambda_val, nSub_complex):
    """Calcule l'admittance effective d'un empilement de manière récursive."""
    k_vac = 2 * np.pi / lambda_val
    n_layers = len(indices_complex)
    if n_layers == 0: return nSub_complex
    
    Y_eff = nSub_complex
    for i in range(n_layers):
        n_curr = indices_complex[i]
        d_curr = ep_physical_nm[i]
        phi = k_vac * n_curr * d_curr
        cos_phi, sin_phi = np.cos(phi), np.sin(phi)
        
        numerator = Y_eff * cos_phi + 1j * n_curr * sin_phi
        denominator = cos_phi + 1j * (Y_eff / n_curr) * sin_phi
        
        if abs(denominator) < 1e-12: return 1e30 + 0j
        Y_eff = numerator / denominator
        
    return Y_eff

def calcul_empilement(n1_r, n2_r, nSub_r, l0, emp_factors, l_range, l_step, n_superstrate_real):
    """Calcule les spectres R, T, A pour un empilement donné."""
    l_nm = np.arange(l_range[0], l_range[1] + l_step, l_step)
    indices_complex, ep_physical_nm = get_layer_properties_from_list(n1_r, n2_r, emp_factors, l0)
    R_results, T_results = np.zeros_like(l_nm, dtype=float), np.zeros_like(l_nm, dtype=float)
    
    n_super = np.complex128(n_superstrate_real)
    n_sub = np.complex128(nSub_r)

    for i_l, current_l_nm in enumerate(l_nm):
        if len(indices_complex) == 0:
            r = (n_super - n_sub) / (n_super + n_sub)
            R_results[i_l] = np.abs(r)**2
        else:
            Y_eff = _calcul_admittance_recursive_numba(indices_complex, ep_physical_nm, current_l_nm, n_sub)
            r = (n_super - Y_eff) / (n_super + Y_eff)
            R_results[i_l] = np.abs(r)**2
            
        T_results[i_l] = 1 - R_results[i_l]

    A_results = np.zeros_like(l_nm, dtype=float)
    R_results = np.clip(R_results, 0, 1)
    T_results = np.clip(T_results, 0, 1)

    return {'l': l_nm, 'R': R_results, 'T': T_results, 'A': A_results}, list(ep_physical_nm)


@numba.njit(cache=True)
def _calcul_metrics_et_champ_numba(indices_c1_cn, ep_c1_cn, n1_r, n2_r, nSub_r, l0, lambda_calc, n_super, integral_points):
    """NOUVELLE FONCTION NUMBA CENTRALE pour calculer R et les intégrales de champ."""
    n_layers = len(ep_c1_cn)
    nSub_complex = nSub_r + 0j
    
    # Calcul de R
    Y_eff = _calcul_admittance_recursive_numba(indices_c1_cn, ep_c1_cn, lambda_calc, nSub_complex)
    r = (n_super - Y_eff) / (n_super + Y_eff)
    R = np.abs(r)**2
    t = 1 + r

    # Calcul des champs
    indices_for_efield = indices_c1_cn[::-1]
    ep_for_efield = ep_c1_cn[::-1]
    
    fields = np.zeros((n_layers + 1, 2), dtype=np.complex128)
    fields[0, 0] = t
    fields[0, 1] = n_super * (1 - r)

    k_vac = 2 * np.pi / lambda_calc
    
    for i in range(n_layers):
        n_curr = indices_for_efield[i]
        d_curr = ep_for_efield[i]
        delta = k_vac * n_curr * d_curr
        cos_delta, sin_delta = np.cos(delta), np.sin(delta)
        
        inv_M_i = np.array([
            [cos_delta, (-1j / np.real(n_curr)) * sin_delta],
            [-1j * np.real(n_curr) * sin_delta, cos_delta]
        ], dtype=np.complex128)
        field_at_prev_interface = np.array([[fields[i, 0]], [fields[i, 1]]])
        field_at_next_interface = inv_M_i @ field_at_prev_interface
        fields[i + 1, 0] = field_at_next_interface[0, 0]
        fields[i + 1, 1] = field_at_next_interface[1, 0]

    # Calcul des intégrales et moyennes
    integrals = np.empty(n_layers, dtype=np.float64)
    averages = np.empty(n_layers, dtype=np.float64)
    
    for i in range(n_layers):
        E_interface, H_interface = fields[i, 0], fields[i, 1]
        n_curr, d_curr = indices_for_efield[i], ep_for_efield[i]
        
        Ai = 0.5 * (E_interface + H_interface / n_curr)
        Bi = 0.5 * (E_interface - H_interface / n_curr)

        z_integral = np.linspace(0, d_curr, integral_points)
        e_values = Ai * np.exp(-1j * k_vac * n_curr * z_integral) + Bi * np.exp(1j * k_vac * n_curr * z_integral)
        e2_values = np.abs(e_values)**2
        
        integral_val = np.trapz(e2_values, z_integral)
        integrals[i] = integral_val
        averages[i] = integral_val / d_curr if d_curr > 1e-9 else 0.0

    return R, averages[::-1], integrals[::-1], fields, indices_for_efield, ep_for_efield

def calcul_opt_metrics(n1_r, n2_r, nSub_r, l0, emp_factors_list, n_super, integral_points):
    """Fonction de calcul optimisée utilisant le nouveau cœur Numba."""
    indices_c1_cn, ep_c1_cn = get_layer_properties_from_list(n1_r, n2_r, emp_factors_list, l0)
    if indices_c1_cn.size == 0:
        R0 = np.abs((n_super - nSub_r) / (n_super + nSub_r))**2
        return {'R': R0, 'ratio_average': 0, 'max_avg_1': 0, 'max_avg_2': 0}

    R, averages, _, _, _, _ = _calcul_metrics_et_champ_numba(
        indices_c1_cn, ep_c1_cn, n1_r, n2_r, nSub_r, l0, l0, n_super, integral_points
    )
    
    averages_n1 = [avg for i, avg in enumerate(averages) if i % 2 == 0]
    averages_n2 = [avg for i, avg in enumerate(averages) if i % 2 != 0]

    max_average_n1 = max(averages_n1) if averages_n1 else 0
    max_average_n2 = max(averages_n2) if averages_n2 else 0
    ratio_avg = max_average_n1 / max_average_n2 if max_average_n2 > 1e-9 else np.inf

    return {
        'R': R, 
        'ratio_average': ratio_avg,
        'max_avg_1': max_average_n1, 
        'max_avg_2': max_average_n2
    }
    
def calcul_champ_electrique(n1_r, n2_r, nSub_r, l0, lambda_calc, emp_factors, n_superstrate_real, integral_points):
    """Calcule la distribution du champ électrique en utilisant le nouveau cœur Numba."""
    indices_c1_cn, ep_c1_cn = get_layer_properties_from_list(n1_r, n2_r, emp_factors, l0)
    n_layers = len(ep_c1_cn)
    if n_layers == 0:
        return np.array([0]), np.array([1]), [], [], []

    _, averages, integrals, fields, indices_for_efield, ep_for_efield = _calcul_metrics_et_champ_numba(
        indices_c1_cn, ep_c1_cn, n1_r, n2_r, nSub_r, l0, lambda_calc, n_superstrate_real, integral_points
    )

    # Reconstitution du profil de champ détaillé pour le tracé (en Python)
    z_coords_final, E2_values_final = [], []
    current_z_start = 0
    k_vac = 2 * np.pi / lambda_calc
    
    for i in range(n_layers):
        E_interface, H_interface = fields[i, 0], fields[i, 1]
        n_curr, d_curr = indices_for_efield[i], ep_for_efield[i]
        
        Ai = 0.5 * (E_interface + H_interface / n_curr)
        Bi = 0.5 * (E_interface - H_interface / n_curr)

        z_local = np.linspace(0, d_curr, 50, endpoint=True)
        E_layer_func = lambda z: Ai * np.exp(-1j * k_vac * n_curr * z) + Bi * np.exp(1j * k_vac * n_curr * z)
        
        z_coords_final.extend(current_z_start + z_local)
        E2_values_final.extend(np.abs(E_layer_func(z_local))**2)
        current_z_start += d_curr
        
    return np.array(z_coords_final), np.array(E2_values_final), list(ep_c1_cn), list(integrals), list(averages)


def top_level_objective_function(p, n1_r, n2_r, nSub_r, l0, seuil_int_1, seuil_int_2, alpha, integral_points, n_super):
    """Fonction coût pour l'optimisation (mode pénalité uniquement)."""
    emp_factors_list = list(p)
    
    metrics = calcul_opt_metrics(n1_r, n2_r, nSub_r, l0, emp_factors_list, n_super, integral_points)
    if not metrics: return 1e30

    metric_1, metric_2 = metrics['max_avg_1'], metrics['max_avg_2']
        
    cost_R = 1.0 - metrics['R']
    
    penalty_1 = 0
    if seuil_int_1 > 1e-9:
        violation_1 = metric_1 - seuil_int_1
        penalty_1 = alpha * (max(0, violation_1) / seuil_int_1)**2
    elif metric_1 > seuil_int_1:
        penalty_1 = 1e30 

    penalty_2 = 0
    if seuil_int_2 > 1e-9:
        violation_2 = metric_2 - seuil_int_2
        penalty_2 = alpha * (max(0, violation_2) / seuil_int_2)**2
    elif metric_2 > seuil_int_2:
        penalty_2 = 1e30
    
    return cost_R + penalty_1 + penalty_2

# ==============================================================================
# --- FONCTIONS DE TRACÉ ---
# ==============================================================================
def plot_spectral_data(res, l0_val, show_r, autoscale_y):
    """Trace le spectre de réflectance."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=90)
    if show_r: ax.plot(res['l'], res['R'], label='Réflectance (R)')
    ax.set_xlabel("Longueur d'onde (nm)")
    ax.set_ylabel('Réflectance')
    ax.grid(True, linestyle='--')
    ax.legend()
    if not autoscale_y: ax.set_ylim(-0.05, 1.05)
    
    if len(res['l']) > 1:
        r_at_l0 = np.interp(l0_val, res['l'], res['R'])
        ax.axvline(x=l0_val, color='r', linestyle=':', lw=1)
        ax.text(0.95, 0.95, f'R({l0_val:.0f}nm) = {r_at_l0:.3f}', va='top', ha='right', transform=ax.transAxes, color='red', fontsize=12, bbox=dict(fc='yellow', alpha=0.7))
    
    fig.tight_layout()
    return fig

def plot_stack_visualization(ep, n_super_val, n1_r, n2_r, nSub_r, emp_factors):
    """Visualise le profil d'indice de l'empilement."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 4), dpi=90)
    if np.size(emp_factors) == 0 or not ep:
        ax.text(0.5, 0.5, "Vide", ha='center', va='center')
        return fig
    
    n_layers = len(ep)
    indices_reels = [n1_r if i % 2 == 0 else n2_r for i in range(n_layers)]
    
    z, n = [], []
    current_z = 0
    for i in range(n_layers):
        z.extend([current_z, current_z + ep[i]])
        n.extend([indices_reels[i], indices_reels[i]])
        ax.text(current_z + ep[i] / 2, indices_reels[i] + 0.05, f"{ep[i]:.1f} nm", ha='center', va='bottom', fontsize=7, color='navy')
        current_z += ep[i]

    z_plot = [-50, 0] + z + [current_z, current_z + 50]
    n_plot = [nSub_r, nSub_r] + n + [n_super_val, n_super_val]
    ax.plot(z_plot, n_plot, drawstyle='steps-post')
    
    ax.text(-25, nSub_r, "SUBSTRAT", ha='center', va='center', fontsize=9, alpha=0.7, bbox=dict(fc='white', alpha=0.5, ec='none'))
    ax.text(current_z+25, n_super_val, "SUPERSTRAT", ha='center', va='center', fontsize=9, alpha=0.7, bbox=dict(fc='white', alpha=0.5, ec='none'))

    ax.set_xlabel('Profondeur z (nm)')
    ax.set_ylabel("Indice réel")
    ax.grid(True, linestyle=':')
    fig.tight_layout()
    return fig

def plot_electric_field(z_coords, E2_values, ep_physical, lambda_efield, integrals, averages):
    """Trace le champ électrique avec le substrat à gauche et affiche uniquement les max par matériau."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=90)

    if z_coords is None or z_coords.size == 0 or not ep_physical:
        ax.text(0.5, 0.5, "Erreur de calcul |E|² ou empilement vide", ha='center')
        return fig
    
    total_thickness = sum(ep_physical)
    z_coords_flipped = total_thickness - z_coords
    sort_idx = np.argsort(z_coords_flipped)
    
    ax.plot(z_coords_flipped[sort_idx], E2_values[sort_idx], color='crimson', label=f'$|E/E_0|^2$ à {lambda_efield:.0f} nm', zorder=10)

    averages_n1 = [avg for i, avg in enumerate(averages) if i % 2 == 0]
    averages_n2 = [avg for i, avg in enumerate(averages) if i % 2 != 0]

    idx_max_n1, max_avg_n1 = (-1, 0)
    if averages_n1:
        max_avg_n1 = max(averages_n1)
        idx_max_n1 = averages.index(max_avg_n1)

    idx_max_n2, max_avg_n2 = (-1, 0)
    if averages_n2:
        max_avg_n2 = max(averages_n2)
        idx_max_n2 = averages.index(max_avg_n2)

    ylim = ax.get_ylim()
    y_positions = [ylim[1] * 0.15, ylim[1] * 0.30] 

    current_pos = 0
    for i, thickness in enumerate(ep_physical):
        color = '#D6EAF8' if i % 2 == 0 else '#FEF9E7'
        ax.axvspan(current_pos, current_pos + thickness, color=color, alpha=0.7, ec=None, zorder=0)

        if i == idx_max_n1:
            label = f"Max Mat 1\n<|E|²>={max_avg_n1:.2f}"
            ax.text(current_pos + thickness / 2, y_positions[0], label, ha='center', va='center', fontsize=8, bbox=dict(boxstyle='round,pad=0.3', fc='lightcoral', alpha=0.9, ec='darkred'))

        if i == idx_max_n2:
            label = f"Max Mat 2\n<|E|²>={max_avg_n2:.2f}"
            ax.text(current_pos + thickness / 2, y_positions[1], label, ha='center', va='center', fontsize=8, bbox=dict(boxstyle='round,pad=0.3', fc='lightblue', alpha=0.9, ec='darkblue'))

        current_pos += thickness
    
    xlim = ax.get_xlim()
    ax.text(xlim[0], ylim[1]*0.95, "SUBSTRAT", ha='left', va='top', fontsize=9, alpha=0.7)
    ax.text(xlim[1], ylim[1]*0.95, "SUPERSTRAT", ha='right', va='top', fontsize=9, alpha=0.7)
    
    ax.set_xlabel("Profondeur z (nm) [0 = interface substrat/C1]")
    ax.set_ylabel("$|E/E_0|^2$ (Intensité normalisée)")
    ax.grid(True, linestyle=':')
    ax.legend(loc='upper right')
    fig.tight_layout()
    return fig


def plot_cost_history(cost_history):
    """Trace l'historique de la fonction coût."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=90)
    ax.set_title("Évolution de la Fonction Coût")
    ax.set_xlabel("Itération")
    ax.set_ylabel("Coût")
    if cost_history:
        valid_indices = [i for i, cost in enumerate(cost_history) if np.isfinite(cost) and cost < 1e20]
        if valid_indices:
            costs_to_plot = [cost_history[i] for i in valid_indices]
            ax.plot(range(len(costs_to_plot)), costs_to_plot, marker='.', linestyle='-')
            ax.set_yscale('log')
    ax.grid(True)
    return fig

# ==============================================================================
# --- GESTION DE L'ÉTAT ET DE L'APPLICATION STREAMLIT ---
# ==============================================================================

def initialize_state():
    """Initialise les variables de session si elles n'existent pas."""
    defaults = {
        'n1_r': 2.3, 'n2_r': 1.48, 'nSub_r': 1.52, 'l0': 550.0,
        'l_range_deb': 450.0, 'l_range_fin': 650.0, 'l_step': 1.0,
        'num_layers': 16, 'integral_points': 50,
        'is_optimizing': False, 'stop_optimization_flag': False,
        'cost_history': [], 'opt_iteration': 0,
        'optimization_log': "Prêt.", 'best_factors_so_far': None,
        'best_cost_so_far': float('inf'),
        'opt_results': None, 'show_r': True, 'autoscale_y': False
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value

    if 'efield_lambda' not in st.session_state:
        st.session_state.efield_lambda = float(st.session_state.l0)

    for i in range(st.session_state.num_layers):
        if f'factor_{i}' not in st.session_state:
            st.session_state[f'factor_{i}'] = 1.0

def apply_optimization_results():
    """Applique les nouveaux facteurs issus d'une optimisation."""
    if 'new_factors_from_opt' in st.session_state and st.session_state.new_factors_from_opt is not None:
        new_factors = st.session_state.new_factors_from_opt
        st.session_state.num_layers = len(new_factors)
        for i, factor_val in enumerate(new_factors):
            st.session_state[f'factor_{i}'] = factor_val
        del st.session_state.new_factors_from_opt

def reset_parameters():
    """Réinitialise tous les paramètres aux valeurs par défaut."""
    keys_to_reset = list(st.session_state.keys())
    for key in keys_to_reset:
        del st.session_state[key]
    st.success("Paramètres réinitialisés aux valeurs par défaut.")

def update_factors_from_num_layers():
    """Initialise les nouveaux facteurs lorsque le nombre de couches change."""
    num = st.session_state.num_layers
    old_factors = get_factors_from_state()
    for i in range(num):
        if i >= len(old_factors):
            st.session_state[f'factor_{i}'] = 1.0

def get_factors_from_state():
    """Rassemble la liste des facteurs à partir de l'état de la session."""
    num_layers = st.session_state.get('num_layers', 0)
    return [st.session_state.get(f'factor_{i}', 1.0) for i in range(num_layers)]

def get_excel_download(params, results):
    """Génère un fichier Excel en mémoire pour le téléchargement."""
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        param_data = {
            'Paramètre': list(params.keys()),
            'Valeur': [str(v) for v in params.values()]
        }
        df_params = pd.DataFrame(param_data)
        df_results = pd.DataFrame(results)
        
        df_params.to_excel(writer, sheet_name='Paramètres', index=False)
        df_results.to_excel(writer, sheet_name='Résultats Spectraux', index=False)
    
    return output.getvalue()

# ==============================================================================
# --- INTERFACE PRINCIPALE ---
# ==============================================================================

st.set_page_config(layout="wide", page_title="Simulateur Optique")
st.title("Simulateur Optique Avancé")

initialize_state()
apply_optimization_results()

# --- SIDEBAR DE CONFIGURATION ---
with st.sidebar:
    st.header("Configuration")
    
    st.button("🔄 Réinitialiser Tout", on_click=reset_parameters)

    with st.expander("Matériaux et Spectre", expanded=True):
        st.number_input("Matériau 1 (n1)", 0.1, 10.0, step=0.01, format="%.2f", key='n1_r')
        st.number_input("Matériau 2 (n2)", 0.1, 10.0, step=0.01, format="%.2f", key='n2_r')
        st.number_input("Substrat (indice)", 0.1, 10.0, step=0.01, format="%.2f", key='nSub_r')
        st.slider("Longueur d'onde de centrage λ0 (nm)", 200, 1200, step=1, key='l0')

    with st.expander("Paramètres Spectraux", expanded=True):
        st.number_input("Début intervalle spectral (nm)", 1.0, 2000.0, step=10.0, key='l_range_deb')
        st.number_input("Fin intervalle spectral (nm)", 1.0, 2000.0, step=10.0, key='l_range_fin')
        st.number_input("Pas spectral (nm)", 0.01, 100.0, step=0.1, format="%.2f", key='l_step')

    with st.expander("Épaisseurs des Couches (Facteurs QWOT)", expanded=True):
        st.number_input(
            "Nombre de couches", 1, 30,
            key='num_layers', on_change=update_factors_from_num_layers
        )
        for i in range(st.session_state.num_layers):
            st.number_input(
                f"Facteur Couche {i+1}", 0.0, 10.0, step=0.01, format="%.4g", key=f"factor_{i}"
            )

    # --- Section Optimisation ---
    if SCIPY_AVAILABLE:
        with st.expander("Optimisation", expanded=True):
            st.info("Mode d'optimisation : Coût = (1-R) + Pénalités")
            seuil_int_1 = st.number_input("Seuil Max <|E|²> (Mat. 1)", 0.0, 1000.0, 0.1, 0.01, key="seuil_1")
            seuil_int_2 = st.number_input("Seuil Max <|E|²> (Mat. 2)", 0.0, 1000.0, 10.0, 0.1, key="seuil_2")
            alpha = st.number_input("Poids Pénalité α", 0.0, 10000.0, 1.0, 0.1, key="alpha_penalty")
            st.number_input("Points pour Intégrale |E|²", 10, 500, step=10, key='integral_points')
            
            opt_cols = st.columns(2)
            if opt_cols[0].button("🚀 Lancer Opt. Globale", disabled=st.session_state.is_optimizing, use_container_width=True):
                st.session_state.is_optimizing = True
                st.session_state.stop_optimization_flag = False
                st.session_state.cost_history = []
                st.session_state.opt_iteration = 0
                st.session_state.best_cost_so_far = float('inf')
                st.session_state.best_factors_so_far = None
                st.session_state.optimization_log = "Lancement de l'optimisation globale..."
                st.session_state.opt_method = "global"
                st.rerun()
            
            if opt_cols[1].button("⚡ Lancer Raffinage Local", disabled=st.session_state.is_optimizing, use_container_width=True):
                st.session_state.is_optimizing = True
                st.session_state.stop_optimization_flag = False
                st.session_state.cost_history = []
                st.session_state.opt_iteration = 0
                st.session_state.best_cost_so_far = float('inf')
                st.session_state.best_factors_so_far = None
                st.session_state.optimization_log = "Lancement du raffinage local..."
                st.session_state.opt_method = "local"
                st.rerun()
            
            if st.session_state.is_optimizing:
                if st.button("🛑 Arrêter l'Optimisation", use_container_width=True):
                    st.session_state.stop_optimization_flag = True
                    st.session_state.optimization_log = "Arrêt de l'optimisation demandé..."
            
            log_placeholder = st.empty()
            log_placeholder.info(st.session_state.optimization_log)

    else:
        st.sidebar.warning("Le module 'scipy' est requis pour l'optimisation.")

# --- PANNEAU PRINCIPAL POUR LES GRAPHIQUES ---
tab_spectral, tab_stack, tab_efield, tab_cost = st.tabs([
    "Graphique Spectral", "Visualisation Empilement", "Champ Électrique |E|²", "Fonction Coût"
])

# --- Logique de calcul et d'affichage ---
if st.session_state.is_optimizing:
    s = st.session_state
    
    spectral_placeholder = tab_spectral.empty()
    stack_placeholder = tab_stack.empty()
    efield_placeholder = tab_efield.empty()
    cost_placeholder = tab_cost.empty()

    initial_factors = get_factors_from_state()
    bounds = [(0.1, 2.5)] * len(initial_factors)
    
    args = (
        s.n1_r, s.n2_r, s.nSub_r, s.l0,
        s.seuil_1, s.seuil_2, s.alpha_penalty,
        s.integral_points, 1.0
    )
    
    def _update_all_plots_during_opt(factors):
        """Met à jour tous les graphiques pendant l'optimisation."""
        spectral_results, ep = calcul_empilement(s.n1_r, s.n2_r, s.nSub_r, s.l0, factors, (s.l_range_deb, s.l_range_fin), s.l_step, 1.0)
        z, E2, ep_phys, integrals, averages = calcul_champ_electrique(s.n1_r, s.n2_r, s.nSub_r, s.l0, s.l0, factors, 1.0, s.integral_points)
        spectral_placeholder.pyplot(plot_spectral_data(spectral_results, s.l0, s.show_r, s.autoscale_y), clear_figure=True)
        stack_placeholder.pyplot(plot_stack_visualization(ep, 1.0, s.n1_r, s.n2_r, s.nSub_r, factors), clear_figure=True)
        efield_placeholder.pyplot(plot_electric_field(z, E2, ep_phys, s.l0, integrals, averages), clear_figure=True)
        cost_placeholder.pyplot(plot_cost_history(s.cost_history), clear_figure=True)


    def callback_scipy(xk, convergence=None):
        """Callback pour mettre à jour l'interface pendant l'optimisation."""
        s.opt_iteration += 1
        cost = top_level_objective_function(xk, *args)
        s.cost_history.append(cost)
        
        if cost < s.best_cost_so_far:
            s.best_cost_so_far = cost
            s.best_factors_so_far = list(xk)
        
        log_message = f"Itération {s.opt_iteration}: Coût={cost:.5f} (Meilleur={s.best_cost_so_far:.5f})"
        log_placeholder.info(log_message)
        
        update_interval = 50
        if s.opt_iteration % update_interval == 0 and s.best_factors_so_far:
            _update_all_plots_during_opt(s.best_factors_so_far)
        
        if s.stop_optimization_flag:
            raise StopIteration("Arrêt manuel par l'utilisateur")

    was_stopped = False
    res = None
    try:
        log_placeholder.info(f"Optimisation en cours... Itération {s.opt_iteration}")
        if s.opt_method == "global":
            res = differential_evolution(top_level_objective_function, bounds, args=args, callback=callback_scipy, maxiter=1500, polish=True)
        else:
            res = minimize(top_level_objective_function, initial_factors, args=args, method='L-BFGS-B', bounds=bounds, callback=callback_scipy, options={'maxiter': 500, 'ftol': 1e-9})
    
    except StopIteration:
        was_stopped = True
        s.optimization_log = f"Optimisation arrêtée à l'itération {s.opt_iteration}."
        log_placeholder.warning(s.optimization_log)
    
    except Exception as e:
        s.optimization_log = f"Erreur durant l'optimisation: {e}"
        log_placeholder.error(s.optimization_log)

    final_factors = None
    if res and res.success:
        final_factors = res.x
        s.optimization_log = f"Optimisation terminée ! Coût final : {res.fun:.5f}"
        log_placeholder.success(s.optimization_log)
    elif s.best_factors_so_far is not None:
        final_factors = s.best_factors_so_far
        if was_stopped:
            s.optimization_log = "Optimisation arrêtée. Le meilleur résultat trouvé a été appliqué."
            log_placeholder.warning(s.optimization_log)
        else:
            s.optimization_log = "Optimisation non convergée. Le meilleur résultat trouvé a été appliqué."
            log_placeholder.warning(s.optimization_log)
    else:
        s.optimization_log = "L'optimisation a échoué sans trouver de solution."
        log_placeholder.error(s.optimization_log)

    if final_factors is not None:
        s.new_factors_from_opt = list(final_factors)

    s.is_optimizing = False
    s.stop_optimization_flag = False
    st.rerun()

else:
    s = st.session_state
    current_factors = get_factors_from_state()

    spectral_results, ep = calcul_empilement(
        s.n1_r, s.n2_r, s.nSub_r, s.l0, current_factors, 
        (s.l_range_deb, s.l_range_fin), s.l_step, 1.0
    )
    
    with tab_spectral:
        c1, c2 = st.columns(2)
        c1.checkbox("Afficher Réflectance (R)", key='show_r')
        c2.checkbox("Échelle Y Automatique", key='autoscale_y')
        st.pyplot(plot_spectral_data(spectral_results, s.l0, s.show_r, s.autoscale_y), clear_figure=True)

    with tab_stack:
        st.pyplot(plot_stack_visualization(ep, 1.0, s.n1_r, s.n2_r, s.nSub_r, current_factors), clear_figure=True)

    with tab_efield:
        st.slider(
             "λ pour calcul |E|² (nm)", 
             min_value=float(s.l_range_deb), 
             max_value=float(s.l_range_fin), 
             step=float(s.l_step), 
             key='efield_lambda'
        )
        z, E2, ep_phys, integrals, averages = calcul_champ_electrique(
            s.n1_r, s.n2_r, s.nSub_r, s.l0, s.efield_lambda, 
            current_factors, 1.0, s.integral_points
        )
        st.pyplot(plot_electric_field(z, E2, ep_phys, s.efield_lambda, integrals, averages), clear_figure=True)

    with tab_cost:
        st.pyplot(plot_cost_history(s.cost_history), clear_figure=True)

    st.sidebar.markdown("---")
    st.sidebar.header("Export")
    
    current_params = {
        'n1_r': s.n1_r, 'n2_r': s.n2_r, 'nSub_r': s.nSub_r,
        'l0': s.l0, 'l_range_deb': s.l_range_deb,
        'l_range_fin': s.l_range_fin, 'l_step': s.l_step,
        'num_layers': s.num_layers,
        'factors': ','.join([f"{f:.4f}" for f in current_factors])
    }
    
    excel_data = get_excel_download(current_params, spectral_results)
    
    st.sidebar.download_button(
        label="📥 Télécharger les Résultats (Excel)",
        data=excel_data,
        file_name=f"simu_optique_{datetime.datetime.now():%Y%m%d_%H%M%S}.xlsx",
        mime="application/vnd.ms-excel"
    )
