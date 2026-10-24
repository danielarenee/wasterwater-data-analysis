"""
Script de EDA — Metales/Elementos (ICP-MS/ICP-OES, Actlabs)
tenemos 6 fuentes (L1-L5 lavanderías + L6 El Asturcón) con ST/T, y 1 blanco de campo (BLK)
"""
#%%

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.sans-serif'] = 'Helvetica, Arial, DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#cccccc'
plt.rcParams['axes.linewidth'] = 0.8


# CARGA DE DATOS

def cargar_datos(ruta_csv='outputs/metales_actlabs.csv'):
    """
    separa metadatos de analitos, códigos Actlabs y genera matrices de
    concentraciones numéricas y estados de detección (metales_actlabs.csv
    ya trae LD y unidad por analito, a diferencia del universo orgánico)
    """
    df_full = pd.read_csv(ruta_csv, encoding='utf-8-sig')

    # fila 0 = trazabilidad (código Actlabs); el resto son analitos
    fila_codigos = df_full.iloc[0]
    df_raw = df_full.iloc[1:].reset_index(drop=True).copy()

    df_raw['id'] = pd.to_numeric(df_raw['id'], errors='coerce')
    df_raw['LD'] = pd.to_numeric(df_raw['LD'], errors='coerce')
    df_raw['analito'] = df_raw['analito'].astype(str).str.strip()

    sample_ids = [c for c in df_full.columns if c not in ('id', 'analito', 'unidad', 'LD', 'metodo')]
    codigo_actlabs_map = fila_codigos[sample_ids].to_dict()

    # metadatos de muestras estructurados
    df_muestras_info = pd.DataFrame({'muestra_id': sample_ids})
    df_muestras_info['codigo_actlabs'] = df_muestras_info['muestra_id'].map(codigo_actlabs_map)

    def extraer_lavanderia(x):
        if x == 'BLK':
            return 'BLK'
        return re.match(r'(L\d+)', x).group(1)

    def extraer_tipo(x):
        if x == 'BLK':
            return 'BLK'  # Blanco de campo
        elif '_ST' in x:
            return 'ST'   # Sin tratar
        elif '_T' in x:
            return 'T'    # Tratada
        return 'Otro'

    df_muestras_info['lavanderia'] = df_muestras_info['muestra_id'].apply(extraer_lavanderia)
    df_muestras_info['tipo'] = df_muestras_info['muestra_id'].apply(extraer_tipo)
    df_muestras_info['tipo_desc'] = df_muestras_info['tipo'].map({
        'ST': 'Agua Residual Sin Tratar',
        'T': 'Agua Tratada',
        'BLK': 'Blanco de Campo'
    })

    # Matriz de detección booleana (True = Detectado/Cuantificado, False = Censurado <LD)
    # Nota: ">100000" (saturación ICP-MS, solo ocurre en Na) se trata como detectado,
    # usando 100000 como estimación de concentración.
    df_detected = pd.DataFrame(index=df_raw.index)
    df_concentraciones = pd.DataFrame(index=df_raw.index)

    for col in sample_ids:
        # Na_OES viene con celdas vacías (NaN) en las muestras donde no se midió
        # (ICP-OES solo se reporta cuando ICP-MS satura); "no medido" debe quedar
        # como NA, no confundirse con "detectado" (no contiene '<') ni con "censurado".
        is_missing = df_raw[col].isna()
        raw_vals = df_raw[col].astype(str).str.strip().str.replace(' ', '')
        is_censored = raw_vals.str.contains('<') & ~is_missing

        df_detected[col] = pd.array(np.where(is_missing, pd.NA, ~is_censored), dtype='boolean')

        numeric_vals = pd.to_numeric(raw_vals.str.replace('<', '').str.replace('>', ''), errors='coerce')
        # LD/sqrt(2) (método EPA, Hornung & Reed 1990) — mismo método que imputacion_inorg.py
        # y que imputacion_org.py para el universo orgánico; ver docs/hallazgos_generales.md
        imputed_vals = np.where(is_censored, df_raw['LD'] / np.sqrt(2.0), numeric_vals)
        df_concentraciones[col] = imputed_vals

    return df_raw, df_muestras_info, df_detected, df_concentraciones, sample_ids


# 2. FRECUENCIA DE DETECCIÓN

def calcular_resumen_deteccion(df_raw, df_detected, df_muestras_info, sample_ids):
    """
    por tipo de matriz (Agua Cruda ST, Tratada T); BLK se reporta aparte como
    control de calidad (blanco de campo), no se mezcla con las estadísticas
    de las muestras de proceso.
    """
    columnas_st = df_muestras_info[df_muestras_info['tipo'] == 'ST']['muestra_id'].tolist()
    columnas_t = df_muestras_info[df_muestras_info['tipo'] == 'T']['muestra_id'].tolist()
    columnas_proceso = columnas_st + columnas_t  # excluye BLK

    resumen = pd.DataFrame({
        'id': df_raw['id'],
        'analito': df_raw['analito'],
        'unidad': df_raw['unidad'],
        'LD': df_raw['LD'],
        'metodo': df_raw['metodo'],
        # denominadores dinámicos (notna): la mayoría de analitos se miden en todas las
        # muestras, pero Na_OES solo se reporta donde ICP-MS satura (celdas NA el resto).
        'det_ST': df_detected[columnas_st].sum(axis=1),
        'total_muestras_ST': df_detected[columnas_st].notna().sum(axis=1),
        'pct_det_ST': (df_detected[columnas_st].mean(axis=1).astype(float) * 100).round(1),
        'det_T': df_detected[columnas_t].sum(axis=1),
        'total_muestras_T': df_detected[columnas_t].notna().sum(axis=1),
        'pct_det_T': (df_detected[columnas_t].mean(axis=1).astype(float) * 100).round(1),
        'det_proceso': df_detected[columnas_proceso].sum(axis=1),
        'total_muestras_proceso': df_detected[columnas_proceso].notna().sum(axis=1),
        'pct_det_proceso': (df_detected[columnas_proceso].mean(axis=1).astype(float) * 100).round(1),
        'detectado_en_blk': df_detected['BLK'].fillna(False) if 'BLK' in df_detected.columns else False,
    })

    return resumen


# 3. GENERACIÓN DE VISUALIZACIONES
# Sin título embebido (va en el caption del reporte) y a dpi=150 para mantener
# archivos ligeros; ver docs/contexto_figuras_inorg.txt para título y caption sugeridos.

def graficar_tasa_deteccion(resumen, output_dir='figuras/inorganico'):
    """
    barras horizontales con la frecuencia de detección de cada analito en ST;
    panel único (no hay una familia química equivalente a PFAS/Restringidas
    para los metales, forzar un split no aporta).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Na_OES se excluye: es el mismo elemento que Na (ICP-OES solo reporta cuando
    # ICP-MS satura), su denominador de muestras es distinto y duplicaría a Na.
    resumen = resumen[resumen['analito'] != 'Na_OES']
    detectados = resumen[resumen['det_ST'] > 0].sort_values(by='pct_det_ST', ascending=True)
    n_muestras_st_max = int(resumen['total_muestras_ST'].max())

    fig, ax = plt.subplots(figsize=(9, 13))
    colores = ['#b30000' if blk else '#08519c' for blk in detectados['detectado_en_blk']]
    bars = ax.barh(np.arange(len(detectados)), detectados['pct_det_ST'], color=colores, alpha=0.85, height=0.65)

    ax.set_yticks(np.arange(len(detectados)))
    ax.set_yticklabels(detectados['analito'], fontsize=9)
    ax.set_xlabel(f'% Muestras de Agua Cruda con Detección (> LD) — n={n_muestras_st_max}', fontsize=11, fontweight='bold')
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.tick_params(axis='x', labelsize=10)
    ax.grid(axis='x', linestyle='--', alpha=0.6)

    for bar, val, n_muestras in zip(bars, detectados['pct_det_ST'], detectados['total_muestras_ST']):
        width = bar.get_width()
        n_abs = int(round(width * n_muestras / 100))
        ax.text(width + 1.5, bar.get_y() + bar.get_height() / 2, f'{val:.0f}% ({n_abs}/{n_muestras})',
                ha='left', va='center', fontsize=8, color='#222222')

    from matplotlib.patches import Patch
    leyenda = [Patch(facecolor='#08519c', label='No detectado en BLK'),
               Patch(facecolor='#b30000', label='También detectado en BLK (QC)')]
    ax.legend(handles=leyenda, loc='lower right', fontsize=9, framealpha=0.9)

    plt.tight_layout(pad=1.5)
    ruta = os.path.join(output_dir, '01_tasa_deteccion_metales.png')
    plt.savefig(ruta, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfica guardada: {ruta}")


def graficar_huella_quimica_heatmap(df_raw, df_detected, df_muestras_info, resumen, output_dir='figuras/inorganico'):
    """
    mapa de calor binario (Detectado vs < LD) de los analitos en agua cruda (ST),
    ordenado por frecuencia de detección; filas con analitos también vistos en
    BLK se marcan con un asterisco para referencia de control de calidad.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Na_OES excluido (ver nota en graficar_tasa_deteccion): además, al tener celdas
    # NA no es representable en el mapa binario sin introducir una tercera categoría.
    # Analitos con 0 detecciones en ST también se excluyen: no aportan nada al mapa
    # (serían una fila 100% gris) y solo abultan la figura.
    analitos_con_deteccion = set(resumen.loc[resumen['det_ST'] > 0, 'analito'])
    idx_plot = df_raw.index[(df_raw['analito'] != 'Na_OES') & (df_raw['analito'].isin(analitos_con_deteccion))]

    muestras_st = df_muestras_info[df_muestras_info['tipo'] == 'ST'].sort_values(by='lavanderia')
    cols_st = muestras_st['muestra_id'].tolist()

    df_st = df_detected.loc[idx_plot, cols_st].copy()
    df_st.index = df_raw.loc[idx_plot, 'analito']

    marca_blk = resumen.set_index('analito')['detectado_en_blk']
    etiquetas = [f"{a}*" if marca_blk.get(a, False) else a for a in df_st.index]
    df_st.index = etiquetas

    orden = df_st.sum(axis=1).sort_values(ascending=False).index
    df_st = df_st.loc[orden]

    plt.figure(figsize=(7, 13))
    cmap = sns.color_palette(['#eceff1', '#0b5394'])

    ax = sns.heatmap(
        df_st.astype(int),
        cmap=cmap,
        cbar_kws={'label': 'Condición Analítica', 'ticks': [0.25, 0.75], 'shrink': 0.3},
        linewidths=0.5,
        linecolor='#cfd8dc'
    )
    cbar = ax.collections[0].colorbar
    cbar.ax.set_yticklabels(['< LD (Censurado)', '> LD (Detectado)'], fontsize=9)
    cbar.ax.set_ylabel('Condición Analítica', fontsize=9)

    plt.xlabel('Muestra de Agua Cruda (ST)', fontsize=10, fontweight='bold', labelpad=8)
    plt.ylabel('Analito (* = también detectado en BLK)', fontsize=10, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(fontsize=7.5)

    plt.tight_layout(pad=1.2)
    ruta = os.path.join(output_dir, '02_huella_quimica_metales_agua_cruda.png')
    plt.savefig(ruta, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfica guardada: {ruta}")


def _preparar_consistencia(df_raw, df_detected, df_muestras_info):
    """Calcula en cuántas de las 6 fuentes (L1-L6) se detecta cada analito en ST. Uso interno."""
    idx_plot = df_raw.index[df_raw['analito'] != 'Na_OES']  # ver nota en graficar_tasa_deteccion

    muestras_st = df_muestras_info[df_muestras_info['tipo'] == 'ST']
    cols_st = muestras_st['muestra_id'].tolist()  # 1 muestra ST por fuente, sin réplicas

    num_fuentes_presente = df_detected.loc[idx_plot, cols_st].sum(axis=1)
    df_consistencia = pd.DataFrame({
        'analito': df_raw.loc[idx_plot, 'analito'],
        'num_fuentes': num_fuentes_presente
    })
    return df_consistencia[df_consistencia['num_fuentes'] > 0].sort_values(
        by='num_fuentes', ascending=False
    )


def graficar_consistencia_pie(df_raw, df_detected, df_muestras_info, output_dir='figuras/inorganico'):
    """Pastel: distribución de los analitos detectados según cuántas de las 6 fuentes los tienen."""
    os.makedirs(output_dir, exist_ok=True)
    df_consistencia = _preparar_consistencia(df_raw, df_detected, df_muestras_info)

    conteo_dist = df_consistencia['num_fuentes'].value_counts().sort_index(ascending=False)
    colores_pie = ['#08519c', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7']
    etiquetas = [f'En las {n} fuentes' if n == 6 else f'En {n} fuentes' for n in conteo_dist.index]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.pie(
        conteo_dist, labels=etiquetas, autopct='%1.1f%%', startangle=140,
        colors=colores_pie[:len(conteo_dist)],
        wedgeprops={'edgecolor': 'white', 'linewidth': 2}, textprops={'fontsize': 11}
    )
    plt.tight_layout(pad=1.5)
    ruta = os.path.join(output_dir, '03a_consistencia_pie.png')
    plt.savefig(ruta, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfica guardada: {ruta}")


def graficar_consistencia_barras(df_raw, df_detected, df_muestras_info, output_dir='figuras/inorganico'):
    """Barras horizontales: en cuántas de las 6 fuentes aparece cada analito (solo los parcialmente detectados)."""
    os.makedirs(output_dir, exist_ok=True)
    df_consistencia = _preparar_consistencia(df_raw, df_detected, df_muestras_info)
    # solo analitos con detección parcial (1-5 de 6); los 100% y 0% ya se ven en figuras 01/02
    parciales = df_consistencia[df_consistencia['num_fuentes'] < 6].sort_values(by='num_fuentes', ascending=True)

    if len(parciales) == 0:
        print("-> Sin analitos de detección parcial; se omite 03b_consistencia_barras.png")
        return

    y_pos = np.arange(len(parciales))
    fig, ax = plt.subplots(figsize=(7, max(3, len(parciales) * 0.35)))
    ax.barh(y_pos, parciales['num_fuentes'], color='#3182bd', height=0.65)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(parciales['analito'], fontsize=10)
    ax.set_xlabel('Número de Fuentes con Presencia (máx. 6)', fontsize=11, fontweight='bold')
    ax.set_xlim(0, 6.8)
    ax.set_xticks(range(1, 7))
    ax.grid(axis='x', linestyle='--', alpha=0.6)
    for i, val in enumerate(parciales['num_fuentes']):
        ax.text(val + 0.1, i, f'{val}/6', va='center', fontsize=9, color='#333333')

    plt.tight_layout(pad=1.2)
    ruta = os.path.join(output_dir, '03b_consistencia_barras.png')
    plt.savefig(ruta, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfica guardada: {ruta}")


def graficar_comparacion_blanco(df_raw, df_concentraciones, df_muestras_info, resumen, output_dir='figuras/inorganico'):
    """
    Para los analitos también detectados en BLK: compara el valor en el blanco contra
    el rango (min-mediana-max) en las muestras de agua cruda (ST), escala log.
    Permite juzgar si la detección en BLK es relevante en magnitud o despreciable.
    """
    os.makedirs(output_dir, exist_ok=True)
    cols_st = df_muestras_info[df_muestras_info['tipo'] == 'ST']['muestra_id'].tolist()

    flagged = resumen[resumen['detectado_en_blk']].copy()
    idx = flagged.index
    flagged['blk_val'] = df_concentraciones.loc[idx, 'BLK']
    flagged['st_min'] = df_concentraciones.loc[idx, cols_st].min(axis=1)
    flagged['st_mediana'] = df_concentraciones.loc[idx, cols_st].median(axis=1)
    flagged['st_max'] = df_concentraciones.loc[idx, cols_st].max(axis=1)
    flagged['ratio_mediana_blk'] = flagged['st_mediana'] / flagged['blk_val'].replace(0, np.nan)

    flagged = flagged.sort_values(by='ratio_mediana_blk', ascending=True)
    y_pos = np.arange(len(flagged))

    fig, ax = plt.subplots(figsize=(8, max(4, len(flagged) * 0.38)))
    ax.hlines(y_pos, flagged['st_min'], flagged['st_max'], color='#9ecae1', linewidth=5, alpha=0.8,
              label='Rango en ST (min-max)')
    ax.plot(flagged['st_mediana'], y_pos, 'o', color='#08519c', markersize=7, label='Mediana en ST')
    ax.plot(flagged['blk_val'], y_pos, 'D', color='#b30000', markersize=7, label='Valor en BLK')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(flagged['analito'], fontsize=9)
    ax.set_xscale('log')
    ax.set_xlabel('Concentración (escala log, µg/L salvo Na_OES en mg/L)', fontsize=10, fontweight='bold')
    ax.grid(axis='x', linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', fontsize=8.5, framealpha=0.9)

    plt.tight_layout(pad=1.2)
    ruta = os.path.join(output_dir, '03c_comparacion_blanco_campo.png')
    plt.savefig(ruta, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfica guardada: {ruta}")
    return flagged[['analito', 'blk_val', 'st_min', 'st_mediana', 'st_max', 'ratio_mediana_blk']]


def graficar_boxplots(df_concentraciones, df_raw, df_muestras_info, resumen, output_dir='figuras/inorganico', umbral_mayoritario=1000.0):
    """
    Boxplots de concentración en ST (escala log). La separación mayoritarios/traza es
    empírica (mediana en ST >= umbral_mayoritario µg/L), no una categoría forzada:
    hay un salto natural en los datos entre ~150 y ~7,255 µg/L.
    """
    os.makedirs(output_dir, exist_ok=True)
    cols_st = df_muestras_info[df_muestras_info['tipo'] == 'ST']['muestra_id'].tolist()

    # Na_OES se excluye (mismo elemento que Na, medido por método alterno) y también
    # los analitos con 0 detecciones en ST (Sc, Hg, Te, Ag, Be): su "concentración"
    # imputada es solo LD/sqrt(2) constante, no aportan nada a un boxplot.
    analitos_con_deteccion = set(resumen.loc[resumen['det_ST'] > 0, 'analito'])
    df_plot = df_raw[(df_raw['analito'] != 'Na_OES') & (df_raw['analito'].isin(analitos_con_deteccion))].copy()
    medianas = df_concentraciones.loc[df_plot.index, cols_st].median(axis=1)

    grupos = {
        'mayoritarios': df_plot.index[medianas >= umbral_mayoritario],
        'traza': df_plot.index[medianas < umbral_mayoritario],
    }

    for nombre_grupo, idx in grupos.items():
        if len(idx) == 0:
            continue
        sub = df_concentraciones.loc[idx, cols_st].copy()
        sub['analito'] = df_raw.loc[idx, 'analito']
        orden = sub.set_index('analito')[cols_st].median(axis=1).sort_values(ascending=False).index
        melt = pd.melt(sub, id_vars=['analito'], value_vars=cols_st, var_name='muestra', value_name='concentracion')

        fig_h = max(4, len(idx) * 0.4)
        plt.figure(figsize=(9, fig_h))
        ax = sns.boxplot(
            data=melt, y='analito', x='concentracion', order=orden,
            hue='analito', legend=False, palette='Blues_r' if nombre_grupo == 'mayoritarios' else 'GnBu_r',
            fliersize=0, linewidth=1.0, boxprops=dict(alpha=0.85)
        )
        sns.stripplot(
            data=melt, y='analito', x='concentracion', order=orden,
            color='#08306b', alpha=0.65, size=5, jitter=0.2
        )
        ax.set_xscale('log')
        ax.set_xlabel('Concentración en ST (µg/L) — Escala Logarítmica', fontsize=11, fontweight='bold')
        ax.set_ylabel('Analito', fontsize=11, fontweight='bold')
        ax.tick_params(labelsize=9.5)

        plt.tight_layout(pad=1.2)
        ruta = os.path.join(output_dir, f'04_boxplots_{nombre_grupo}.png')
        plt.savefig(ruta, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"-> Gráfica guardada: {ruta}")


if __name__ == '__main__':
    df_raw, df_muestras_info, df_detected, df_concentraciones, sample_ids = cargar_datos()

    muestras_st = df_muestras_info[df_muestras_info['tipo'] == 'ST']['muestra_id'].tolist()
    total_celdas_st = len(df_raw) * len(muestras_st)
    censuradas_st = (~df_detected[muestras_st]).values.sum()

    print("[OK] Datos cargados correctamente:")
    print(f"     - Total de analitos: {len(df_raw)}")
    print(f"     - Muestras: {len(sample_ids)} ({', '.join(sample_ids)})")
    print(f"     - Agua Cruda (ST): {len(muestras_st)} muestras ({', '.join(muestras_st)})")
    print(f"     - Censura en ST: {censuradas_st}/{total_celdas_st} ({censuradas_st/total_celdas_st*100:.1f}%)")
    print(df_muestras_info.to_string(index=False))

    resumen = calcular_resumen_deteccion(df_raw, df_detected, df_muestras_info, sample_ids)

    detectados_st = resumen[resumen['det_ST'] > 0]
    no_detectados_st = resumen[resumen['det_ST'] == 0]
    en_blk = resumen[resumen['detectado_en_blk']]

    print("\n" + "=" * 75)
    print(" HALLAZGOS PRINCIPALES DE DETECCIÓN EN AGUA CRUDA (ST):")
    print("=" * 75)
    print(f"• De {len(resumen)} analitos monitoreados:")
    print(f"  - {len(detectados_st)} ({len(detectados_st)/len(resumen)*100:.1f}%) se detectan en AL MENOS 1 muestra de agua cruda.")
    print(f"  - {len(no_detectados_st)} ({len(no_detectados_st)/len(resumen)*100:.1f}%) son 100% CENSURADOS (<LD) en agua cruda.")

    print(f"\n• Analitos presentes en el blanco de campo (BLK) — posible contaminación de muestreo: {len(en_blk)}")
    if len(en_blk) > 0:
        print(f"  {', '.join(en_blk['analito'].tolist())}")

    print("\n• Top analitos más frecuentes en Agua Cruda (ST):")
    cols_mostrar = ['id', 'analito', 'unidad', 'det_ST', 'pct_det_ST', 'detectado_en_blk']
    print(detectados_st.sort_values(by='pct_det_ST', ascending=False)[cols_mostrar].head(15).to_string(index=False))

    print("\n" + "=" * 75)
    print(" GENERANDO VISUALIZACIONES...")
    print("=" * 75)
    graficar_tasa_deteccion(resumen, output_dir='figuras/inorganico')
    graficar_huella_quimica_heatmap(df_raw, df_detected, df_muestras_info, resumen, output_dir='figuras/inorganico')
    graficar_consistencia_pie(df_raw, df_detected, df_muestras_info, output_dir='figuras/inorganico')
    graficar_consistencia_barras(df_raw, df_detected, df_muestras_info, output_dir='figuras/inorganico')
    tabla_blk = graficar_comparacion_blanco(df_raw, df_concentraciones, df_muestras_info, resumen, output_dir='figuras/inorganico')
    graficar_boxplots(df_concentraciones, df_raw, df_muestras_info, resumen, output_dir='figuras/inorganico')

    print("\n• Comparación BLK vs. ST (mediana/BLK, ordenado de más a menos preocupante):")
    print(tabla_blk.sort_values(by='ratio_mediana_blk').to_string(index=False))

    print("\n[ÉXITO] EDA de metales completado. Revisa las figuras en './figuras/inorganico'.")
