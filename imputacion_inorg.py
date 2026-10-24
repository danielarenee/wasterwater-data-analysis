"""
Script de Imputación EPA (LD / sqrt(2)) — Metales/Elementos (Actlabs)
Mismo método que imputacion_org.py, aplicado a la base completa de metales
(no hay un "universo prioritario" reducido para metales, ver docs/hallazgos_generales.md).
"""

import os
import csv
import numpy as np
import pandas as pd

import eda_inorg as eda

ANALITO_EXCLUIDO = 'Na_OES'  # mismo elemento que Na, medido por método alterno; no se imputa


def calcular_gsd_resumen(df_raw, df_detected, sample_ids):
    """
    GSD de los valores detectados para cada analito con censura parcial, para
    justificar (o señalar excepciones a) la sustitución LD/sqrt(2) — Hornung & Reed (1990).
    """
    filas = []
    for i, row in df_raw.iterrows():
        if row['analito'] == ANALITO_EXCLUIDO:
            continue
        det_col = df_detected.loc[i, sample_ids]
        raw_strs = row[sample_ids].astype(str).str.replace('<', '').str.replace('>', '').str.replace(' ', '')
        numeric_vals = pd.to_numeric(raw_strs, errors='coerce')
        det_vals = numeric_vals[det_col.values.astype(bool)].dropna()
        n_censurados = int((~det_col).sum())
        n_detectados = int(det_col.sum())

        if n_censurados == 0:
            continue  # sin censura, no requiere imputación

        gsd = np.exp(np.std(np.log(det_vals), ddof=1)) if n_detectados >= 2 else np.nan
        filas.append((row['analito'], n_censurados, n_detectados, gsd))

    return pd.DataFrame(filas, columns=['analito', 'n_censurados', 'n_detectados', 'GSD']).sort_values(
        by='GSD', ascending=False, na_position='first'
    )


def imputar_dataset_epa(df_raw, df_detected, sample_ids):
    """
    Imputa concentraciones censuradas (<LD) usando LD/sqrt(2). Las celdas "no medidas"
    (Na_OES donde no aplicó ICP-OES) se conservan como NaN, no se imputan.
    """
    df_imputado = pd.DataFrame(index=df_raw.index)

    for i, row in df_raw.iterrows():
        ld = row['LD']
        det_col = df_detected.loc[i, sample_ids]
        raw_strs = row[sample_ids].astype(str).str.replace('<', '').str.replace('>', '').str.replace(' ', '')
        numeric_vals = pd.to_numeric(raw_strs, errors='coerce')

        for col in sample_ids:
            estado = det_col[col]
            if pd.isna(estado):
                df_imputado.loc[i, col] = np.nan  # no medido (Na_OES)
            elif estado:
                df_imputado.loc[i, col] = numeric_vals[col]  # detectado, valor tal cual
            else:
                df_imputado.loc[i, col] = ld / np.sqrt(2.0)  # censurado <LD

    return df_imputado


def exportar_dataset_imputado(df_raw, df_full_original, sample_ids, df_imputado, output_csv='metales_actlabs_imputado_epa.csv', output_dir='outputs'):
    """
    Exporta la base completa de metales con los valores censurados ya imputados,
    conservando la estructura de outputs/metales_actlabs.csv (fila de códigos Actlabs,
    columnas id/analito/unidad/LD/metodo + muestras).
    """
    os.makedirs(output_dir, exist_ok=True)

    meta_cols = ['id', 'analito', 'unidad', 'LD', 'metodo']
    fila_codigos = df_full_original.iloc[0]  # fila de trazabilidad (codigo_actlabs)

    fila_codigos_out = ['' if pd.isna(fila_codigos[c]) else fila_codigos[c] for c in meta_cols] + \
                        [fila_codigos[c] for c in sample_ids]
    filas = [fila_codigos_out]
    for i, row in df_raw.iterrows():
        vals = [f"{df_imputado.loc[i, c]:.6g}" if pd.notna(df_imputado.loc[i, c]) else '' for c in sample_ids]
        meta_vals = [int(row['id']) if pd.notna(row['id']) else ''] + [row[c] for c in meta_cols[1:]]
        filas.append(meta_vals + vals)

    ruta = os.path.join(output_dir, output_csv)
    with open(ruta, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(meta_cols + sample_ids)
        writer.writerows(filas)
    print(f"-> Base de metales imputada (EPA LD/sqrt2) exportada: {ruta}")


def main():
    df_raw, df_muestras_info, df_detected, df_concentraciones, sample_ids = eda.cargar_datos()
    df_full_original = pd.read_csv('outputs/metales_actlabs.csv', encoding='utf-8-sig')

    cols_proceso = df_muestras_info[df_muestras_info['tipo'].isin(['ST', 'T'])]['muestra_id'].tolist()
    idx = df_raw.index[df_raw['analito'] != ANALITO_EXCLUIDO]

    total_celdas = df_detected.loc[idx, cols_proceso].size
    censuradas = (~df_detected.loc[idx, cols_proceso]).values.sum()
    print("=" * 75)
    print(" DIAGNÓSTICO DE CENSURA (muestras de proceso ST+T, excluyendo Na_OES)")
    print("=" * 75)
    print(f"Celdas censuradas (<LD): {censuradas}/{total_celdas} ({censuradas/total_celdas*100:.1f}%)")
    print("Muy por debajo del 76.0% observado en el universo orgánico; hay margen para")
    print("usar LD/sqrt(2) en vez de una sustitución más simple (LD/2).")

    print("\n" + "=" * 75)
    print(" GSD DE ANALITOS CON CENSURA PARCIAL (justificación Hornung & Reed, 1990)")
    print("=" * 75)
    gsd_resumen = calcular_gsd_resumen(df_raw, df_detected, cols_proceso)
    print(gsd_resumen.to_string(index=False))

    excepciones = gsd_resumen[gsd_resumen['GSD'] >= 3]
    print(f"\n[AVISO] {len(excepciones)} analito(s) con GSD >= 3 (alta dispersión; la sustitución")
    print("LD/sqrt(2) es menos robusta para estos, aunque se aplica igual por consistencia):")
    print(excepciones.to_string(index=False) if len(excepciones) else "  (ninguno)")

    muy_censurados = gsd_resumen[gsd_resumen['n_detectados'] < 2]
    if len(muy_censurados):
        print(f"\n[AVISO] Analitos con <2 valores detectados (GSD no estimable, imputación poco informativa):")
        print(muy_censurados.to_string(index=False))

    # Analitos sin ninguna detección en muestras de proceso (ST+T): su "concentración"
    # sería solo LD/sqrt(2) constante en todas las muestras, sin ninguna información
    # cuantitativa real. Se excluyen de la base exportada (igual que ya se excluyen de
    # las figuras de eda_inorg.py), no solo de Na_OES.
    det_proceso = df_detected.loc[idx, cols_proceso].sum(axis=1)
    sin_deteccion = df_raw.loc[idx].loc[det_proceso == 0, 'analito'].tolist()
    print(f"\n[INFO] Analitos sin ninguna detección en proceso (excluidos de la base exportada): {sin_deteccion}")

    idx_exportar = idx[~df_raw.loc[idx, 'analito'].isin(sin_deteccion)]
    df_raw_exportar = df_raw.loc[idx_exportar]

    print("\n" + "=" * 75)
    print(" IMPUTANDO BASE FINAL (LD/sqrt(2))...")
    print("=" * 75)
    df_imputado = imputar_dataset_epa(df_raw_exportar, df_detected, sample_ids)
    exportar_dataset_imputado(df_raw_exportar, df_full_original, sample_ids, df_imputado)

    print("\n[ÉXITO] Imputación EPA completada.")


if __name__ == '__main__':
    main()
