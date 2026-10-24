# %% [1] Configuración
# Lee el reporte de Actlabs (A26-15546, metales/elementos por ICP-MS) y construye
# un CSV con la misma lógica que el de CIATEC: analitos en filas, muestras en columnas,
# valores como texto ("<0.1", ">100000") para conservar la censura.
from pathlib import Path
import pandas as pd

RUTA_XLSX = Path("Información/A26-15546Final.xlsx")
RUTA_CSV  = Path("outputs/metales_actlabs.csv")

# Código Actlabs -> ID con la convención del proyecto L{n}_{ST|T}_{réplica}.
# E (entrada) = sin tratar, S (salida) = tratada. El sufijo es la fecha de muestreo
# (mmdd) y coincide con el muestreo CIATEC del 28-30/10/2025.
# OJO: XL06 es la textilera El Asturcón (muestreada 07/09), no una de las 5 lavanderías.
MAPA_MUESTRAS = {
    "XL01-E-1028": "L1_ST_1", "XL01-S-1028": "L1_T_1",
    "XL02-E-1029": "L2_ST_1", "XL02-S-1029": "L2_T_1",
    "XL03-E-1029": "L3_ST_1", "XL03-S-1029": "L3_T_1",
    "XL04-E-1029": "L4_ST_1", "XL04-S-1029": "L4_T_1",
    "XL05-E-1030": "L5_ST_1", "XL05-S-1030": "L5_T_1",
    "XL06-E-0709": "L6_ST_1", "XL06-S-0709": "L6_T_1",
    "BLK": "BLK",   # blanco entregado con las muestras (no es el blanco de método del lab)
}

# %% [2] Leer la hoja de resultados sin encabezados
# La hoja trae 2 filas de título y luego filas de metadatos identificadas por su
# etiqueta en la columna 0; las buscamos por nombre en vez de por posición.
raw = pd.read_excel(RUTA_XLSX, sheet_name="Results", header=None, dtype=str)
etiqueta = raw[0].str.strip()

def fila(nombre):
    return raw.loc[etiqueta == nombre].iloc[0, 2:]   # cols 0-1 = código y descripción

simbolos = fila("Analyte Symbol").str.strip()
unidades = fila("Unit Symbol").str.strip()
lds      = fila("Detection Limit").str.strip()
metodos  = fila("Analysis Method").str.strip()

# Na aparece dos veces: ICP-MS (ug/L, satura en >100000) e ICP-OES (mg/L, solo
# reportado donde ICP-MS saturó). Renombramos para no tener símbolos duplicados.
simbolos = simbolos.where(metodos != "ICP-OES", simbolos + "_OES")

# %% [3] Matriz de valores (muestras en filas por ahora)
muestras = raw.loc[etiqueta.isin(MAPA_MUESTRAS.keys())]
valores = muestras.iloc[:, 2:].copy()
valores.columns = simbolos.values
valores.index = muestras[0].str.strip()

# Normaliza el texto: "< 0.1" -> "<0.1", "> 100000" -> ">100000"; celdas vacías quedan NaN
valores = valores.apply(lambda c: c.str.replace(r"\s+", "", regex=True))

# %% [4] Armar el CSV final: analitos en filas, muestras en columnas
meta = pd.DataFrame({
    "id": range(1, len(simbolos) + 1),
    "analito": simbolos.values,
    "unidad": unidades.values,
    "LD": pd.to_numeric(lds.values),
    "metodo": metodos.values,
})
datos = valores.T.reset_index(drop=True)
datos.columns = [MAPA_MUESTRAS[c] for c in valores.index]
tabla = pd.concat([meta, datos], axis=1)

# Fila extra con el código original de Actlabs, para trazabilidad (como en CIATEC)
codigos = pd.DataFrame([[""] * meta.shape[1] + list(valores.index)], columns=tabla.columns)
codigos.iloc[0, 1] = "codigo_actlabs"
tabla = pd.concat([codigos, tabla], ignore_index=True)

RUTA_CSV.parent.mkdir(exist_ok=True)
tabla.to_csv(RUTA_CSV, index=False, encoding="utf-8-sig")
print(f"{RUTA_CSV}: {len(simbolos)} analitos x {datos.shape[1]} muestras")