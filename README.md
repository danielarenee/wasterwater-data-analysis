# PFAS y Sustancias Restringidas en Efluentes Textiles de Lavanderías Industriales

Análisis de datos del proyecto de servicio social. Se estudian muestras de agua residual de
5 lavanderías industriales (L1–L5): agua **sin tratar (ST)**, agua **tratada (T)** y una
muestra de lodo. Hay dos análisis, uno por laboratorio:

| Análisis | Laboratorio | Qué se mide |
|---|---|---|
| **Orgánico** | CIATEC | 97 compuestos: 29 PFAS + 68 sustancias restringidas |
| **Inorgánico** | Actlabs (informe A26-15546) | 60 metales/elementos por ICP-MS (Na también por ICP-OES) |

Los dos siguen el mismo flujo: **cargar datos → análisis exploratorio (EDA) → imputación de
valores censurados → figuras y tablas**.

---

## Estructura del repositorio

```
serviciosocial/
├── eda_org.py              # EDA orgánico
├── imputacion_org.py       # Imputación orgánica
├── preparar_metales.py     # Excel de Actlabs → CSV
├── eda_inorg.py            # EDA inorgánico
├── imputacion_inorg.py     # Imputación inorgánica
│
├── outputs/                # CSV generados por los scripts
├── figuras/
│   ├── organico/
│   └── inorganico/
└── docs/                   # Normatividad, artículo de referencia y notas para el reporte
```

### Inputs

Solo hay **dos archivos de entrada**. Todo lo que está en `outputs/` y `figuras/` se genera
a partir de ellos.

| Archivo | Contenido |
|---|---|
| `muestras.csv` | Resultados de CIATEC. Compuestos en filas y muestras en columnas (`L1_ST_1`, `L1_T`, …). La segunda fila tiene el código de laboratorio (`XA01-1`, …). |
| `Información/A26-15546Final.xlsx` | Informe de Actlabs, hoja `Results`. Lo lee `preparar_metales.py`. |

El resto de `Información/` (el certificado en PDF, `Parámetros muestreo_Compilación_Andrea.xlsx`
y la compilación de datos de muestreo) es material de consulta. **Ningún script lo lee.**

### Outputs (`outputs/`)

| Archivo | Lo genera | Contenido |
|---|---|---|
| `resumen_deteccion_completo.csv` | `eda_org.py` | Tasa de detección de los 97 compuestos |
| `muestras_24.csv`, `muestras_24_st.csv` | `eda_org.py` | Los 24 compuestos prioritarios (todas las muestras / solo ST) |
| `resumen_24_compuestos.csv` | `eda_org.py` | Resumen de los 24 compuestos prioritarios |
| `muestras_24_st_imputado_epa.csv` | `imputacion_org.py` | Los 24 compuestos en ST, con los censurados imputados |
| `metales_actlabs.csv` | `preparar_metales.py` | El Excel de Actlabs en el mismo formato que `muestras.csv` |
| `metales_actlabs_imputado_epa.csv` | `imputacion_inorg.py` | Metales con los censurados imputados |

---

## Cómo correrlo

### 1. Instalar dependencias

Se probó con Python 3.13.

```bash
pip install pandas numpy matplotlib seaborn openpyxl
```

### 2. Ejecutar los scripts en orden

Los scripts se corren **desde la raíz del repo**, porque las rutas son relativas.

```bash
# Orgánico
python eda_org.py
python imputacion_org.py      # usa outputs/muestras_24_st.csv, que genera eda_org.py

# Inorgánico
python preparar_metales.py    # crea outputs/metales_actlabs.csv
python eda_inorg.py
python imputacion_inorg.py
```

Cada script imprime en la terminal un resumen de los hallazgos: conteos de detección,
porcentaje de censura y avisos. Además, guarda sus figuras y CSV.

Los scripts de imputación importan funciones del EDA correspondiente
(`imputacion_org.py` usa `eda_org.py` e `imputacion_inorg.py` usa `eda_inorg.py`), así que
tienen que estar en la misma carpeta.

---

## Conceptos y decisiones clave

**Nomenclatura de muestras:** `L{n}_{tipo}_{réplica}`.
- `ST` = sin tratar (agua cruda), `T` = tratada, `L` = lodo.
- Ejemplo: `L3_ST_2` es la segunda réplica de agua cruda de la lavandería 3.

**Valores censurados:** cuando el laboratorio no puede cuantificar un compuesto, reporta
`<LPC` (CIATEC, límite práctico de cuantificación) o `<LD` (Actlabs, límite de detección).
Por ejemplo, `<0.1`. En los metales también aparece `>100000` cuando el equipo se satura.
Los CSV guardan estos valores como texto para no perder esa información.

**Imputación US EPA (LPC/√2 o LD/√2):** los censurados se sustituyen por el límite dividido
entre √2, según Hornung & Reed (1990). Los scripts de imputación generan los histogramas y
el cálculo de GSD que justifican el método.

**Los 24 compuestos prioritarios (orgánico):** de los 97 compuestos solo se conservan los que
se detectan en **3 o más lavanderías**. En agua cruda, el 76% de los datos orgánicos están
censurados, así que con los demás no se puede hacer estadística.

**Metales:** la censura es mucho menor (alrededor del 19% en ST), así que se imputa la base
completa. Hay dos exclusiones:
- `Na_OES`, porque es el mismo sodio medido con otro método;
- los analitos que no se detectaron en ninguna muestra de proceso.

**L6 y BLK (solo en metales):**
- **L6** no es una de las 5 lavanderías. Es la textilera *El Asturcón*, muestreada el 9 de
  julio en otra campaña. Antes de compararla con L1–L5, hay que tomar esto en cuenta.
- **BLK** es un blanco de campo que se entregó junto con las muestras. Si un elemento aparece
  en BLK, parte de lo que se mide en las muestras podría venir del muestreo o del manejo.
  `eda_inorg.py` compara la concentración en las muestras con la del blanco.

---

## Pendientes

- Comparar las concentraciones con los límites de la NOM-001-SEMARNAT-2021 y de la
  Declaratoria del río Atoyac.
- Cruzar los resultados orgánicos con los inorgánicos
