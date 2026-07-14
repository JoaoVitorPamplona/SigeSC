# Data Sources

All input data are publicly available from the Brazilian Institute of Geography and Statistics (IBGE).

## Generated datasets (GitHub Release)

The pipeline-generated datasets required to run the code are available as zip files in the [GitHub Release (`v1.0`)](../../releases/tag/v1.0).

Download the four zip files and extract them into the `output/` directory:

```
output/
├── stage2_attributes/
│   └── correspondencia_setores.csv
├── stage3_population/
│   ├── dataset_familias_2025.csv
│   └── dataset_pessoas_com_Ordem_final.csv
└── stage4_update_2025/
    └── residencias_2025.csv
```

| File | Zip size | Description |
|------|----------|-------------|
| `dataset_pessoas_com_Ordem_final.csv` | 48 MB | Synthetic persons with household order (k = 1–15) |
| `dataset_familias_2025_v2.csv` | 10 MB | Synthetic families projected to 2025 |
| `residencias_2025.csv` | 45 MB | Georeferenced residence registry (2022 + new constructions) |
| `correspondencia_setores.csv` | 176 KB | Census tract code matching (CNEFE ↔ aggregates) |

## IBGE public data

One additional file is sourced directly from IBGE:

### Census aggregates by tract

- **File**: `Agregados_por_setores_caracteristicas_domicilio2_SC.csv`
- **Description**: Aggregate counts of occupied private dwellings per census tract, stratified by dwelling subtype, number of bathrooms, and head gender/ethnicity.
- **Source**: https://www.ibge.gov.br/estatisticas/sociais/saude/22827-censo-demografico-2022.html (Aggregates section → Santa Catarina)
- **Place in**: `data/raw/aggregates/`
- **Used by**: `10_allocation_model.jl`, `benchmark_lp_vs_milp.jl`

## Municipality code list

- **File**: `municipios_SC_IBGE2022.csv` (included in this repository under `data/`)
- **Description**: IBGE codes for the 295 municipalities of Santa Catarina (2022 Census).
