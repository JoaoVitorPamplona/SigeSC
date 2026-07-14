# Data Sources

All input data are publicly available from the Brazilian Institute of Geography and Statistics (IBGE). This file documents each required dataset, its source URL, and where to place it in the repository.

## Directory Structure

After downloading, your `data/raw/` directory should look like:

```
data/raw/
├── cnefe/
│   └── 42_SC.csv                          # CNEFE address registry
├── aggregates/
│   ├── Agregados_por_setores_basico_SC.csv
│   └── Agregados_por_setores_caracteristicas_domicilio2_SC.csv
├── sidra/
│   ├── feminino_idade_raça.xlsx            # SIDRA Table 9606 (Female)
│   ├── masculino_idade_raça.xlsx           # SIDRA Table 9606 (Male)
│   ├── quantidadedemoradores_idade_F.xlsx  # SIDRA Table 9880 (Female)
│   ├── quantidadedemoradores_idade_M.xlsx  # SIDRA Table 9880 (Male)
│   ├── municipios_SC_IBGE2022.csv          # Municipality codes
│   └── TabelaNacionalidade.xlsx            # SIDRA Table 10157
├── microdata_2010/
│   └── pop_filtrado.csv                    # 2010 Census microdata (filtered)
└── pop_estimates/
    └── Extensao2025.xlsx                   # IBGE 2025 population estimates
```

## Dataset Details

### 1. CNEFE — National Address Registry (Residence Registry)

- **Description**: Georeferenced dwelling units for Santa Catarina from the 2022 Census.
- **Source**: https://www.ibge.gov.br/estatisticas/sociais/populacao/38734-cadastro-nacional-de-enderecos-para-fins-estatisticos.html
- **File**: `42_SC.csv` (state code 42 = Santa Catarina)
- **Place in**: `data/raw/cnefe/`
- **Used by**: `01_residence_registry.jl`

### 2. Census Aggregates by Tract

- **Description**: Aggregate counts of households by census tract, including dwelling subtypes, bathroom counts, and head demographics.
- **Source**: https://www.ibge.gov.br/estatisticas/sociais/saude/22827-censo-demografico-2022.html (Aggregates section)
- **Files**:
  - `Agregados_por_setores_basico_SC.csv` — Basic tract-level counts
  - `Agregados_por_setores_caracteristicas_domicilio2_SC.csv` — Dwelling characteristics
- **Place in**: `data/raw/aggregates/`
- **Used by**: `01_residence_registry.jl`, `02_residence_attributes.jl`, `10_allocation_model.jl`

### 3. SIDRA Tables — Population and Family Marginals

#### Table 9606: Population by age, sex, and ethnicity
- **Source**: https://sidra.ibge.gov.br/tabela/9606
- **Files**: `feminino_idade_raça.xlsx`, `masculino_idade_raça.xlsx`
- **Place in**: `data/raw/sidra/`
- **Used by**: `03_population_families.jl`

#### Table 9880: Households by number of residents, head demographics
- **Source**: https://sidra.ibge.gov.br/tabela/9880
- **Files**: `quantidadedemoradores_idade_F.xlsx`, `quantidadedemoradores_idade_M.xlsx`
- **Place in**: `data/raw/sidra/`
- **Used by**: `03_population_families.jl`

#### Table 10157: Population by nationality
- **Source**: https://sidra.ibge.gov.br/tabela/10157
- **File**: `TabelaNacionalidade.xlsx`
- **Place in**: `data/raw/sidra/`
- **Used by**: `06_nationality.jl`

### 4. 2010 Census Microdata Sample

- **Description**: Individual-level microdata from the 2010 Census sample, used as the structural seed for IPF.
- **Source**: https://www.ibge.gov.br/estatisticas/sociais/populacao/9662-censo-demografico-2010.html
- **File**: `pop_filtrado.csv` (pre-filtered for Santa Catarina; see note below)
- **Place in**: `data/raw/microdata_2010/`
- **Used by**: `04_head_assignment.jl`, `05_ipf_member_allocation.jl`
- **Note**: The raw IBGE microdata file is large (~10 GB nationally). The file `pop_filtrado.csv` is a pre-filtered extract containing only Santa Catarina records with the columns used in this pipeline. A filtering script will be provided in a future release.

### 5. IBGE 2025 Population Estimates

- **Description**: Municipal population estimates for 2025, used to update the synthetic population from 2022 to 2025.
- **Source**: https://ftp.ibge.gov.br/Estimativas_de_Populacao/Estimativas_2025/
- **File**: `Extensao2025.xlsx`
- **Place in**: `data/raw/pop_estimates/`
- **Used by**: `07_population_update.jl`

### 6. Municipality Code List

- **Description**: IBGE codes for the 295 municipalities of Santa Catarina.
- **File**: `municipios_SC_IBGE2022.csv`
- **Place in**: `data/raw/sidra/`
- **Used by**: `03_population_families.jl`, `06_nationality.jl`

## Notes

- All SIDRA tables should be downloaded for the state of Santa Catarina (code 42), at the municipality level.
- The 2022 Census microdata sample had not been released by IBGE as of July 2026. When available, it should replace the 2010 microdata to improve intra-household correlation fidelity.
