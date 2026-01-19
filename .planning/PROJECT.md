# HA-IE2025: Aiyagari Continuo con Trabajo Endógeno y Dualidad Formal/Informal

## Visión del Proyecto

Desarrollar un modelo de **Agentes Heterogéneos (HA) en tiempo continuo** tipo Aiyagari, con:
- **Oferta de trabajo endógena** (labor supply)
- **Dos firmas/sectores**: Formal e Informal
- **Elección intensiva**: El agente distribuye su tiempo de trabajo entre ambos sectores

### Fundamentos Teóricos

El modelo se basa en:
1. **Achdou, Han, Lasry, Lions & Moll (2017)** - Framework matemático HACT
2. **Moll's Labor Supply Extension** - Huggett con trabajo endógeno [labor_supply.pdf](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/MOLL/modelo_base/labor_market_moll/labor_supply.pdf)
3. **Aiyagari (1994)** - Modelo base con acumulación de capital

### Estrategia Incremental

```mermaid
graph LR
    A[Huggett + Endog Labor<br/>Moll Base ✓] --> B[Aiyagari 1 Firma<br/>+ Endog Labor 🔄]
    B --> C[Aiyagari 2 Firmas<br/>Formal/Informal ⏳]
```

| Fase | Descripción | Estado |
|------|-------------|--------|
| Base Moll | Huggett con labor endógeno | ✅ Código verificado |
| Aiyagari 1 Firma | Adaptar lógica Moll a Aiyagari | 🔄 En desarrollo |
| Aiyagari 2 Firmas | Dualidad formal/informal | ⏳ Pendiente |

---

## Estructura del Repositorio

```
HA-IE2025/
├── .planning/                    # Documentación de proyecto
│   ├── PROJECT.md               # Este archivo
│   ├── REQUIREMENTS.md          # Requerimientos v1/v2
│   ├── ROADMAP.md               # Fases y progreso
│   ├── STATE.md                 # Decisiones y bloqueos
│   └── research/                # Conocimiento matemático
├── Code/CONTINUOUS_TIME/
│   ├── Aiyagari_firmas/         # Implementación principal
│   │   ├── base_Moll/           # Código base Aiyagari (sin labor endog)
│   │   └── try_endog_labor_1_firm/  # TRABAJO EN CURSO
│   └── Hugget_bonos/
│       └── base_Moll/
│           └── endog_labor_1_firm/  # Referencia de Moll
└── pdfs/MOLL/modelo_base/       # Papers de referencia
```

---

## Tecnología

- **Lenguaje principal**: MATLAB (.m)
- **Alternativa**: Jupyter Notebooks (.ipynb) para experimentación
- **Método numérico**: Diferencias finitas, esquema upwind, método implícito

---

## Referencias Clave

| Documento | Contenido | Ruta |
|-----------|-----------|------|
| `Moll_teoria_HA.pdf` | Teoría HACT completa | [pdfs/MOLL/modelo_base/](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/MOLL/modelo_base/Moll_teoria_HA.pdf) |
| `HACT_Numerical_Appendix.pdf` | Detalles numéricos | [pdfs/MOLL/modelo_base/](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/MOLL/modelo_base/HACT_Numerical_Appendix.pdf) |
| `labor_supply.pdf` | Trabajo endógeno (Huggett) | [pdfs/MOLL/modelo_base/labor_market_moll/](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/MOLL/modelo_base/labor_market_moll/labor_supply.pdf) |
| `achdou_income_wealth_distribution.pdf` | Paper base HACT | [pdfs/MOLL/modelo_base/](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/MOLL/modelo_base/achdou_income_wealth_distribution.pdf) |
