# Roadmap: Aiyagari HACT con Labor Endógeno y 2 Firmas

## Fase 1: Validación del Código Base ✅

**Objetivo**: Verificar que el código de Moll (Huggett + labor endógeno) funciona correctamente.

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Ejecutar `HJB_labor_supply.m` | [Hugget_bonos/.../endog_labor_1_firm/](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/Code/CONTINUOUS_TIME/Hugget_bonos/base_Moll/endog_labor_1_firm/HJB_labor_supply.m) | ✅ |
| Entender `lab_solve.m` | [lab_solve.m](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/Code/CONTINUOUS_TIME/Hugget_bonos/base_Moll/endog_labor_1_firm/lab_solve.m) | ✅ |
| Verificar gráficos de políticas | N/A | ✅ |

---

## Fase 2: Aiyagari 1 Firma + Labor Endógeno 🔄

**Objetivo**: Adaptar la lógica de Moll (Huggett) al modelo de Aiyagari con firma representativa.

### Milestone 2.1: Curva de Oferta de Activos ✅

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Implementar firma (FOC capital/trabajo) | [codigo_prueba_v5.m](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/Code/CONTINUOUS_TIME/Aiyagari_firmas/try_endog_labor_1_firm/asset_supply_poisson/codigo_prueba_v5.m) | ✅ |
| Loop sobre `r_grid` | `codigo_prueba_v5.m` | ✅ |
| Gráfico $S(r)$ vs $K^D(r)$ | `codigo_prueba_v5.m` | ✅ |

### Milestone 2.2: Equilibrio General ⏳

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Implementar bisección en $r$ | TBD | ⏳ |
| Verificar clearing de mercados | TBD | ⏳ |
| Añadir equilibrio en mercado de trabajo $L$ | TBD | ⏳ |

### Milestone 2.3: Verificación y Debugging 🔄

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Revisar condiciones de borde | `codigo_prueba_v5.m` L132-142 | 🔄 |
| Verificar estabilidad numérica de `fzero` | `codigo_prueba_v5.m` | 🔄 |
| Comparar resultados con Aiyagari base (sin labor) | [aiyagari_poisson_steadystate.m](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/Code/CONTINUOUS_TIME/Aiyagari_firmas/base_Moll/aiyagari_poisson_steadystate.m) | ⏳ |

---

## Fase 3: Aiyagari 2 Firmas (Formal/Informal) ⏳

**Objetivo**: Extender el modelo a 2 sectores con decisión intensiva de trabajo.

### Milestone 3.1: Diseño Teórico

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Definir utilidad con 2 tipos de trabajo | research/MATHEMATICAL_FOUNDATIONS.md | ⏳ |
| Derivar FOCs para $(\ell_F, \ell_I)$ | research/MATHEMATICAL_FOUNDATIONS.md | ⏳ |
| Decidir diferenciación de sectores | .planning/STATE.md | ⏳ |

### Milestone 3.2: Implementación

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Modificar `lab_solve.m` para 2 sectores | TBD | ⏳ |
| Modificar HJB para 2 variables de control | TBD | ⏳ |
| Añadir 2 ecuaciones de firma | TBD | ⏳ |
| Resolver equilibrio con $(w_F, w_I)$ | TBD | ⏳ |

### Milestone 3.3: Análisis

| Tarea | Archivo | Estado |
|-------|---------|--------|
| Distribución de trabajo por sector y riqueza | TBD | ⏳ |
| Comparative statics: cambios en $A_F/A_I$ | TBD | ⏳ |
| Welfare analysis | TBD | ⏳ |

---

## Timeline Sugerido

```
Semana 1-2: Completar Fase 2 (Milestone 2.2, 2.3)
Semana 3:   Diseño Teórico Fase 3 (Milestone 3.1)
Semana 4-5: Implementación Fase 3 (Milestone 3.2)
Semana 6:   Análisis y documentación (Milestone 3.3)
```

---

## Riesgos y Mitigaciones

| Riesgo | Impacto | Mitigación |
|--------|---------|------------|
| Corner solutions en FOC 2 firmas | Alto | Introducir fricciones/costos diferenciados |
| Inestabilidad numérica con 4+ estados | Medio | Usar grilla más fina, preconditioners |
| Convergencia lenta en equilibrio general | Medio | Warm starts, bisección adaptativa |
