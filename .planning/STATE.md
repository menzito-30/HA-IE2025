# Estado del Proyecto: Aiyagari HACT 2 Firmas

> **Última actualización**: 2026-01-18

---

## Posición Actual

**Fase**: 2 - Aiyagari 1 Firma + Labor Endógeno  
**Milestone**: 2.2 (Transición de asset supply a equilibrio general)  
**Archivo activo**: [codigo_prueba_v5.m](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/Code/CONTINUOUS_TIME/Aiyagari_firmas/try_endog_labor_1_firm/asset_supply_poisson/codigo_prueba_v5.m)

---

## Decisiones Tomadas

| Fecha | Decisión | Justificación |
|-------|----------|---------------|
| 2025 | Usar proceso de Poisson (2 estados) | Simplifica la implementación inicial |
| 2025 | Grilla equiespaciada de activos | Consistente con código base de Moll |
| 2025 | Método implícito para HJB | Estabilidad numérica |

---

## Decisiones Pendientes → **RESUELTAS**

| ID | Pregunta | Decisión Final |
|----|----------|----------------|
| **PD1** | ¿Cómo diferenciar trabajo formal/informal? | **Salarios diferentes** $w_F > w_I$ + **costos fijos opcionales** |
| **PD2** | ¿Desutilidad del trabajo conjunta o separada? | **Separada**: $\psi_F \ell_F^{...} + \psi_I \ell_I^{...}$ (evita corner solutions) |
| **PD3** | ¿El shock $z$ afecta ambos sectores? | **z afecta menos al sector informal** (ej: $z_I = \theta z$ con $\theta < 1$) |
| **PD4** | ¿Mercado de capital único o segmentado? | Pendiente (simplificar con mercado único primero) |

### Referencias Adicionales (Agente Representativo → adaptar a HA)

- [restrepo_EER_R&R.pdf](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/otros_autores/restrepo_EER_R&R.pdf): 2 firmas
- [horvath bussines cycles.pdf](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/otros_autores/horvath%20bussines%20cycles.pdf): 2 firmas
- [financial frictions Giovanni Merlin-Vladimir Teles.pdf](file:///c:/Users/johnb/Documents/GitHub/HA-IE2025/pdfs/otros_autores/financial%20frictions%20Giovanni%20Merlin-Vladimir%20Teles.pdf)

> [!WARNING]
> Estos papers son de agente representativo. La lógica debe adaptarse con cuidado a HA (distribución, agregación, equity premium).

---

## Bloqueos Actuales

| ID | Descripción | Tipo | Acción Requerida |
|----|-------------|------|------------------|
| **B1** | Verificar estabilidad numérica de `codigo_prueba_v5.m` | Técnico | Ejecutar y analizar resultados |
| **B2** | Definir forma funcional para 2 firmas | Teórico | Decisión del usuario (PD1, PD2) |

---

## Próximos Pasos Inmediatos

1. [ ] Ejecutar `codigo_prueba_v5.m` y verificar gráficos
2. [ ] Comparar resultados con Aiyagari base (sin labor endógeno)
3. [ ] **Resolver decisiones PD1-PD3** antes de proceder a Fase 3
4. [ ] Implementar bisección para equilibrio general en 1 firma

---

## Contexto para Sesiones Futuras

### Archivos Clave

- **Base Moll (referencia)**: `Code/CONTINUOUS_TIME/Hugget_bonos/base_Moll/endog_labor_1_firm/`
- **Trabajo en progreso**: `Code/CONTINUOUS_TIME/Aiyagari_firmas/try_endog_labor_1_firm/asset_supply_poisson/`
- **Aiyagari base (sin labor endog)**: `Code/CONTINUOUS_TIME/Aiyagari_firmas/base_Moll/`

### Diferencias Clave Huggett vs Aiyagari

| Aspecto | Huggett | Aiyagari |
|---------|---------|----------|
| Activo | Bonos | Capital |
| $r$ | Exógeno o equilibrio en bonos | $r = F_K - \delta$ |
| Firma | No hay | Cobb-Douglas con $K$ y $L$ |
| Restricción | $a \geq \underline{a}$ (puede ser negativo) | $a \geq 0$ (no short capital) |

### Estructura del Código `codigo_prueba_v5.m`

```
Líneas 1-60:    Parámetros y grillas
Líneas 60-120:  Inicialización, boundary labor
Líneas 120-260: Loop interno HJB (upwind + implícito)
Líneas 260-300: KFE y distribución
Líneas 300-350: Gráficos y guardado
```
