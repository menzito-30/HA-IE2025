# Research: Modelos de 2 Firmas (Formal/Informal)

Este documento resume los hallazgos de los papers de referencia para adaptar la lógica de modelos de 2 sectores (agente representativo) al modelo de agentes heterogéneos en tiempo continuo.

---

## 1. Restrepo-Echavarría (EER 2014): "Macroeconomic Volatility: The Role of The Informal Economy"

### Resumen

- **Tipo de modelo**: DSGE dinámico de equilibrio general
- **Contexto**: Conexión entre economías formal e informal a través de asignación de trabajo
- **Resultado clave**: El sector informal actúa como "amortiguador" de shocks económicos

### Estructura del Modelo

**Hogares deciden cómo asignar trabajo:**
$$
L = L_F + L_I
$$

Donde:

- $L_F$: Trabajo en sector formal
- $L_I$: Trabajo en sector informal

**Firmas:**

- **Sector Formal**: Producción con capital y trabajo, sujeto a regulaciones
- **Sector Informal**: Producción (posiblemente solo trabajo), menos regulaciones

### Ideas Clave para Adaptar a HA

1. **Asignación endógena de trabajo**: El hogar elige cuánto trabajar en cada sector basándose en salarios relativos
2. **Diferencias salariales**: $w_F > w_I$ típicamente
3. **El informal como "buffer"**: Absorbe trabajo cuando hay shocks negativos al sector formal

---

## 2. Horvath (2018): "Business Cycles and the Informal Economy"

### Resumen

- **Tipo de modelo**: Two-sector Real Business Cycle (RBC) de economía pequeña abierta
- **Enfoque**: Países emergentes
- **Resultado clave**: Relación positiva entre volatilidad de consumo/output y tamaño del sector informal

### Hallazgos Principales

1. **Aumento de tasas de interés**:
   - Contrae output, consumo, inversión y horas formales
   - **Expande el sector informal** (efecto sustitución)

2. **Volatilidad**:
   - Países con mayor sector informal experimentan mayor volatilidad en output, inversión y consumo

3. **Informalidad contracíclica**:
   - Cuando el sector formal se contrae, el informal se expande

### Estructura Típica

**Preferencias (hogar representativo):**
$$
U = \sum_{t=0}^{\infty} \beta^t \left[ \frac{c_t^{1-\sigma}}{1-\sigma} - \chi_F \frac{\ell_{F,t}^{1+\phi}}{1+\phi} - \chi_I \frac{\ell_{I,t}^{1+\phi}}{1+\phi} \right]
$$

**Producción:**

- Formal: $Y_F = A_F K^{\alpha} L_F^{1-\alpha}$
- Informal: $Y_I = A_I L_I$ (solo trabajo, típicamente)

---

## 3. Adaptación a Modelo HA en Tiempo Continuo

### Diferencias Clave: Agente Representativo vs Heterogéneo

| Aspecto | Agente Representativo | Agentes Heterogéneos |
|---------|----------------------|----------------------|
| Distribución | No hay | $g(a,z)$ estacionaria |
| Decisión de trabajo | Agregada | Por $(a,z)$ |
| Agregación | Trivial | $L_j = \int z \ell_j(a,z) g(a,z) \, da \, dz$ |
| Aseguración | Completa (implícita) | Incompleta (auto-aseguración) |

### Modificaciones Necesarias

#### 3.1 HJB con 2 Sectores

$$
\rho v(a,z) = \max_{c, \ell_F, \ell_I} \left\{ u(c, \ell_F, \ell_I) + v'(a,z) \cdot \dot{a} + \lambda_z [v(a,z') - v(a,z)] \right\}
$$

Con:
$$
\dot{a} = w_F z \ell_F + w_I \theta z \ell_I + ra - c
$$

#### 3.2 FOCs Intra-temporales

De la función de utilidad separada:
$$
u(c, \ell_F, \ell_I) = \frac{c^{1-\gamma}}{1-\gamma} - \psi_F \frac{\ell_F^{1+1/\varphi_F}}{1+1/\varphi_F} - \psi_I \frac{\ell_I^{1+1/\varphi_I}}{1+1/\varphi_I}
$$

FOC para $\ell_F$:
$$
\ell_F(a,z) = \left( \frac{v'(a,z) \cdot w_F z}{\psi_F} \right)^{\varphi_F}
$$

FOC para $\ell_I$:
$$
\ell_I(a,z) = \left( \frac{v'(a,z) \cdot w_I \theta z}{\psi_I} \right)^{\varphi_I}
$$

#### 3.3 Equilibrio con 2 Firmas

**Firma Formal** (usa capital):
$$
Y_F = A_F K^{\alpha} L_F^{1-\alpha}
$$

FOCs:

- $r + \delta = \alpha A_F (L_F/K)^{1-\alpha}$
- $w_F = (1-\alpha) A_F (K/L_F)^{\alpha}$

**Firma Informal** (solo trabajo):
$$
Y_I = A_I L_I
$$

FOC:

- $w_I = A_I$ (constante)

**Clearing de mercados:**

- Capital: $K = \int a \cdot g(a,z) \, da \, dz$
- Trabajo Formal: $L_F = \int z \cdot \ell_F(a,z) \cdot g(a,z) \, da \, dz$
- Trabajo Informal: $L_I = \int \theta z \cdot \ell_I(a,z) \cdot g(a,z) \, da \, dz$

---

## 4. Parámetros Sugeridos (Basado en Literatura)

| Parámetro | Símbolo | Valor Típico | Fuente |
|-----------|---------|--------------|--------|
| Atenuación shock informal | $\theta$ | 0.5 - 0.8 | Calibrar |
| Productividad informal relativa | $A_I/A_F$ | 0.3 - 0.6 | Restrepo, Horvath |
| Peso desutilidad informal | $\psi_I/\psi_F$ | < 1 (trabajo informal "más fácil") | A calibrar |
| Elasticidad Frisch | $\varphi$ | 0.5 - 2.0 | Literatura macro |

---

## 5. Finuras Matemáticas al Adaptar

> [!WARNING]
> **Cuidados al pasar de Agente Representativo a HA:**

1. **No linealidad de agregación**:
   $$L_F \neq \ell_F \cdot \bar{z}$$
   En HA, la distribución $g(a,z)$ afecta cómo se agregan las decisiones

2. **Precautionary savings**:
   En HA, la incertidumbre idiosincrática genera ahorro precautorio ausente en modelos representativos

3. **Distribución endógena**:
   La forma de $g(a,z)$ depende de las políticas $\ell_F, \ell_I, c$, que a su vez dependen de los precios

4. **Equilibrio de punto fijo**:
   - Dado $(r, w_F, w_I)$: resolver HJB → políticas → KFE → distribución
   - Agregación: $K, L_F, L_I$ → verificar FOCs de firmas
   - Ajustar precios hasta clearing

---

## 6. Próximos Pasos de Implementación

1. [ ] Modificar `codigo_prueba_v5.m` para incluir 2 variables de trabajo
2. [ ] Añadir cálculo de $\ell_F, \ell_I$ en cada punto de la grilla
3. [ ] Modificar drift: $s = w_F z \ell_F + w_I \theta z \ell_I + ra - c$
4. [ ] Añadir firma informal con $w_I = A_I$ (constante)
5. [ ] Actualizar equilibrio: buscar $(r, w_F)$ tal que clearing se cumpla
