# Fundamentos Matemáticos: Aiyagari HACT con Labor Endógeno

Este documento resume las "finuras" matemáticas de los papers de referencia aplicadas al caso de **Aiyagari en tiempo continuo con oferta de trabajo endógena**.

---

## 1. The Household Problem

### 1.1 Ecuación de Bellman (HJB)

El problema del hogar es maximizar la utilidad esperada descontada:

$$
\rho v(a,z) = \max_{c,\ell} \left\{ u(c,\ell) + \frac{\partial v}{\partial a}(a,z) \cdot \dot{a} + \lambda_z \left[ v(a,z') - v(a,z) \right] \right\}
$$

Donde:

- $v(a,z)$: función de valor
- $a$: activos (capital en Aiyagari)
- $z \in \{z_1, z_2\}$: productividad idiosincrática (proceso de Poisson)
- $\ell$: oferta de trabajo
- $\lambda_z$: intensidad de transición entre estados $z$

### 1.2 Función de Utilidad con Labor

$$
u(c,\ell) = \frac{c^{1-\gamma}}{1-\gamma} - \frac{\ell^{1+\frac{1}{\varphi}}}{1+\frac{1}{\varphi}}
$$

Donde:

- $\gamma$: coeficiente de aversión al riesgo
- $\varphi$: elasticidad de Frisch de la oferta laboral

### 1.3 Restricción Presupuestaria

$$
\dot{a} = wz\ell + ra - c
$$

Con restricción de endeudamiento: $a \geq \underline{a}$

### 1.4 Condiciones de Primer Orden (FOC)

**FOC para consumo:**
$$
u_c = c^{-\gamma} = v'(a,z)
$$

**FOC intratemporal (consumo-trabajo):**
$$
\ell^{\frac{1}{\varphi}} \cdot c^{\gamma} = wz
$$

> [!IMPORTANT]
> Combinando ambas FOC:
> $$\ell = (v'(a,z) \cdot wz)^{\varphi}$$
> Esta es la ecuación clave que conecta la derivada de la función de valor con la oferta laboral.

---

## 2. Extensión a 2 Firmas (Formal/Informal) — ESPECIFICACIÓN FINAL

### 2.1 Nueva Restricción Presupuestaria

$$
\dot{a} = w_F z \ell_F + w_I \theta z \ell_I + ra - c
$$

Donde:

- $\ell_F$: trabajo en sector formal (salario $w_F$)
- $\ell_I$: trabajo en sector informal (salario $w_I < w_F$)
- $\theta \in (0,1)$: factor de atenuación del shock $z$ en sector informal

> [!IMPORTANT]
> **Decisión clave**: El shock de productividad $z$ afecta **menos** al sector informal. Esto captura que en la informalidad la productividad idiosincrática puede ser menos relevante.

### 2.2 Nueva Función de Utilidad (Desutilidad Separada)

$$
u(c,\ell_F,\ell_I) = \frac{c^{1-\gamma}}{1-\gamma} - \psi_F \frac{\ell_F^{1+\frac{1}{\varphi_F}}}{1+\frac{1}{\varphi_F}} - \psi_I \frac{\ell_I^{1+\frac{1}{\varphi_I}}}{1+\frac{1}{\varphi_I}}
$$

Donde:

- $\psi_F, \psi_I > 0$: pesos de desutilidad (normalizar $\psi_F = 1$)
- $\varphi_F, \varphi_I$: elasticidades de Frisch (pueden ser iguales inicialmente)

> [!NOTE]
> La desutilidad **separada** evita corner solutions y permite soluciones interiores $(\ell_F > 0, \ell_I > 0)$.

### 2.3 FOC Intra-temporal (2 Firmas)

**FOC trabajo formal:**
$$
\psi_F \ell_F^{\frac{1}{\varphi_F}} = c^{-\gamma} \cdot w_F z = v'(a,z) \cdot w_F z
$$

**FOC trabajo informal:**
$$
\psi_I \ell_I^{\frac{1}{\varphi_I}} = c^{-\gamma} \cdot w_I \theta z = v'(a,z) \cdot w_I \theta z
$$

**Despejando las ofertas laborales:**
$$
\ell_F = \left( \frac{v'(a,z) \cdot w_F z}{\psi_F} \right)^{\varphi_F}
$$
$$
\ell_I = \left( \frac{v'(a,z) \cdot w_I \theta z}{\psi_I} \right)^{\varphi_I}
$$

### 2.4 Costo Fijo de Entrada (Opcional)

Si se desea añadir un costo fijo $\kappa$ para trabajar en el sector formal:

$$
\dot{a} = w_F z \ell_F - \kappa \cdot \mathbf{1}_{\ell_F > 0} + w_I \theta z \ell_I + ra - c
$$

> [!CAUTION]
> El costo fijo introduce **no-convexidades**. Se recomienda primero implementar el modelo SIN costo fijo y añadirlo como extensión.

### 2.5 Firmas

**Firma Formal (usa capital y trabajo):**
$$
Y_F = A_F K^{\alpha} L_F^{1-\alpha}
$$

**Firma Informal (solo trabajo):**
$$
Y_I = A_I L_I
$$

> [!IMPORTANT]
> El sector informal **solo usa trabajo**, no capital. Esto simplifica el equilibrio ya que solo la firma formal demanda capital.

**FOCs Firma Formal:**
$$
r + \delta = \alpha A_F \left(\frac{L_F}{K}\right)^{1-\alpha}
$$
$$
w_F = (1-\alpha) A_F \left(\frac{K}{L_F}\right)^{\alpha}
$$

**FOC Firma Informal:**
$$
w_I = A_I
$$

> [!NOTE]
> $w_I = A_I$ es constante (no depende de $L_I$). Esto es típico de modelos con sector informal de "autoempleo" o productividad marginal constante.

---

## 3. Stationary Equilibrium

### 3.1 Definición Formal

Un equilibrio estacionario consiste en:

1. **Función de valor** $v(a,z)$ y **políticas** $c(a,z), \ell(a,z)$ que resuelven el problema del hogar dado $(r,w)$

2. **Distribución estacionaria** $g(a,z)$ que satisface la **Ecuación de Kolmogorov Forward (KFE)**:
   $$
   0 = -\frac{\partial}{\partial a}[s(a,z) \cdot g(a,z)] + \lambda_{z'} g(a,z') - \lambda_z g(a,z)
   $$

3. **Precios de firma** que satisfacen FOCs:
   $$
   r = F_K(K,L) - \delta = \alpha A K^{\alpha-1} L^{1-\alpha} - \delta
   $$
   $$
   w = F_L(K,L) = (1-\alpha) A K^{\alpha} L^{-\alpha}
   $$

4. **Clearing de mercados**:
   - Capital: $K = \int a \cdot g(a,z) \, da \, dz$
   - Trabajo: $L = \int z \cdot \ell(a,z) \cdot g(a,z) \, da \, dz$

### 3.2 Equilibrio con 2 Firmas

Para el caso de 2 firmas, se añade:

- **Firma Formal**: $Y_F = A_F K_F^{\alpha_F} L_F^{1-\alpha_F}$
- **Firma Informal**: $Y_I = A_I K_I^{\alpha_I} L_I^{1-\alpha_I}$

**Clearing de trabajo**:
$$
L_F = \int z \cdot \ell_F(a,z) \cdot g(a,z) \, da \, dz
$$
$$
L_I = \int z \cdot \ell_I(a,z) \cdot g(a,z) \, da \, dz
$$

---

## 4. Computational Algorithm

### 4.1 Paso a Paso (1 Firma con Labor Endógeno)

```
OUTER LOOP: Sobre r (bisección)
│
├── 1. Dado r, calcular precios de firma:
│      K_D = (αA/(r+δ))^(1/(1-α)) * L
│      w = (1-α)A(K_D/L)^α
│
├── 2. Pre-calcular l0(a,z): labor con ahorro cero
│      Resolver: l - (wzl + ra)^(-γφ)(wz)^φ = 0  ∀(a,z)
│
└── INNER LOOP: Iteración de función de valor
    │
    ├── 3. Calcular derivadas (diferencias finitas)
    │      v'_F = (v_{i+1} - v_i)/Δa    (forward)
    │      v'_B = (v_i - v_{i-1})/Δa    (backward)
    │
    ├── 4. Políticas candidatas:
    │      c_F = (v'_F)^(-1/γ),  l_F = (v'_F · wz)^φ,  s_F = wzl_F + ra - c_F
    │      c_B = (v'_B)^(-1/γ),  l_B = (v'_B · wz)^φ,  s_B = wzl_B + ra - c_B
    │      c_0 = wzl_0 + ra,     l_0 = pre-calculado
    │
    ├── 5. Esquema UPWIND:
    │      Si s_F > 0: usar (c_F, l_F)     [ahorrando]
    │      Si s_B < 0: usar (c_B, l_B)     [desahorrando]
    │      Si cambio de signo: usar (c_0, l_0)
    │
    ├── 6. Construir matriz de transición A
    │      A = A_drift + A_switch
    │
    ├── 7. Resolver sistema implícito:
    │      B·v^{n+1} = u^n + v^n/Δt
    │      con B = (1/Δt + ρ)I - A
    │
    └── 8. Verificar convergencia: ||v^{n+1} - v^n|| < ε

POST-CONVERGENCE:
├── 9. Resolver KFE: A'g = 0 (normalizar)
├── 10. Calcular oferta de capital: K_S = ∫a·g(a,z)dadz
└── 11. Actualizar r por bisección hasta K_S = K_D
```

### 4.2 Modificaciones para 2 Firmas

1. **Pre-cálculo**: Resolver sistema 2x2 para $(l_{0,F}, l_{0,I})$ en cada $(a,z)$
2. **FOCs**: Calcular $\ell_F, \ell_I$ simultáneamente dada $v'$
3. **Equilibrio**: Loop adicional sobre $(w_F, w_I)$ o $(K_F, K_I)$

---

## 5. Grid Construction

### 5.1 Grilla de Activos

```matlab
a = linspace(amin, amax, I)';  % Equiespaciada
da = (amax - amin)/(I - 1);
```

> [!TIP]
> Para más precisión cerca de $\underline{a}$, usar grilla cuadrática:
>
> ```matlab
> a = amin + (amax-amin)*linspace(0,1,I)'.^2;
> ```

### 5.2 Parámetros Recomendados

| Parámetro | Valor Típico | Nota |
|-----------|--------------|------|
| $I$ (grid points) | 500-1000 | Mayor para 2 firmas |
| $\underline{a}$ | -0.15 (Huggett), 0 (Aiyagari) | Natural borrowing limit |
| $\bar{a}$ | 20 | Ajustar según $A_{prod}$ |
| $\Delta$ (step size) | 1000 | Para estabilidad implícita |
| $\epsilon$ (tolerancia) | $10^{-6}$ | Convergencia HJB |

---

## 6. Labor Supply Solution

### 6.1 Sistema a Resolver (Ahorro Cero)

En cada punto $(a_i, z_j)$, encontrar $\ell$ tal que:

$$
\ell = (wz\ell + ra)^{-\gamma\varphi} (wz)^{\varphi}
$$

Equivalentemente, $f(\ell) = 0$ donde:

$$
f(\ell) = \ell - (wz\ell + ra)^{-\gamma\varphi} (wz)^{\varphi}
$$

### 6.2 Implementación en MATLAB

```matlab
function eq = lab_solve(l, params)
    a = params(1); z = params(2); w = params(3);
    r = params(4); gamma = params(5); frisch = params(6);
    
    eq = l - (w*z*l + r*a)^(-gamma*frisch) * (w*z)^frisch;
end

% Uso:
for i = 1:I
    for j = 1:2
        params = [a(i), z(j), w, r, gamma, frisch];
        l0(i,j) = fzero(@(l) lab_solve(l, params), 1);
    end
end
```

---

## 7. Transition Matrices

### 7.1 Matriz de Transición por Drift (A_drift)

La discretización del operador $\mathcal{A}v = s(a,z) \frac{\partial v}{\partial a}$ produce:

$$
A_{drift} = \begin{pmatrix} A_1 & 0 \\ 0 & A_2 \end{pmatrix}
$$

Donde $A_j$ es tridiagonal con:

- **Diagonal inferior** ($X$): $-\min(s_B, 0)/\Delta a$
- **Diagonal principal** ($Y$): $-\max(s_F, 0)/\Delta a + \min(s_B, 0)/\Delta a$
- **Diagonal superior** ($Z$): $\max(s_F, 0)/\Delta a$

### 7.2 Matriz de Transición por Saltos Poisson (A_switch)

$$
A_{switch} = \begin{pmatrix} -\lambda_1 I & \lambda_1 I \\ \lambda_2 I & -\lambda_2 I \end{pmatrix}
$$

> [!IMPORTANT]
> **Propiedad fundamental**: Las filas de $A = A_{drift} + A_{switch}$ deben sumar cero.
>
> ```matlab
> if max(abs(sum(A,2))) > 1e-9
>     error('Improper Transition Matrix');
> end
> ```

### 7.3 Extensión para 2 Firmas

La estructura de $A$ no cambia sustancialmente, pero el **drift** $s(a,z)$ ahora depende de $\ell_F$ y $\ell_I$:

$$
s(a,z) = w_F z \ell_F(a,z) + w_I z \ell_I(a,z) + ra - c(a,z)
$$

---

## 8. Condiciones de Convexidad y Regularidad

### 8.1 Condiciones Suficientes para Solución Única

1. **Utilidad estrictamente cóncava**: $u_{cc} < 0$, $u_{\ell\ell} < 0$
2. **Borrowing constraint natural**: $wz_{\min}\ell + r\underline{a} \geq 0$
3. **Impaciencia suficiente**: $\rho > r$ (en equilibrio)

### 8.2 Regularidad de la Solución

- $v(a,z)$ es creciente y cóncava en $a$
- $v'(a,z) > 0$ para todo $(a,z)$
- $c(a,z), \ell(a,z)$ son funciones suaves de $a$

> [!CAUTION]
> Si $v' \leq 0$ en algún punto, indica problemas numéricos. Aplicar:
>
> ```matlab
> dV = max(dV, 1e-10);  % Floor para estabilidad
> ```

---

## 9. Preguntas Abiertas para Diseño

1. **Desutilidad del trabajo**: ¿Una función conjunta $\ell_F + \ell_I$ o separada?
2. **Diferenciación sectorial**: ¿Por salarios $w_F \neq w_I$, elasticidades $\varphi_F \neq \varphi_I$, o costos fijos?
3. **Productividad por sector**: ¿El shock $z$ afecta igual a ambos sectores?
4. **Capital por sector**: ¿$K = K_F + K_I$ con mercado único de capital, o mercados segmentados?
