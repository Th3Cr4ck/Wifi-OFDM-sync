Referencia:[[/home/isaac/Documents/Cinvestav/Tesis/Thesis References/S&C.pdf|S&C]]
## 1. Problema: Carrier Frequency Offset (CFO)

En un sistema OFDM, el receptor puede tener una diferencia de frecuencia respecto al transmisor. Esta diferencia se conoce como **Carrier Frequency Offset (CFO)**.

Una señal recibida puede modelarse como:
$$
r[n] = s[n]e^{j2\pi \Delta f nT_s}+w[n]
$$
donde:

- $s[n]$: señal transmitida.    
- $\Delta f$: diferencia de frecuencia entre transmisor y receptor.
- $T_s$: periodo de muestreo.
- $w[n]$: ruido.

El término
$$
e^{j2\pi \Delta f nT_s}
$$
representa una **rotación progresiva de fase** causada por el CFO.

En OFDM, esta rotación es problemática porque puede destruir la ortogonalidad entre las subportadoras, provocando:

- Inter-Carrier Interference (ICI).
- Rotación de la constelación.
- Degradación de la detección de símbolos.

Por lo tanto, el receptor necesita estimar $\Delta f$ para posteriormente compensarlo.

---
## 2. Estructura especial del preámbulo

El símbolo se construye de manera que sus dos mitades sean idénticas:
$$
s[n+L]=s[n]
$$
donde $L$ es la longitud de cada mitad.

Por lo tanto, si no existiera CFO, dos muestras separadas exactamente $L$ posiciones serían iguales.

---
# 3. Efecto del CFO sobre las dos mitades

La señal recibida puede modelarse como:
$$
r[n]=s[n]e^{j2\pi\Delta f nT_s}
$$
Consideremos una muestra de la segunda mitad:
$$
r[n+L]
=
s[n+L]e^{j2\pi\Delta f(n+L)T_s}
$$
Como las dos mitades son iguales:
$$
s[n+L]=s[n]
$$
entonces:
$$
r[n+L]
=
s[n]e^{j2\pi\Delta f(n+L)T_s} 
$$
$$
r[n+L]
=
s[n]
e^{j2\pi\Delta f nT_s}
e^{j2\pi\Delta f LT_s}
$$
Pero la señal recibida en la primera mitad es:
$$
r[n]
=
s[n]e^{j2\pi\Delta f nT_s}
$$
Por lo tanto:
$$
\boxed{
r[n+L]
=
r[n]e^{j2\pi\Delta f LT_s}
}
$$
Esto significa que las muestras correspondientes de las dos mitades tienen la misma información de la señal, pero la segunda mitad presenta una **rotación de fase constante** respecto a la primera:
$$
\boxed{
\phi=2\pi\Delta f LT_s
}
$$
donde $\phi$ es la diferencia de fase introducida por el CFO entre las dos mitades.

Esta relación es la base de la estimación de frecuencia de Schmidl & Cox: si podemos medir $\phi$, podemos despejar el CFO:
$$
\Delta f
=
\frac{\phi}{2\pi LT_s}
$$
y posteriormente utilizar esta estimación para compensar la rotación de fase producida por el CFO.

---

## 4. Correlación entre las dos mitades

De acuerdo con el punto anterior:
$$
r[n+L]
=
r[n]e^{j2\pi\Delta f LT_s}
$$
Por lo tanto, podemos comparar las dos mitades mediante una correlación.

Schmidl & Cox definen:
$$
P(d)
=
\sum_{n=0}^{L-1}
r^*[d+n]r[d+n+L]
$$
donde:
* $d$ es una posible posición dentro de la señal recibida.
* $L$ es la longitud de cada mitad del símbolo.
* $r^*[n]$ es el conjugado complejo de la señal.

Para entender qué contiene esta correlación, consideremos un término:
$$
r^*[n]r[n+L]
$$
Usando la relación obtenida anteriormente:
$$
r[n+L]
=
r[n]e^{j2\pi\Delta f LT_s}
$$
tenemos:
$$
r^*[n]r[n+L]
=
r^*[n]r[n]
e^{j2\pi\Delta f LT_s}
$$
$$
r^*[n]r[n+L]
=
|r[n]|^2e^{j2\pi\Delta f LT_s}
$$
Entonces:
$$
P(d)
=
\sum_{n=0}^{L-1}
|r[d+n]|^2
e^{j2\pi\Delta f LT_s}
$$
Como la exponencial no depende de $n$, podemos sacarla de la suma:
$$
P(d)
=
\left(
\sum_{n=0}^{L-1}|r[d+n]|^2
\right)
e^{j2\pi\Delta f LT_s}
$$
El primer término es una cantidad real y positiva, por lo que no aporta fase.

Por lo tanto:
$$
\boxed{
\angle P(d)
=
2\pi\Delta f LT_s
}
$$
Es decir, **la fase de la correlación contiene directamente la información del CFO**.

---

## 5. Estimación del CFO

Partimos de:
$$
\angle P(d)
=
2\pi\Delta f LT_s
$$
Despejando $\Delta f$:
$$
\boxed{
\hat{\Delta f}
=
\frac{\angle P(d)}
{2\pi LT_s}
}
$$
Esta expresión proporciona el CFO estimado en Hz.

Sin embargo, en OFDM es común **normalizar el CFO respecto al espaciamiento entre subportadoras**.

Si el tamaño de la FFT es $N$, el espaciamiento entre subportadoras es:
$$
\Delta f_{\text{sub}}
=
\frac{1}{NT_s}
$$
Definimos el CFO normalizado como:
$$
\epsilon
=
\frac{\Delta f}{\Delta f_{\text{sub}}}
$$
Entonces:
$$
\Delta f
=
\frac{\epsilon}{NT_s}
$$
Sustituyendo esta expresión en:
$$
\angle P(d)
=
2\pi\Delta f LT_s
$$
obtenemos:
$$
\angle P(d)
=
2\pi
\left(
\frac{\epsilon}{NT_s}
\right)
LT_s
$$
$$
\angle P(d)
=
2\pi\epsilon\frac{L}{N}
$$
Despejando $\epsilon$:
$$
\boxed{
\hat{\epsilon}
=
\frac{N}{2\pi L}
\angle P(d)
}
$$
Para el caso en que las dos mitades tienen longitud:
$$
L=\frac{N}{2}
$$
la expresión se simplifica a:
$$
\boxed{
\hat{\epsilon}
=
\frac{1}{\pi}
\angle P(d)
}
$$
Por lo tanto, para este caso, el CFO normalizado se obtiene directamente a partir de la fase de la correlación.

---

## 6. Limitación: ambigüedad del CFO

El estimador anterior tiene una limitación importante.

La función de fase de un número complejo, $\angle(\cdot)$, devuelve normalmente un valor dentro del intervalo:
$$
-\pi < \angle P(d) \leq \pi
$$
Como:
$$
\angle P(d)
=
2\pi\epsilon\frac{L}{N}
$$
debemos tener:
$$
-\pi
<
2\pi\epsilon\frac{L}{N}
\leq
\pi
$$
Dividiendo entre $2\pi L/N$:
$$
-\frac{N}{2L}
<
\epsilon
\leq
\frac{N}{2L}
$$
Por lo tanto, el rango de CFO normalizado que puede estimarse sin ambigüedad es:
$$
\boxed{
-\frac{N}{2L}
<
\epsilon
\leq
\frac{N}{2L}
}
$$
Si:
$$
L=\frac{N}{2}
$$
entonces:
$$
\boxed{
-1<\epsilon\leq1
}
$$
Por lo tanto, el estimador basado únicamente en la periodicidad de dos mitades permite estimar aproximadamente un CFO dentro de:
$$
|\epsilon|\leq1
$$
es decir, dentro de aproximadamente **un espaciamiento entre subportadoras**.

### ¿Por qué aparece la ambigüedad?

La razón es que una fase de $2\pi$ es indistinguible de una fase de $0$:
$$
e^{j(\phi+2\pi k)}
=
e^{j\phi}
$$
para cualquier entero $k$.

Por ejemplo:
$$
e^{j0.3(2\pi)}
=
e^{j(1.3)(2\pi)}
$$
Por lo tanto, el estimador basado solamente en la correlación no puede distinguir entre:
$$
\epsilon=0.3
$$
y
$$
\epsilon=1.3
$$
porque ambos producen la misma fase después de considerar la periodicidad de $2\pi$.

Por esta razón, Schmidl & Cox necesitan **información adicional del entrenamiento** para resolver esta ambigüedad y obtener una estimación del CFO con un rango mayor.

En términos conceptuales:
$$
\boxed{
\text{Correlación entre mitades}
\rightarrow
\text{CFO fraccionario}
}
$$
mientras que la información adicional del entrenamiento permite determinar la **componente entera** del CFO.

---
## 7. Estimación y corrección de la componente entera del CFO

El CFO puede expresarse como:
$$
\boxed{
\Delta f
=
\Delta f_{\mathrm{frac}}
+
\frac{2g}{T}
}
$$
donde:
* $\Delta f_{\mathrm{frac}}$ es la componente fraccionaria.
* $g\in\mathbb{Z}$ es un entero.
* $T$ es la duración útil del símbolo OFDM, sin el intervalo de guarda.

En términos del CFO normalizado:
$$
\boxed{
\epsilon
=
\epsilon_{\mathrm{frac}}+2g
}
$$
para la estructura de entrenamiento utilizada por Schmidl & Cox.

La componente fraccionaria se obtiene mediante el primer símbolo:
$$
\hat{\epsilon}_{\mathrm{frac}}
=
\frac{\angle P(d)}{\pi}
$$
y primero se utiliza esta estimación para realizar una **corrección parcial** de la señal recibida:
$$
r_c[n]
=
r[n]
e^{-j2\pi\hat{\Delta f}_{\mathrm{frac}}nT_s}
$$
Después de esta corrección, todavía puede permanecer un CFO correspondiente a un número entero de desplazamientos de subportadora.

### Segundo símbolo de entrenamiento

Para resolver esta ambigüedad, Schmidl & Cox utilizan un segundo símbolo de entrenamiento construido de manera diferente al primero.

Mientras que el primer símbolo contiene una estructura periódica para estimar la parte fraccionaria, el segundo contiene una secuencia conocida con una relación diferencial entre sus subportadoras.

Después de corregir parcialmente los dos símbolos utilizando:
$$
\hat{\Delta f}_{\mathrm{frac}}
$$
se realiza la FFT de ambos:
$$
X_1[k]=FFT\{r_{c,1}[n]\}
$$
$$
X_2[k]=FFT\{r_{c,2}[n]\}
$$
Ahora el problema consiste en determinar cuánto se desplazaron las subportadoras debido a la componente entera del CFO.

### Desplazamiento producido por un CFO entero

Un CFO entero de una separación de subportadora produce un desplazamiento de los índices de las subportadoras.

Si el CFO normalizado contiene una componente:
$$
2g
$$
entonces una subportadora transmitida en el índice $k$ aparece en el receptor desplazada hacia:
$$
k+2g
$$
Por eso, para diferentes valores posibles de $g$, podemos intentar "desplazar" los índices recibidos y comprobar cuál produce la mayor correspondencia con la secuencia conocida transmitida.
### Función de decisión

Schmidl & Cox definen una métrica $B(g)$ que mide la correspondencia entre los dos símbolos de entrenamiento para un determinado desplazamiento entero $g$.

De forma general:
$$
B(g)
=
\frac{
\left|
\displaystyle\sum_{k\in X}
X_1^*[k+2g]\,v_k^*\,X_2[k+2g]
\right|^2
}{
2
\left(
\displaystyle\sum_{k\in X}|X_2[k]|^2
\right)^2
}
$$
donde:

* $X_1[k]$: FFT del primer símbolo de entrenamiento después de la corrección fraccionaria.
* $X_2[k]$: FFT del segundo símbolo.
* $v_k$: relación conocida entre las subportadoras de los dos símbolos de entrenamiento.
* $X$: conjunto de subportadoras utilizadas para la estimación.
* $g$: posible desplazamiento entero.

La idea es probar diferentes valores de $g$.

Por ejemplo:
$$
g=-2,-1,0,1,2,\ldots
$$
y calcular:
$$
B(-2),B(-1),B(0),B(1),B(2),\ldots
$$
El valor que produzca la mayor métrica se toma como la estimación del desplazamiento entero:
$$
\boxed{
\hat g
=
\arg\max_g B(g)
}
$$
---

## CFO total

Una vez obtenidas las dos componentes:
$$
\hat{\epsilon}_{\mathrm{frac}},\quad
\hat g
$$
el CFO total estimado es:

$$
\boxed{
\hat{\epsilon}
=
\hat{\epsilon}_{\mathrm{frac}}
+
2\hat g
}
$$
Por lo tanto, la estimación completa consiste en:
$$
\boxed{
\text{CFO}
=
\text{parte fraccionaria}
+
\text{parte entera}
}
$$
o:
$$
\boxed{
\hat{\epsilon}
=
\frac{\angle P(d)}{\pi}
+
2\hat g
}
$$