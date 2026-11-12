Core de la [[Arquitectura general del sistema HDL|arquitectura del receptor]].

Su funcionalidad principal es tomar las muestras provenientes del paquete generado por el [[Paquetizer-1|núcleo anterior]] y aplicar un proceso de **normalización bajo criterios de energı́a y magnitud**.

Se procede al cálculo del valor absoluto al cuadrado de cada muestra, representado mediante la siguiente expresión:

$$abs^2 = Re^2 + Im^2$$

Con el objetivo de **evitar el uso de bloques de procesamiento de alta complejidad**, como divisores, que incrementan significativamente el uso de recursos y latencias y no representan condiciones crı́ticas para esta etapa, se adopta un **enfoque alternativo basado en el uso de desplazamientos lógicos**. Para ello:

1. Se detecta el valor máximo del módulo entre todas las muestras almacenadas en el paquete. Este valor, que es estrictamente positivo, es acumulado internamente con la finalidad de calcular la energı́a total del paquete recibido.
2. Si la energı́a no supera un umbral mı́nimo, definido por el parámetro ENERGY UMBRAL, el paquete no es normalizado. En su lugar, se reenvı́a directamente al bloque siguiente, comportándose como un buffer tipo bypass.
3. En cambio, si el paquete sı́ supera al umbral energético, se procede a normalizar las muestras utilizando como referencia el valor máximo encontrado. Para este proceso se emplean únicamente desplazamientos lógicos en lugar de divisiones.

