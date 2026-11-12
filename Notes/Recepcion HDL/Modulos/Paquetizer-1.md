Core de la [[Arquitectura general del sistema HDL|arquitectura del receptor]].

Este bloque funcional ha sido diseñado para recibir las muestras digitales provenientes de la etapa anterior y organizarlas en paquetes de tamaño definido, conforme a una señal de conteo interna regulada por el parámetro configurable *NPACKAGE*.

El Paquetizer almacena temporalmente las muestras en una **memoria interna** que actúa como búfer, permitiendo la acumulación de datos hasta alcanzar el tamaño especificado del paquete.