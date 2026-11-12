Módulo central de comunicación que controla la interacción coordinada entre los cores del [[Arquitectura general del sistema HDL|receptor]].

Garantiza un acceso ordenado y eficiente a los recursos compartidos, así como el control del flujo de información a lo largo de la arquitectura.

## Formas de implementación

- Mediante un **firmware** ejecutándose sobre un núcleo RISC-V bajo el sistema operativo en  tiempo real Zephyr, lo que permite una solución embebida completamente integrada.
- Mediante una computadora personal que opere bajo un **sistema operativo convencional**, como Windows o Linux, ofreciendo una plataforma de control y supervisión más flexible durante las etapas de desarrollo o prueba.
