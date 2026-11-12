- [.] CFO Grueso
- [.] CFO Fino
- [.] Implementacion ambos

## Arquitecturas correccion grueso + fino
- [.] Correcion Grueso -> Correccion Fina
- [.] Estimacion Grueso + Fino -> Correcion Total
- [ ] Correccion Fina -> FFT -> Correccion Gruesa por desplazamiento (AGWN ✔, TDL ✘ falta Estimacion de canal)

## Metodo de estimacion STO
Actualmente, las funciones utilizadas devuelen el inicio del siguiente simbolo OFDM (CP+Payload), pero en un canal multitrayectoria esto puede presentar problemas de ISI. Cambiar indice devuelto o sumar un factor (p. ej. CP/2) para evitarlo.
- [ ] Revisar indice STO para simbolos OFDM, posibles casos ([[/home/isaac/Documents/Libros/Comunicaciones/OFDM/MIMO-OFDM Wireless Communications with MATLAB (Yong Soo Cho, Jaekwon Kim, Won Young Yang etc.) (z-library.sk, 1lib.sk, z-lib.sk).pdf|MIMO-OFDM Wireless Communications with MATLAB, page 154]])
- [.] Por umbral
- [x] Por derivadas