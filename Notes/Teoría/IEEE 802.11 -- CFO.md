Ver tambien: [[IEEE 802.11 PHY]]
Referencia: [[/home/isaac/Documents/Cinvestav/Tesis/Thesis References/frequency-offset-estimation-and-correction-in-the-ieee-80211a-wl.pdf|Souror et al]]

### CFO representation
$\alpha$ represents the **rotation per sample** due to CFO
$$
\alpha = 2\pi \Delta_f T_s
$$
where $\Delta_f$ is the frequency error between the received carrier and the LO.
To model the effect of CFO, sample n of the received frame is multiplied by $e^{jn\alpha}$.

## Estimation Using the Short Training (Coarse estimation)
![[IEEE 802.11 frame structure.png]]

With a short training consisting of 10 repeated parts with 16 complex-value samples each part ([[IEEE 802.11 PHY]]) , lets denote de ST complex samples of the last 5 parts as $S_m$, where $m=0,1,\dots 79$. 

The autocorrelation for a window of 16 samples:
$$
P = \sum_{m=0}^{63} S^*_m S_{m+16}
$$
has an angle of:
$$
\angle P = 16\alpha_{ST}
$$
See [[Schmidl & Cox — CFO#4. Correlación entre las dos mitades|Angle of an autocorrelation process]] for its demonstration, even though it is obvious from the definition of $\alpha$.

Estimating $\alpha_{ST}$:
$$ \boxed{
\hat\alpha_{ST} = \frac{1}{16} \angle P
}$$
Due to the restriction $-\pi \lt \angle(.) \le \pi$ :
$$-\pi \lt 16\alpha \le \pi$$
$$\boxed{
-\frac{\pi}{16} \lt \alpha \le \frac{\pi}{16}
}$$
and avoid ambiguity.

Ambiguity over subcarrier spacing:
$$-\frac{\pi}{16} \lt 2\pi \Delta_f T_s \le \frac{\pi}{16}$$
$$-\frac{\pi}{16} \lt 2\pi (\frac{\epsilon}{NT_s}) T_s \le \frac{\pi}{16}$$
With $N=64$ ([[IEEE 802.11 PHY]]):
$$-\frac{\pi}{16} \lt 2\pi \frac{\epsilon}{64} \le \frac{\pi}{16}$$
$$\boxed{-2 \lt \epsilon \le 2}$$

## Estimation Using the Long Training (Fine estimation)
The two identical parts $T_0$ and $T_1$ of LT have 128 samples denoted here as $S_m$, $m=0,1,\dots 127$.

Calculating as same as the short training:
$$
P = \sum_{m=0}^{63} S^*_m S_{m+64}
$$
$$
\angle P = 64\alpha_{LT}
$$
$$ \boxed{
\hat\alpha_{LT} = \frac{1}{64} \angle P
}$$
$$\boxed{
-\frac{\pi}{64} \lt \alpha \le \frac{\pi}{64}
}$$
$$-\frac{\pi}{64} \lt 2\pi \frac{\epsilon}{64} \le \frac{\pi}{64}$$
$$\boxed{-\frac{1}{2} \lt \epsilon \le \frac{1}{2}}$$
## Estimation and Correction Process with Training Symbols
1. Coarse estimation $\alpha_{ST}$.
2. Coarse correction: The 128 samples of LT are corrected as $S_m = S_m e^{-jm \hat\alpha_{ST}}$, $m=0,1\dots 127$.
3. Fine estimation $\alpha_{LT}$.
4. Fine correction: The samples are corrected one more time as $S_m = S_m e^{-jm \hat\alpha_{LT}}$, $m=0,1\dots 127$.

The final phase estimation $\hat\alpha = \hat\alpha_{ST} + \hat\alpha_{LT}$ is used for carrier frequency correction for the remaining part of the frame after LT.

## Estimation and Correction Process with Pilot Symbols
### Channel estimation
The parts $T_0$ and $T_1$ are known by the receptor, and are used for channel estimation. We denote the channel gain estimates as $H_k$, $k=-26 \dots 26$. The Pilot Symbols are given in the subcarrier positions $k=\{-21,-7,7,21\}$. We represent the **channel gains at the Pilot Symbols locations** as:
$$Q_{-2}=H_{-21}$$
$$Q_{-1}=H_{-7}$$
$$Q_{1}=H_{7}$$
$$Q_{2}=H_{21}$$

At the receiver, the FFT output is denoted as $X_{l,k}$, where $k$ is the subcarrier number, $k=-26 \dots 26$, and $l$ is the OFDM Symbol numbre starting from $l=1$. For simple notation, **Pilot Symbols at the l'th OFDM symbol** are referred as:
$$P_{l,-2}=X_{l,-21}$$
$$P_{l,-1}=X_{l,-7}$$
$$P_{l,1}=X_{l,7}$$
$$P_{l,2}=X_{l,21}$$

### Residual CFO model
In the presence of **residual frequency offset and sampling frequency offset** the FFT outputs, including the Pilot Symbols, will incur phase rotation. The Pilot Symbols $P_{l,i}$  and the corresponding channel gains $Q_{i},i=-2,-1,1,2$, are used to estimate this offset. It is shown that due to frequency offset, each FFT output $X_{l,k}$ incurs extra phase rotation $\varphi_{l,k}$ given by:
$$
\varphi_{l,k}=2\pi l \left(\epsilon_r T_u f_c - \frac{N+N_g}{N} \epsilon k \right)
$$
where
- $N=64$ is the FFT size
- $N_g=16$ is the guard samples length
- $T_u=4$ usec is the OFDM Symbol duration
- $\epsilon_r = \epsilon - \hat\epsilon_0$ is the residual carrier frequency offset after LT correction, which is nonzero if $\hat\alpha \ne \alpha$

The **first part in the right hand side** of the equation is due to the residual carrier frequency offset not corrected after LT correction. This term is equal for all subcarriers of the same OFDM Symbol.
The **second part** is due to the continuous drift in the sampled waveform of the OFDM Symbol due to sampling offset. This phase depends linearly on the subcarrier number.

### Steps to correct offsets
#### 1. Correct Sampling Frequency Offset
we first compensate for the second term in the right hand side of the equation. Using the available frequency offset estimate $\hat\epsilon_{l-1}$ the l’th OFDM Symbol FFT output is updated as:
$$X_{l,k}=X_{l,k} \ exp  \left(j 2\pi l \ \frac{N+N_g}{N} \ \hat\epsilon_{l-1}k \right) \quad ,k=-26\dots 26$$
Note that this correction includes the Pilot Symbols. When processing the Signal field (l=1) the value of $\hat\epsilon_0$ is the one generated from the LT processing. As more OFDM Symbols are received the frequency offset estimate $\hat\epsilon_l$ improves.
#### 2. Correct Residual Carrier Frequency Offset
First term in the right hand side of the equation is due to the accumulated phase rotation, at the l’th OFDM Symbols, when $\hat\alpha \ne \alpha$. A good estimate of it is given by:
$$\hat\beta_l = \angle \left( \sum_{i=-2,-1,1,2} P_{l,i} \ Q^*_i \right)$$
This is the correlation between the Pilot Symbols, corrected in the previous step, and the channel gain at the corresponding subcarrier. Again, see [[Schmidl & Cox — CFO#4. Correlación entre las dos mitades|Angle of an autocorrelation process]] if it is not clear.
Now the output of the previous step is corrected one more time as follows:
$$X_{l,k}=X_{l,k} \ e^{-j\hat\beta_l} \quad ,k=-26\dots 26, \quad k\ne -21,-7,7,21$$
The Pilot Symbols are not corrected since they will be used in the next step for finer estimation of the frequency offset estimate $\hat\epsilon_l$ .
#### 3. Fine estimation of Frequency Offset
The goal is to improve $\hat\epsilon_l$ to be used in step 1 for the following OFDM Symbol.

We first calculate the variable $W_l$ as follows:$$W_l=\sum_{i=-2,-1,1,2} P_{l,i} \ P^*_{l-1,i}$$We notice that $\angle(W_l)$ is an estimate of the phase increment $2\pi \epsilon_r T_u f_c$ over one OFDM symbol, due to residual carrier frequency offset. 
To improve the accuracy of estimation we could continuously accumulate $W_l$ every new OFDM symbol. However, **continuous addition would be problematic in a fixed-point hardware implementation**. As a trade off we add groups of 4 successive values and exponentially filter the sum.

Addition: $$V_l = \sum^l_{m=l-3} W_m$$Filtering:$$U_l = \rho V_l + (1-\rho) U_{l-4}$$The parameter $\rho$ is taken to be 1/32. 
Every 4'th OFDM Symbol a filtered value $U_l$ is generated, and a corresponding new estimates $\epsilon_r$ and $\hat\epsilon_l$ are updated  as follows:
$$\hat\epsilon_r = \angle(U_l)/(2\pi T_u f_c) \quad ,l=4,8,12,16,\dots$$
$$\hat\epsilon_{l+1} = \hat\epsilon_{l+2} = \hat\epsilon_{l+3} = \hat\epsilon_{l} = \hat\epsilon_{0} + \hat\epsilon_{r} \quad ,l=4,8,12,16,\dots$$ 