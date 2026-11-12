![[MovingAverageFilter.png]]
Smoothing of a noisy sine (blue curve) with a moving average (red curve).

The moving average of streaming data is computed with a finite sliding window:

$$
movAvg=\frac{x[n]+x[n-1]+\cdots+x[n-N]}{N+1}
$$

The moving average filter is a **special case of the regular FIR filter**. Both filters have finite impulse responses. The moving average filter uses a sequence of scaled 1s as coefficients, while the FIR filter coefficients are designed based on the filter specifications. They are not usually a sequence of 1s.