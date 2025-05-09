# Demos of Section 6.7. Multiuser Diversity: System Aspect. 
## This Demo file consists of the following parts. Have fun!

*Author: Seongwook Jung*

- Referred Book: Fundamentals of Wireless Communication (David Tse).
- Last Revised on 2025/05/09.
- If you want to see demos of entire parts altogether, execute `Proportional_fair_scheduler_Demo.m`.


### Part 1. Comparison of Various Schedulers in Asymmetric Environment
- This part compares basic schedulers; **(i) Round-Robin, (ii) Max-Rate Only, (iii) Proportional-Fair Scheduler**.

### Part 2. Proportional Fair Scheduler
- This part simulates PF scheduler and observes the user-selecting behavior.
- For one's taste, please tailor minor parameters such as latency time scale $t_c$ (`LPF`) and `plot range`.

### Part 3. Multiuser Diversity on Various Mobility Scenarios
- For fixed(Rician) and from low-mobility (3km/h Jake's Model, Rayleigh) to high-mobility (70km/h Jake's Model, Rayleigh), we delve into exploitness of multiuser diversity for each channel scenario.
- Jake's approximate version of Clarke's model for implementation is given as $$h(t) \approx \frac{1}{\sqrt{M}} \sum_{l=1}^{M} \exp(j f_D \cos (\alpha_l) t + \phi_l)$$ where $\alpha_l, \phi_l \sim \text{Unif}[0, 2\pi)$
  - Autocorrelation $R(\tau)$: $$R(\tau) = \frac{1}{M} \sum_{l=1}^{M} \exp(j f_D \cos (\alpha_l) \tau) \longrightarrow J_0 (2 \pi f_D \tau)$$
  - Please take note that if $M \rightarrow \infty$ (`M`), it becomes more accurate but requires more computation flops (so tailor this to your purpose)
  - For more detailed information about Jake's model representation of Rayleigh and Rician that includes doppler effect, see https://en.wikipedia.org/wiki/Rayleigh_fading.

### Part 4. Opportunistic Beamforming on Slow Fading Channel
- From Part 3, we have come to know dynamic range of the channel (big fluctuation) induces multiuser diversity for static environment.
- Thus let we analyze the performance of opportunistic beamforming and discuss/compare with the tx beamforming (coherent beamforming) techniques in *Slow* Fading Channel.

### Part 5. Opportunistic Beamforming on Fast Fading Channel
- The natural question that arises is that does opportunistic scheme is useless when it comes to **Fast** Fading channel?
- Generally yes, but if for specific scenarios such that the dynamicity can be added although the channel is already fast-fading, it impinges upon the channel via adding additional fluctuation. Take a look at this.

I welcome any kind of suggestion or revision requests.
