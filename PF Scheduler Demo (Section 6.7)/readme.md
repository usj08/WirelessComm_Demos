# Demos of Section 6.7. Multiuser Diversity: System Aspect. 


## This Demo file consists of the following parts. Have fun!

*Author: Seongwook Jung*

*Referred Book: Fundamentals of Wireless Communication (David Tse).*

*Last Revised on 2025/05/08.*

(1) If you want to see demos of entire parts altogether, execute `Proportional_fair_scheduler_Demo.m`.
(2) Else, execute the parts of interest.

### Part 1. Comparison of Various Schedulers in Asymmetric Environment
- This part compares basic schedulers; **(i) Round-Robin, (ii) Max-Rate Only, (iii) Proportional-Fair Scheduler**.

### Part 2. Proportional Fair Scheduler
- This part simulates PF scheduler and observes the user-selecting behavior.
- For one's taste, please tailor minor parameters such as latency time scale $t_c$ (`LPF`) and `plot range`.

### Part 3. Multiuser Diversity on Various Mobility Scenarios
- For fixed(Rician) and from low-mobility (3km/h Clarke's Model, Correlated Rayleigh) to high-mobility (70km/h Rayleigh), we delve into exploitness of multiuser diversity for each channel scenario.

### Part 4. Opportunistic Beamforming on Slow Fading Channel
- From Part 3, we have come to know dynamic range of the channel (big fluctuation) induces multiuser diversity for static environment.
- Thus let we analyze the performance of opportunistic beamforming and discuss/compare with the tx beamforming (coherent beamforming) techniques in *Slow* Fading Channel.

### Part 5. Opportunistic Beamforming on Fast Fading Channel
- The natural question that arises is that does opportunistic scheme is useless when it comes to **Fast** Fading channel?
- Generally yes, but if for specific scenarios such that the dynamicity can be added although the channel is already fast-fading, it impinges upon the channel via adding additional fluctuation. Take a look at this.

I welcome any kind of suggestion or revision requests.
