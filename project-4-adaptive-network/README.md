# Project IV: An Adaptive Thalamus–Subplate–Layer 4 Network (Depressing Synapses and STDP)

Course project for Computational Neuroscience (EC60007), IIT Kharagpur, Autumn 2026. Implemented in MATLAB.

The project simulates a small spiking network in which two thalamic inputs, a **standard (S)** and a **deviant (D)** stimulus channel, drive a single **subplate (SP)** neuron and a single **layer 4 (L4)** neuron, and SP also projects to L4. All five synapses are short-term depressing, and the three synapses onto L4 are also plastic through spike-timing-dependent plasticity (STDP). The questions are how the network responds to an oddball stimulus stream, and how its synaptic weights develop over an hour of stimulation.

![Oddball responses of SP and L4](results/Figure_3.png)

## What is in this folder

```
project-4-adaptive-network/
├── src/project4_adaptive_network.m   # single script, produces every figure
└── results/Figure_1.png ... Figure_17.png, Figure_30.png (Figures 9, 10 and 18-29 are not used)
```

Run `src/project4_adaptive_network.m` in MATLAB. Each question is a `%%` section. The two 60-minute runs (Q4 and Q5) simulate 3.6 million 1 ms steps each, in chunks of 100 stimuli, and write large `.mat` files (`q4_results.mat`, `q5_results.mat` and state files) next to the script. The Q4 section that follows (the `verNOp` runs) loads `q4_results.mat`. Those files are excluded from the repository by `.gitignore`.

## Model

**Synapses (depressing).** Each synapse has a recovered, effective and inactive fraction `xr, xe, xi`:

```
dxr/dt = −M·xr/τre + xi/τir
dxe/dt = +M·xr/τre − xe/τei
dxi/dt = xe/τei − xi/τir
```

with `M` the presynaptic spike and `τre = 0.9 ms`, `τir = 5000 ms`, and `τei = 10 ms` for the thalamic synapses and `27 ms` for SP → L4. The effective fraction `xe` scales the synaptic weight. The equations are integrated with RK4 at 1 ms steps.

**Neurons (leaky integrate-and-fire with refractory drop).** The input to each neuron is `g(t)·w·xe` summed over its synapses, where `g(t)` is the presynaptic spike train convolved with an exponential kernel (`τsyn = 10 ms`) after a 1 ms synaptic delay. The voltage keeps 90% of its previous value at each step (the 10% leak). A spike occurs when the voltage crosses 0.05. A drop `5·exp(−t/2 ms)` is then subtracted from the voltage for 20 ms (absolute and relative refractory period). The refractory term is applied to the value compared with the threshold and is not stored in the integrated voltage.

**Weights.** Th → SP: 0.2 each, fixed. Th → L4: 0.02 each, range 0.0001–0.4. SP → L4: 0.11, range 0.0001–0.11. A weight at a limit cannot cross it but can move away from it.

**STDP on L4 synapses.** Multiplicative updates with `Δt = t_post − t_pre`:

```
Δw/w = +0.015·exp(−Δt/13 ms)    (Δt > 0, potentiation)
Δw/w = −0.021·exp( Δt/20 ms)    (Δt < 0, depression)
```

Pairs are counted all-to-all, using exponentially decaying spike traces so the cost per time step is constant.

**Input spike trains.** Inhomogeneous Poisson processes with rate vectors at 1 ms resolution. During stimulus S the S neuron fires at 10 sp/s and D at 2.5 sp/s (and the reverse for stimulus D); between stimuli both fire at 0.5 sp/s. Stimuli last 50 ms with a 250 ms gap.

## Code structure

| Function | Purpose |
| --- | --- |
| `generate_poisson` | Inhomogeneous Poisson spike train from a rate vector |
| `generate_rates` | Rate vectors for S and D from a stimulus sequence |
| `sp_neuron` | SP neuron, with its two thalamic depressing synapses |
| `l4_neuron` | L4 neuron, with three depressing synapses and STDP (switchable with a `plastic` flag) |
| `rk4` | One RK4 step of the synaptic depression equations |

`sp_neuron` and `l4_neuron` take an optional state structure (voltage, filtered inputs, depression variables, refractory tail, weights and STDP traces) and return the updated state, so a long run can be done in chunks with exactly the same result as a single run. `plastic = 0` gives the version without long-term plasticity (called `verNOp` in the handout).

## Questions and results

| Question | What the code does | Figures |
| --- | --- | --- |
| Q1 | Inhomogeneous Poisson generator for a random rate around 40–50 sp/s. PSTH from 10 to 640 repeats against the true rate, and the RMSE against the number of repetitions. The RMSE falls from about 64 sp/s at 10 repeats to about 11 sp/s at 320 and about 8 sp/s at 640, roughly as `1/√n`. | 1–2 |
| Q2 | Oddball stream (14 S, one D as the 8th stimulus), 50 repeats, no plasticity: PSTHs of SP and L4 (10 ms bins). Both responses depress over the S stimuli, and the response to D is larger than the S response just before it. | 3 |
| Q3 | The same protocol with `τir` = 1000, 3000 and 10000 ms. Shorter `τir` gives weaker adaptation and a less distinct D response; longer `τir` gives strong adaptation, especially in L4. | 4–6 |
| Q4 | 60 minutes of stimulation with S and D at probabilities 0.9 and 0.1, with plasticity. The weights are tracked over time, and `verNOp` is run at time points that the script detects from the weight curves (start, first rise of the thalamic weights, start of the SP → L4 decline, saturation of D → L4 and of S → L4, SP → L4 reaching its floor, and 5 minutes after that) and marks on the weight figure. A summary figure gives the spike count per trial in the 100 ms after each stimulus onset, for every time point. | 7, 11–17, 30 |
| Q5 | The same 60 minutes with S and D equally likely. | 8 |

### Weight development (Q4)

![Synaptic weights over 60 minutes, p(D) = 0.1](results/Figure_7.png)

In the run shown, D → L4 and S → L4 both start at 0.02 and are potentiated, D → L4 first (it reaches 0.4 at about 5.6 minutes) and S → L4 later (about 7.3 minutes). SP → L4 starts at its upper limit (0.11), stays there until the thalamic weights have grown, and then falls to its lower limit by about 9.5 minutes. After that the weights stay at their limits for the rest of the hour.

A likely reason for the pruning of SP → L4 is that, once the thalamic synapses are strong, L4 tends to fire before SP does, so STDP depresses the SP synapse. This was not tested directly.

The time points, and the weights used at each of them, are printed in the title of each `verNOp` figure. Which of S → L4 and D → L4 leads depends on the random stimulus sequence of the run, and the time points change from run to run.

### Oddball responses after learning (Q4)

The `verNOp` network is run with the weights frozen at each detected time point (Figures 11–17). To compare the time points without reading single PSTH peaks, which are noisy with 50 repeats, Figure 30 plots the spike count per trial in the 100 ms after each stimulus onset against the stimulus number (red line: the deviant, stimulus 8).

![Spikes per trial against stimulus number, one line per time point](results/Figure_30.png)

For the run shown (time points 0.50, 2.59, 5.64, 5.69, 7.35, 9.46 and 14.46 minutes):

- **Start (0.50 and 2.59 min, weights close to the initial ones).** L4 adapts strongly, from about 5 spikes per trial at stimulus 1 to well below 1 by stimulus 5. The deviant response is small.
- **Just after D → L4 has been potentiated (5.64 and 5.69 min, D → L4 ≈ 0.38–0.40, S → L4 ≈ 0.11, SP → L4 ≈ 0.10).** The L4 response to the deviant is clearly above that to the preceding S, and so is the SP response. The enhanced response here comes from D → L4 being strong while S → L4 is not, not only from the D synapse being fresh. This is the clearest case of deviance detection after learning.
- **Both thalamic weights saturated and SP → L4 near zero (7.35 min onward, last panel 14.46 min).** L4 is driven directly by the thalamus, and its responses to S stay high (about 1 to 2 spikes per trial) instead of adapting away. The deviant response no longer stands out from its neighbours.

Counts at single stimuli are noisy (50 repeats), so differences of a few tenths of a spike between neighbouring stimuli should not be over-interpreted.

![Oddball responses at the last time point](results/Figure_17.png)

### Equal stimulus probabilities (Q5)

![Synaptic weights over 60 minutes, p(D) = 0.5](results/Figure_8.png)

When S and D are equally likely, the two thalamic synapses on L4 follow similar trajectories and saturate together, while SP → L4 is again pruned. In the p(D) = 0.1 run (Figure 7), D → L4 leads, so the two runs can be compared directly.

## Figure index

| Figure | Content |
| --- | --- |
| 1 | PSTH estimation error versus number of repeats (10 to 640) |
| 2 | PSTH from 640 repeats against the true rate λ(t) |
| 3 | Oddball PSTHs of SP and L4 (initial weights, τir = 5000 ms) |
| 4, 5, 6 | Oddball PSTHs for τir = 1000, 3000 and 10000 ms |
| 7 | Synaptic weights over 60 minutes, p(D) = 0.1, with the detected time points marked |
| 8 | Synaptic weights over 60 minutes, p(D) = 0.5 |
| 11–17 | Oddball PSTHs (no plasticity) with the weights frozen at the detected time points of the Q4 run, one figure per time point |
| 30 | Spikes per trial in the 100 ms after each stimulus onset, against stimulus number, for SP and L4 and every time point |

## Notes

- Results depend on the random spike trains and stimulus sequence, so the exact timings of the weight changes differ from run to run.
- The assignment handout is not part of this repository.
