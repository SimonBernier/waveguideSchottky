# THz Wave Control via Moving Dielectric Fronts - SLIPSTREAM Platform

2D Finite-Difference Time-Domain (FDTD) simulations of terahertz (THz) wave manipulation using relativistic moving fronts in semiconductor-filled waveguides. This work demonstrates temporal stretching and time-reversal of THz pulses through spatiotemporal control of photoexcited carrier density.

## 🔬 Overview

This repository contains the FDTD simulation code used to generate **Figure 2** in our Nature Communications Physics publication on the SLIPSTREAM (Spacetime Light-Induced Photonic STRucturEs for Advanced Manipulation) platform.

The work explores how moving dielectric perturbations - created by photoexciting mobile charge carriers in a semiconductor waveguide - can manipulate THz light in exotic ways including temporal pulse stretching and time-reversal operations.

**Paper:** [Front-induced transitions control THz waves](https://doi.org/10.1038/s42005-021-00667-4)  
A.W. Schiff-Kearn, L. Gingras, **S. Bernier**, J.-M. Ménard, and D.G. Cooke  
*Communications Physics* **4**, 162 (2021)

---

## 🎯 Physical Concept

### The SLIPSTREAM Platform

**Key Innovation:** By tilting the pulse front of a near-infrared pump laser, we create a moving front of photoexcited carriers in a silicon-filled parallel plate waveguide. This front travels at a controllable velocity v<sub>f</sub> relative to the THz wave velocity c/n<sub>Si</sub>.

**Three Regimes:**

1. **Subluminal (v<sub>f</sub> < c/n<sub>Si</sub>)** - THz pulse stretching with quasi-static plateaus
2. **Luminal (v<sub>f</sub> ≈ c/n<sub>Si</sub>)** - Optimal phase-matched emission
3. **Superluminal (v<sub>f</sub> > c/n<sub>Si</sub>)** - Time-reversal via front-induced transitions

### Physical Mechanism

- **THz Generation:** Built-in Schottky fields at metal-semiconductor interfaces
- **Pulse Shaping:** Spatiotemporal modulation via moving photoexcitation front
- **Control Parameter:** Front velocity tuned by optical pump tilt angle
- **Applications:** Sub-cycle THz control, dispersion compensation, pulse engineering

---

## 💻 My Contributions - Simulation Work

### What I Did (Subluminal Regime Simulations) ✅

I developed and executed **2D-FDTD simulations for the subluminal regime** that successfully:

- **Modeled the parallel plate waveguide geometry** with silicon and conducting boundaries
- **Implemented Schottky field emission** at top and bottom metal-semiconductor interfaces
- **Incorporated Drude dispersion model** for photoexcited silicon with realistic carrier densities (~10¹⁷ cm⁻³)
- **Simulated moving carrier density fronts** with velocity v<sub>f</sub> = 0.86 c/n<sub>Si</sub>
- **Reproduced experimental THz waveforms** showing temporal pulse stretching
- **Validated the quasi-static plateau formation** mechanism

**Key Achievement:** The simulations in Figure 2 of the paper quantitatively matched experimental data for various beam clipping configurations, confirming our physical understanding of the subluminal pulse stretching mechanism.

### Technical Implementation

**FDTD Algorithm Details:**
- 2D spatial grid with perfectly conducting boundaries
- Time-stepping with Courant stability condition
- Drude model: carrier scattering time τ = 0.1 ps
- Schottky field depth: ~1 μm (Debye length)
- TEM mode extraction at fixed position
- Post-processing filter for detection response

---

## 🚧 Challenges and Limitations

### Superluminal Regime - Numerical Difficulties

**The Challenge:** While the subluminal simulations worked excellently, I encountered significant numerical challenges when attempting to simulate the superluminal regime (v<sub>f</sub> ≥ c/n<sub>Si</sub>).

**Technical Issues:**
- Standard FDTD algorithms become unstable near or beyond the phase velocity
- Numerical dispersion errors accumulate for relativistic front velocities
- Courant condition violations for fast-moving dielectric perturbations
- Difficulty capturing front-induced transitions at phase-matched conditions

**What Was Needed:**
- Modified FDTD schemes with moving reference frames
- Specialized boundary conditions for superluminal fronts
- Enhanced numerical stability for relativistic regime
- Time-domain formulation of frequency-shifting processes

### Personal Context

Unfortunately, I became seriously ill (cancer diagnosis and treatment) before I could develop and implement the modified FDTD algorithm necessary for the superluminal regime simulations. The paper was published during my recovery period.

**Impact:** The experimental results for time-reversal (superluminal regime) were not independently confirmed numerically in the publication. While the physics is well-supported by theory and experimental data, full numerical validation would have strengthened the complete picture.

---

## 📂 Repository Contents

This repository contains the **working FDTD code for subluminal simulations** that successfully generated Figure 2 of the paper.

**Included:**
- 2D-FDTD solver for THz propagation
- Moving front carrier density profiles
- Schottky field implementation
- Drude dispersion for photoexcited silicon
- Post-processing and filtering

**Not Included:**
- Superluminal regime solver (requires algorithmic modifications)
- Time-reversal simulations
- Full 3D field distributions

---

## 🔗 Broader Context - Moving Front Research

This experimental work on THz optics connects to my subsequent computational studies on quantum systems:

**Common Thread:** Spatiotemporal quenches and moving fronts

- **This work (2021):** THz photonics with moving dielectric fronts
- **Long-range paper (2023):** Quantum quenches in Ising models with power-law interactions
- **2D paper (2025):** Efficient ground state preparation via moving parameter fronts

**Insight:** The concept of using moving fronts to control systems - whether electromagnetic waves or quantum wavefunctions - proved to be a powerful unifying theme across my research.

---

## 📄 Publication

**[Front-induced transitions control THz waves](https://doi.org/10.1038/s42005-021-00667-4)**  
A.W. Schiff-Kearn, L. Gingras, S. Bernier, J.-M. Ménard, and D.G. Cooke  
*Communications Physics* **4**, 162 (2021)  
**Open Access** - Nature Publishing Group

### My Role

**Simulation contributor** - Developed FDTD simulations for subluminal regime (Figure 2), validating experimental pulse stretching mechanism and demonstrating quantitative agreement between theory and experiment.

---

### Software

- MATLAB implementation (this repository)
- Custom 2D-FDTD solver
- Drude dispersion integration
- Spatiotemporal source modeling

---

## 💡 Lessons and Future Work

### What I Learned

✓ FDTD is powerful for electromagnetic problems but requires care near phase transitions  
✓ Subluminal regime is more numerically stable than superluminal  
✓ Good agreement between simulation and experiment validates physical models  
✗ Superluminal regime needs specialized numerical methods  
✗ Standard FDTD breaks down for relativistic moving boundaries

---

## 🔗 Dependencies

- MATLAB (any recent version)
- Standard numerical libraries
- No special toolboxes required

---

## 📧 Contact

**Simon Bernier**
- Email: simon.bernier@mail.mcgill.ca
- LinkedIn: [simon-bernier-6701a9285](https://www.linkedin.com/in/simon-bernier-6701a9285)

---

## 📝 Citation

If you use this code or build upon this work, please cite:

```bibtex
@article{schiffkearn2021front,
  title={Front-induced transitions control THz waves},
  author={Schiff-Kearn, A.W. and Gingras, L. and Bernier, S. and M{\'e}nard, J.-M. and Cooke, D.G.},
  journal={Communications Physics},
  volume={4},
  pages={162},
  year={2021},
  publisher={Nature Publishing Group},
  doi={10.1038/s42005-021-00667-4}
}
```

---

## 🎯 Bottom Line

**What worked:** FDTD simulations successfully validated the subluminal THz pulse stretching mechanism, providing quantitative agreement with experimental observations.

**What didn't:** Superluminal regime simulations required algorithmic modifications I couldn't complete due to illness.

**Why it matters anyway:** The successful subluminal simulations confirmed our physical understanding and enabled publication of a complete experimental story. Sometimes research challenges come from unexpected places, but the work that *did* get done contributed meaningfully to advancing THz photonics.

---

*This project demonstrates: electromagnetic simulation, FDTD methods, ultrafast optics, THz photonics, spatiotemporal control, and the reality that research doesn't always go according to plan - but valuable contributions can still be made.*
