# RF Matching Network Design (3 Designs)

This repository contains MATLAB code and supporting materials for designing and analyzing RF matching networks. The project focuses on three classical matching techniques

1. **L-Matching Network**
2. **Single-Stub Shunt Tuning**
3. **Single-Stub Series Tuning**

The goal is to perform impedance matching between a load and a transmission line to minimize signal reflection and maximize power transfer.

---

## 📌 Objectives

- To design and implement three types of impedance matching networks.
- To compute and plot the **reflection coefficient (Γ)** over frequency.
- To determine and compare the **fractional bandwidth (FBW)** of each design.
- To evaluate which design offers better bandwidth performance.

---

## 🧠 Theoretical Background

Impedance mismatch in RF systems causes signal reflection, reducing power transfer efficiency. This project applies classical matching techniques to address this:

- **L-Matching Network** uses lumped LC components for narrowband matching.
- **Stub Tuning** uses sections of transmission line (open or shorted) to cancel reactive components.

---

## 📈 Key Metrics

For each design:
- Plot of \ (|Γ|) vs frequency.
- Fractional Bandwidth calculated 
- Evaluation against \( \|Γ| = 0.2 \) threshold.

---

## 🛠️ How to Use

1. Open each folder in MATLAB.
2. Run the `.m` file (e.g., `design1_plot.m`) to generate plots.
3. Adjust parameters (e.g., \( Z_L \), \( Z_0 \), frequency range) as needed for your test cases.

---

## 📚 References

- Watcharapan Suwansantisuk, *EIE/ENE 450: Applied Communications and Transmission Lines*, KMUTT.
- Lecture notes on impedance matching and MATLAB simulations.

---

## 🧑‍💻 Author

This project is part of coursework for EIE/ENE 450 — Applied Communications and Transmission Lines at KMUTT.

---

## 📃 License

MIT License or Educational Use Only (you may specify).
