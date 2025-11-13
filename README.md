# Hexatic Order (ψ₆) Correlation Analysis Pipeline

A high-performance, modular C pipeline for computing the hexatic order parameter (ψ₆) and its spatial correlation function (g₆(r)) from 2D simulation snapshots.



---

## 🚀 Quick Start

1.  **Clone:**
    ```bash
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    cd your-repo-name
    ```
2.  **Compile:**
    ```bash
    make
    ```
3.  **Run:**
    ```bash
    ./hexatic_g6_avg ./data/ 1000 1200 ./out/
    ```

---

## 🌟 Overview

This project provides a complete analysis pipeline, starting from raw particle coordinate files (`time_*.dat`) and ending with a single, averaged `g6_avg.dat` file.

It is designed for analyzing 2D systems (like active matter or molecular simulations) to determine the extent and range of hexatic (six-fold) order. The entire codebase is modular, so you can easily reuse or replace components.

### Core Pipeline

1.  **Read:** Ingests particle `(x, y)` positions from snapshot files.
2.  **Cluster:** Groups nearby particles into clusters using a fast Union-Find algorithm.
3.  **COM:** Calculates the center-of-mass (COM) for each cluster.
4.  **Neighbors:** Identifies all "natural" neighbors for each COM using **Delaunay Triangulation**.
5.  **Compute ψ₆:** Calculates the complex hexatic order parameter (ψ₆) for every COM based on its neighbors.
6.  **Accumulate g₆(r):** Computes all pair correlations `psi6(i) * conj(psi6(j))` and adds them to the correct radial distance bin.
7.  **Average:** After processing all snapshots, it divides by the total pair counts to output the final, averaged `g6(r)` function.

---

## ✨ Key Features

* **Geometrically Correct:** Uses the **Delaunay Triangulation** (via the `Triangle` library) to find true neighbors, not a k-NN approximation.
* **Modular by Design:** Each step (I/O, clustering, psi6, etc.) is in its own `.c`/`.h` module.
* **PBC-Aware:** Correctly handles Periodic Boundary Conditions for clustering, COM calculations, and neighbor finding.
* **Robust Averaging:** Uses pair-weighted averaging to ensure snapshots with more clusters contribute correctly.
* **Efficient:** Written in C with a focus on memory and speed.

---

## 🛠️ Build & Run

### Dependencies

* A C compiler (e.g., `gcc` or `clang`).
* The `make` build tool.
* The `Triangle` library (included as `triangle.c` and `triangle.h` in the repo).

### Build

The included `Makefile` handles all compilation.

* **Standard Release Build:**
    ```bash
    make
    ```
* **Debug Build (slower, with checks):**
    ```bash
    make DEBUG=1
    ```
* **Clean up compiled files:**
    ```bash
    make clean
    ```

### Run

Run the program from the command line, providing the key parameters.

**Usage:**
```bash
./hexatic_g6_avg DATA_DIR START_INDEX END_INDEX OUTPUT_DIR [OPTIONS]
