# 📊 seqMeta Automated Analysis Workflow

A streamlined workflow for performing gene-based association analysis using **seqMeta**, starting from VCF, GEN, or BGEN genotype files.

This framework supports distributed cluster execution, local multiprocessing, and resumable analysis for large-scale sequencing studies.

---

## 🧰 System & Software Dependencies

Ensure the following tools are installed prior to running the pipeline:

### Required Software

- Python 2.7 or Python 3.x
- pandas ≥ 0.23.4
- DRMAA ≥ 0.7.9
- pybgen ≥ 0.7.0
- R ≥ 3.2.3
- seqMeta (R package) ≥ 1.6.7

### Optional (for distributed or high-throughput runs)

- LSF / bsub workload manager

---

## ⚙️ Installation Steps

### Install Required Python Modules

