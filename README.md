# MP-BioPath Integration Pipeline

**Author:** Julian Benson (B.S. Biochemistry & Molecular Biology)
**Target:** Reactome / Open Genome Informatics (GSoC 2026)

## Project Overview
This repository contains the data ingestion and optimization bridge designed to transition [MP-BioPath](https://github.com/OICR/mp-biopath) from a literature-validation tool to a high-throughput predictive engine for genomic data. 

The pipeline utilizes a polyglot architecture:
* **Python/R (`src/module_a_ingestion`):** Processes raw RNA-seq differential expression data into the continuous 0.01 to 100 boundaries required by the MP-BioPath non-linear mathematical models.
* **Python (`src/module_b_graphs`):** Interfaces with the Reactome REST API to dynamically stitch sub-pathway logic graphs.
* **Julia/JuMP (`src/module_c_solver`):** Modifies the core objective function weights to resolve entity-set dilution and automate the 15% biological significance threshold.

## Repository Structure
* `/data`: (Ignored via `.gitignore`) Local staging for MCB112 quantitative datasets and RNA-seq counts.
* `/src`: Core module scripts for ingestion, API translation, and JuMP solver execution.
* `/docs`: Project proposals, mathematical documentation, and architecture blueprints.

## References
Wright, A.J., Orlic-Milacic, M., Rothfels, K. et al. Evaluating the predictive accuracy of curated biological pathways in a public knowledgebase. *Database* (2022). DOI: 10.1093/database/baac009