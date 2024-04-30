# README

## Paper scripts
The scripts used to generate the simulated data from the paper are located in the "paper_scripts" folder. These scripts may not be updated to run on external architectures, so errors may occur during execution.

### 1 - Inference of Parameters
Firstly (00_infer_params.R), two parametrizations were inferred (real_asymmetrical_dist.yml and real_symmetrical_dist.yml) from the log fold change distribution of real data (GSE185985 - TLN468 treatment). Additionally, the number of transcripts per gene, proportional to the true values, was inferred from real data (GSE185985).

### 2- Simulation Pipeline Execution
Secondly, the simulation pipeline (simulation_pipeline.R), which integrates cosimu for read counts generation and DEA, was executed for two different parametrizations: asymmetrical (01_sim_asymmetrical.R) and symmetrical (02_sim_symmetrical.R). After the simulation, the data was formatted for further analysis (03_data_formatting.R). All data can be found in the Zenodo repository: [complete].

### 3- Connectivity Score Estimation
Thirdly, the sbatch commands were run (04_sbatch_scores_asym.sh and 05_sbatch_scores_sym.sh) to estimate seven connectivity scores from the simulated data. These scores were calculated using scripts from the general_scores folder, which call upon scores from the coinf package (https://github.com/cgonzalez-gomez/coinf.git).

The asymmetrical distribution was utilized to generate the benchmark data presented in Figure 5.A (associated with Table S4), and the symmetrical distribution was applied in Figure 5.B and the supplementary Figure S4.

Finally, the file 06_sim_V2_real_data.R presents the correlation analysis conducted to compare experimental data with simulated data.

### Data associated
- `real_asymmetrical_dist.yml`
- `real_symmetrical_dist.yml`
- `real_lfc_data.csv`
- `ptt_real.RDS`
- `gene_width.RDS`

## Simple examples
Here we provide an example of simulation using the YAML parameterization file (example_sim.R). Additionally, example_param.R presents how to generate this type of YAML.

### Data associated
- `example_param.yaml`
- `simple_yaml_structure.yaml` : this file serves as a model of the structure that the YAML file should have for easy loading as parameterization for the simulation.