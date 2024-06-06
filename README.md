# README

## Paper scripts
The scripts used to generate the simulated data from the paper are located in the "paper_scripts" folder. These scripts may not be updated to run on external architectures, so errors may occur during execution.

### 1 - Inference of Parameters
Firstly (00_infer_params.R), two parametrizations were inferred (real_asymmetrical_dist.yml and real_symmetrical_dist.yml) from the log fold change distribution of real data (GSE185985 - TLN468 treatment). Additionally, the number of transcripts per gene, proportional to the true values, was inferred from real data (GSE185985).

### 2- Comparison of real data VS simulated data
Secondly, we compared the properties between real and simulated data. This comparison was done in two stages. On one hand, we analyzed real data biological replicates (01_real_data_replicates_analysis.R) and then compared correlation between experimental replicates VS the correlation between simulated signatures (02_1_sim_VS_real.R). On the other hand, we compared the log2 fold-change distribution of real and simulated signatures (02_2_sim_VS_real.R).

### 3- Simulation Pipeline Execution
Third, the simulation pipeline (simulation_pipeline.R), which integrates cosimu for read counts generation and DEA, was executed for two different parametrizations: asymmetrical (03_sim_asymmetrical.R) and symmetrical (04_sim_symmetrical.R). After the simulation, the data was formatted for further analysis (03_data_formatting.R). Formatted data can be found in the Zenodo repository: https://zenodo.org/records/11195795.

### 4- Connectivity Score Estimation
Then, the sbatch commands were run (04_sbatch_scores_asym.sh and 05_sbatch_scores_sym.sh) to estimate seven connectivity scores from the simulated data. These scores were calculated using scripts from the general_scores folder, which call upon scores from the coinf package (https://github.com/cgonzalez-gomez/coinf.git).

The asymmetrical distribution was utilized to generate the benchmark data presented in Figure 5.A (associated with Table S4), and the symmetrical distribution was applied in Figure 5.B and the supplementary Figure S4.

Finally, the benchmark results were summarized by calculating and plotting the average precision score (08_avg_precision_score.R). The regression plot (Figure 6) was generated with the script 09_scores_regression.R.

### 5- Cosimu properties
The simulation to analyze the properties of the cosimu transition parameters, as well as the visualization of the results (Figure 3), was done using the script 10_cosimu_properties.R.

#### Data associated
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