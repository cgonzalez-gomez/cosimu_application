# README

## Paper scripts
The scripts employed to generate the simulated data from the paper reside in the "paper_scripts" folder. This scripts aren't necessary updated to be run in external architechtures so errors might appear during their execution. 
- First (00_infer_params.R), two parametrizations were inferred (real_asymmetrical_dist.yml and real_symmetrical_dist.yml) from real data log fold change distribution (GSE185985 - TLN468 treatment). Additionally, the proportional to true number of transcripts per gene was also inferred from real data (GSE185985). 
- Second, the simulation pipeline (simulation_pipeline.R) integrating cosimu, the read counts generation and the DEA, was run for two different parametrizations : asymmetrical (01_sim_asymmetrical.R) and symmetrical (02_sim_symmetrical.R). After the simulation, the data was formatted in order to be exploided in the following step (03_data_formatting.R). All the data could be found in Zenodo repository : [complete].
- Third, run the sbatch commands (04_sbatch_scores_asym.sh and 05_sbatch_scores_sym.sh) to estimate seven connectivity scores from the simulated data, using the scripts in the general_scores folder that call the scores from the coinf package (https://github.com/cgonzalez-gomez/coinf.git). 

The asymmetrical distribution was utilized to generate the benchmark data presented in Figure 5.A (associated with Table S4), and the symmetrical distribution was applied in Figure 5.B and the supplementary Figure S4.

The last file 06_sim_V2_real_data.R presents the correlation analysis done to compare experimental VS simulated data.

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
- `simple_yaml_structure.yaml` : this file serves as model of the structure that the YAML file should have for easy loading as parameterization for the simulation.

 This serves as the model structure that the YAML file should adhere to for easy loading as parameterization for the simulation.