# README

## Functions
The "functions" folder contains the simulation pipeline and read count generation scripts. The pipeline invokes "cosimu," the read count generator, and DESeq2 to produce the simulated signatures.

## Paper scripts
The scripts employed to generate the simulated data from the paper reside in the "paper_scripts" folder.
An asymmetrical distribution was utilized to generate the benchmark data presented in Figure 5.A (associated with Table S4). A symmetrical distribution was applied in Figure 5.B and the supplementary Figure S4.

### Data associated
- `real_asymmetrical_dist.yml`
- `real_symmetrical_dist.yml`

## Simple examples
Here we provide an example of simulation using the YAML parameterization file (example_sim.R). Additionally, example_param.R presents how to generate this type of YAML.

### Data associated
- `example_param.yaml`
- `simple_yaml_structure.yaml` : this file serves as model of the structure that the YAML file should have for easy loading as parameterization for the simulation.

 This serves as the model structure that the YAML file should adhere to for easy loading as parameterization for the simulation.