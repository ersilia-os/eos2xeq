# Antibiotic downselection criteria based on similarity to known antibiotics

Predictive and generative efforts to discover antibiotics require post-filtering of compounds to ensure novelty and synthetic accessibility of the molecules. In a study to find de novo antibiotics against N. gonorrhoeae or S. aureus, filters for dissimilarity to ca. 500 known antibiotics, as well as (retro)synthetic analysis combined with structural filtering, was applied to downselect compounds. Here, we do not add the filters for activity, cytotoxicity, and (retro)synthetic accessibility, which were also part of the original study.

This model was incorporated on 2025-09-17.Last packaged on 2025-09-22.

## Information
### Identifiers
- **Ersilia Identifier:** `eos2xeq`
- **Slug:** `antibiotics-downselection`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Property calculation or prediction`
- **Biomedical Area:** `Antimicrobial resistance`
- **Target Organism:** `Staphylococcus aureus`, `Neisseria gonorrhoeae`
- **Tags:** `Antimicrobial activity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `7`
- **Output Consistency:** `Fixed`
- **Interpretation:** This model returns the similarity to known antibiotics (at 0.5 Tanimoto cutoff), as as well as the presence of some antibiotic motifs.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| has_pains | integer | high | The molecule has PAINS alert matches |
| has_brenk | integer | high | The molecule has Brenk filter alert matches |
| is_sim_known_ab | integer | high | The molecule is similar to at least one of 500+ known antibiotics using a Tanimoto similarity of 0.5 |
| nitrofuran_motif | integer | high | The nitrofuran motif is found in the molecule using substructure matching |
| fluoroquinolone_motif | integer | high | The fluoroquinolone motif is found in the molecule using substructure matching |
| carbepenem_motif | integer | high | The carbepenem motif is found in the molecule using substructure matching |
| betalactam_motif | integer | high | The beta-lactam motif is found in the molecule using substructure matching |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos2xeq](https://hub.docker.com/r/ersiliaos/eos2xeq)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos2xeq.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos2xeq.zip)

### Resource Consumption
- **Model Size (Mb):** `1`
- **Environment Size (Mb):** `549`
- **Image Size (Mb):** `495.86`

**Computational Performance (seconds):**
- 10 inputs: `27.48`
- 100 inputs: `28.29`
- 10000 inputs: `426.97`

### References
- **Source Code**: [https://github.com/aartikrish/de-novo-antibiotics/](https://github.com/aartikrish/de-novo-antibiotics/)
- **Publication**: [https://www.cell.com/cell/abstract/S0092-8674(25)00855-4](https://www.cell.com/cell/abstract/S0092-8674(25)00855-4)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2025`
- **Ersilia Contributor:** [miquelduranfrigola](https://github.com/miquelduranfrigola)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos2xeq
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos2xeq
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
