# Serotonin data

This simulated data mimics a PET imaging study where the 5-HT2A receptor
and serotonin transporter (SERT) binding potential has been quantified
into 8 different regions. The 5-HT2A cortical regions are considered
high-binding regions measurements. These measurements can be regarded as
proxy measures of the extra-cellular levels of serotonin in the brain

|       |         |                                                  |
|-------|---------|--------------------------------------------------|
| day   | numeric | Scan day of the year                             |
| age   | numeric | Age at baseline scan                             |
| mem   | numeric | Memory performance score                         |
| depr  | numeric | Depression (mild) status 500 days after baseline |
| gene1 | numeric | Gene marker 1 (HTR2A)                            |
| gene2 | numeric | Gene marker 2 (HTTTLPR)                          |
| cau   | numeric | SERT binding, Caudate Nucleus                    |
| th    | numeric | SERT binding, Thalamus                           |
| put   | numeric | SERT binding, Putamen                            |
| mid   | numeric | SERT binding, Midbrain                           |
| aci   | numeric | 5-HT2A binding, Anterior cingulate gyrus         |
| pci   | numeric | 5-HT2A binding, Posterior cingulate gyrus        |
| sfc   | numeric | 5-HT2A binding, Superior frontal cortex          |
| par   | numeric | 5-HT2A binding, Parietal cortex                  |

## Format

data.frame

## Source

Simulated
