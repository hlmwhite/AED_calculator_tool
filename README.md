# PYTHON COURSE - WORK IN PROGRESS

# AED calculator tool

Gene models can be assessed using mapped transcript or protein evidence using a metric called annotation edit distance (AED). The theory behind this can be found at https://link.springer.com/article/10.1186/1471-2105-10-67. 

This can be a useful measure to assess gene predictor accuracy. Here you can find a simple tool for calculating AED per gene model, given some evidence.

# Installation

```bash
git clone https://github.com/hlmwhite/AED_calculator_tool.git
```

## Usage

For single use cases (i.e. one gene model one evidence model), simply run from the command line:

```bash

python3 example_calculator_v2.py -A genome.gff -E evidence.gff -a "<annotation gene>" -e "<evidence model to assess AED>"

```

where genome.gff (-A) is your annotation, evidence.gff (-E) is something like a stringtie output (gff/gtf), a gene ID of interest in your annotation (-a) and finally the transcript ID (-e) to calculate AED for against the gene ID in -a.

An example using the test files included here is:

```bash

python3 example_calculator_v2.py -A test.genome.gff -E stringtie.test.gff -a "rna-XM_034815643.1" -e "STRG.95.1"

## expected output
Number of overlapping values (true positives): 622
Number of non-overlapping values (false negatives): 34
SN (or sensitivity): 0.948170731707317
Number of overlapping values: 622
SP (or specificity): 1.0
distance: 0.02591463414634143

```

# next steps

1. nextflow integration for working through multiple comparisons.



## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
