# PYTHON COURSE - WORK IN PROGRESS

# AED calculator tool

Gene models can be assessed using mapped transcript or protein evidence using a metric called annotation edit distance (AED). The theory behind this can be found at https://link.springer.com/article/10.1186/1471-2105-10-67. 

This can be a useful measure to assess gene predictor accuracy. Here you can find a simple tool for calculating AED per gene model, given some evidence.

# Installation

```bash
git clone https://github.com/hlmwhite/AED_calculator_tool.git
```

```

## Usage

For single use cases (i.e. one gene model one evidence model), simply run:

```bash

python3 example_calculator_v2.py -A genome.gff -E evidence.gff -a "<annotation gene>" -e "<evidence model to assess AED>"

```

An example using the test files inlcuded here is:

```bash

python3 example_calculator_v2.py -A test.genome.gff -E stringtie.test.gff -a "rna-XM_034815643.1" -e "STRG.95.1"

```

# next steps

1. nextflow integration for working through multiple comparisons.



## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
