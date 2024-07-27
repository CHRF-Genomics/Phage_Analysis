# Phage Specific Protein Sequence Extractor 

This script extracts the translation sequences of the tail fiber protein gene and the terminase large subunit gene from GBK files. It then joins these sequences and writes them to an output file.

## Prerequisites

- Python 3.x
- Biopython
- An input directory containing GBK files

## Installation

1. Install Biopython if you haven't already:

```bash
pip install biopython

```

## Usage
Place your GenBank (.gbk) files in a directory.
Update the input_dir_gbk variable in the script to point to this directory.
Run the script.

### Libraries
```python

from Bio import SeqIO
from Bio.Seq import translate
import os

```

### Define Input Directory and List GBK Files
```python

input_dir_gbk = "/path/to/your/gbk_files"
gbk_files = [f for f in os.listdir(input_dir_gbk) if f.endswith(".gbk")]
print(gbk_files)

```
### Extract Translation Sequence of Tail Fiber Protein Gene
```python

def extract_translation_seq_of_tail_fiber_gene(input_dir_gbk, gbk_files):
    gbk_dic_tail = {}
    for gbk_file in gbk_files:
        file_path = os.path.join(input_dir_gbk, gbk_file)
        file_name = os.path.basename(file_path)
        phage_name = os.path.splitext(file_name)[0].split('_pharokka')[0]

        for rec in SeqIO.parse(file_path, "gb"):
            for feature in rec.features:
                for key, val in feature.qualifiers.items():
                    if "tail fiber protein" in val:
                        translation_str = str(feature.qualifiers['translation'])
                        gbk_dic_tail[phage_name] = translation_str

    gbk_dic_tail_keys_sorted = sorted(list(gbk_dic_tail.keys()))
    gbk_dic_tail_keys_sorted_dict = {i: gbk_dic_tail[i] for i in gbk_dic_tail_keys_sorted}

    return gbk_dic_tail_keys_sorted_dict

```
### Extract Translation Sequence of Terminase Large Subunit Gene
```python

def extract_translation_seq_of_terminase_gene(input_dir_gbk, gbk_files):
    gbk_dic_terminase = {}
    for gbk_file in gbk_files:
        file_path = os.path.join(input_dir_gbk, gbk_file)
        file_name = os.path.basename(file_path)
        phage_name = os.path.splitext(file_name)[0].split('_pharokka')[0]

        for rec in SeqIO.parse(file_path, "gb"):
            for feature in rec.features:
                for key, val in feature.qualifiers.items():
                    if "terminase large subunit" in val:
                        translation_str = str(feature.qualifiers['translation'])
                        gbk_dic_terminase[phage_name] = translation_str

    gbk_dic_terminase_keys_sorted = sorted(list(gbk_dic_terminase.keys()))
    gbk_dic_terminase_keys_sorted_dict = {i: gbk_dic_terminase[i] for i in gbk_dic_terminase_keys_sorted}

    return gbk_dic_terminase_keys_sorted_dict

```

### Example Usage
```python
extract_translation_seq_of_terminase_gene(input_dir_gbk, gbk_files)
extract_translation_seq_of_tail_fiber_gene(input_dir_gbk, gbk_files)


```
This will create two separate files that contain the sequences of the tail fiber protein and 
terminase large subunit genes for each phage.




