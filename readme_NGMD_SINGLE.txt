
=====================
NGMD_SINGLE.PY README
=====================

Description:
------------
NGMD_SINGLE.py is a Python utility script designed to streamline a series of bioinformatics tasks. The primary functionalities include:
- Running HybPiper for target gene enrichment analysis.
- Using get_organelle for organellar genome assembly.
- Performing sequence alignment and variant calling.
- Extracting flanking sequences around variants.

Requirements:
-------------
- Python 3
- Biopython
- HybPiper
- get_organelle

Usage:
------
The script uses argparse to provide a command-line interface and can be invoked with different sub-commands.

1. Running HybPiper for assembly:
   ./NGMD_SINGLE.py run_hybpiper_assemble [hybpiper_assemble_config_file] [target_type]
   - hybpiper_assemble_config_file: Configuration file for HybPiper.
   - target_type: Specify the target type (either 'dna' or 'aa').

2. Retrieving sequences with HybPiper:
   ./NGMD_SINGLE.py run_hybpiper_retrieve [sample_names_file] [target_file]

3. Running get_organelle for organellar genome assembly:
   ./NGMD_SINGLE.py run_getorganelle [getorganelle_config_file] [target_database]

4. Sequence alignment and variant calling:
   ./NGMD_SINGLE.py run_variant_calling [fasta_list]

5. Extracting flanking sequences:
   ./NGMD_SINGLE.py run_extract_flanking_sequences [files_list]

For more detailed instructions and options, refer to the script's internal documentation or use the `-h` flag with the script.

Feedback:
---------
For any feedback or issues, please contact the developer at heyuyang@tju.edu.cn.
