
====================
NGMD_ALL.PY README
====================

Description:
------------
NGMD_all.py is a comprehensive Python utility script designed to perform a series of bioinformatics tasks in an integrated manner. The main functionalities include:
- Running HybPiper for target gene enrichment analysis, followed by sequence retrieval.
- Using get_organelle for organellar genome assembly.
- Preparing sequences for alignment and variant calling.

Requirements:
-------------
- Python 3
- Biopython
- HybPiper
- get_organelle

Usage:
------
The script offers a command-line interface with sub-commands for different pipelines:

1. Running HybPiper for assembly and retrieval:
   ```bash
   ./NGMD_all.py run_hybpiper_assemble_and_retrieve [hybpiper_config_file] [output_dir]
   ```
   - hybpiper_config_file: Configuration file for HybPiper.
   - output_dir: Directory for saving results.

2. Running get_organelle for organellar genome assembly:
   ```bash
   ./NGMD_all.py run_getorganelle [getorganelle_config_file] [target_database] [output_dir]
   ```

For more detailed instructions and options, refer to the script's internal documentation or use the `-h` flag with the script.

Feedback:
---------
For any feedback or issues, please contact the developer at heyuyang@tju.edu.cn.
