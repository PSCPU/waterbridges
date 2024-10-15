<img src="logo.png" alt="My Logo" width="800" />


# Waterbridge: A python code for detection of water bridges in RNA–protein complexes

Water plays an important role in the assembly, stability, and function of RNA–protein complexes. Water forms directional interactions, and its dynamic clusters contribute to molecular recognition. To gain a deeper understanding, we developed a graph theory-based classification scheme for water-mediated amino acid-ribonucleotide motifs, categorizing them into triplets, quartets, and quintet bridging topologies, with further sub-topologies. This categorization not only enhances insights into biomolecular dynamics but also informs the rational design of RNA–protein complexes, providing a framework for potential applications in bioinformatics and therapeutics.

<!--- BADGES: START --->
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
![Static Badge](https://img.shields.io/badge/build-MIT-brightgreen?style=flat&logo=gitbook&logoColor=black&logoSize=auto&label=License&labelColor=%23a9a9a9&color=brightgreen&link=https%3A%2F%2Fgithub.com%2FRamanCompChem%2Fwaterbridges%2Fblob%2Fmain%2FLICENSE&link=https%3A%2F%2Fgithub.com%2FRamanCompChem%2Fwaterbridges%2Fblob%2Fmain%2FLICENSE)
![Static Badge](https://img.shields.io/badge/PyPI-v0.2.1-brightgreen?logo=pypi&logoColor=darkblue&labelColor=lightgrey&link=https%3A%2F%2Fpypi.org%2Fproject%2Fwaterbridge%2F&link=https%3A%2F%2Fpypi.org%2Fproject%2Fwaterbridge%2F)
![Static Badge](https://img.shields.io/badge/Installable%20via%20pip-v0.2.1-orange?logo=pypi&logoColor=blue&link=https%3A%2F%2Fpypi.org%2Fproject%2Fwaterbridge%2F&link=https%3A%2F%2Fpypi.org%2Fproject%2Fwaterbridge%2F)
<!--- BADGES: END --->

## Features

- **Water Bridge Identification:** Identifies water bridges in RNA–protein complexes based on the output of HBPLUS.
- **Topological Classification:** Analyzes and classifies the topology of identified water bridges, such as A1:w:N1, A1:w:N2, and so on.
- **Naming Water Bridges:** Writes nomenclature for the identified topologies.
- **Automatic and Fast:** Automatically processes all .hbplus files in the current directory.
- **Installation:** Simple to install by pip .

## Table of Contents
- [Installation](#installation)
- [How to run](#how-to-run)
- [Dataset](#dataset)
- [Output Files](https://github.com/RamanCompChem/waterbridges/tree/main/Result)
- [License](#license)
- [Contact](#contact)

## Installation
  
1.  On Linux machine:
      Directly via pip
      ```bash
      pip install waterbridge
      ```
2.  On Windows machine:
      Save the repository to your local machine and python scripts can be run as they are.

## How to run
1.  On Linux machine:
   
      Installation through pip will automatically adds 'waterbridge' to your environment. You can simply use the `waterbridge` command in your terminal alnog with arguments in workring directory. For this place HBPLUS outputs and respective PDB files to working directory.
     ### Command-Line Usage
     ```bash
     waterbridge -1 | -2 | -3
     ```
      - `-1`: Identify water bridges in RNA–protein complexes. For this, ensure that the output files from HBPLUS are placed in the current directory.
      - `-2`: Analyze the output of option 1 and identify the topology of the water bridges, such as A1:w:N1, A1:w:N2, and others.
      - `-3`: Provide nomenclature for the identifed topology only.


2.  On Windows machine:
   
   - Place all three codes (‘1_Identification of water bridges.py’, ‘2_Topological classification of water bridges.py’, and ‘3_Nomenclature of waterbridges.py’) into a folder along with ‘.hb2’ files obtained from HBPLUS and respective PDB files.

   - Run three python codes in sequence:
      - 1_Identification of water bridges.py (identify water bridges)
      - 2_Topological classification of water bridges.py (classify water bridges)
      - 3_Nomenclature of waterbridges.py (assign a name to each identified water bridge)


## Dataset

A dataset of RNA–protein complexes was created by extracting all 3D structures containing at least one water molecule and having resolution better than 2.5 Å available at the Protein Data Bank (PDB) before March 15, 2022. To remove redundancy, complexes with over 30% sequence identity were filtered out using the [CDHIT suite](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3516142/), resulting in a final nonredundant dataset of 329 crystal structures. These 329 crystal structures were analyzed using [HBPLUS](https://pubmed.ncbi.nlm.nih.gov/8182748/) to identify hydrogen bonds within RNA–protein complexes, resulting in the generation of ‘.hb2’ files.

The complete dataset of 329 PDB files, along with their corresponding .hb2 files, is available in the [DATASET](https://github.com/RamanCompChem/waterbridges/tree/main/Dataset) directory of this repository.


## Sample output
<img src="sample output.png" alt="Sample" width="900" />


## License
Waterbridge is licensed under the [MIT License](https://opensource.org/license/mit), Version 2.0. You may obtain a copy of the License at [https://github.com/RamanCompChem/waterbridges/blob/main/LICENSE](https://github.com/RamanCompChem/waterbridges/tree/main/Dataset). 


## Contact

For queries or support, please feel free to contact us:

- **Dr. Purshotam Sharma**: [psharmapuchd@gmail.com](mailto:psharmapuchd@gmail.com)
