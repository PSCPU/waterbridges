<img src="logo.png" alt="My Logo" width="800" />


This code is written in python3. To run this code, go through the following steps:
1.	Place all three codes (‘1_Identification of water bridges.py’, ‘2_Topological classification of water bridges.py’, and ‘3_Nomenclature of A1-w-N1 waterbridges.py’) into a folder along with ‘.hb2’ files obtained from HBPLUS and respective PDB files.
   
2.	Double click ‘1_Identification of water bridges.py’, the code will run automatically and generate three text files:- 
        ‘waterbridges_ALL.txt’ – contains all identified water bridges.
 	    ‘waterbridges_Base.txt’ – contains all identified water bridges involving at least one ribonucleotide utilizing its nucleobase moiety to form hydrogen bonds with water.
        ‘waterbridges_RiboPhos.txt’ – contains all other identified water bridges involving only phosphate moiety or ribose moiety.

3.	For demonstration, we used ‘waterbridges_Base.txt’. Double click on ‘2_Topological classification of water bridges.py’, it will classify water bridges according to their topology and generate many test files, one for each topology (‘A1-w-N1_Base.txt’, ‘A1-w-B2_Base’)
  
4. Now, run ‘3_Nomenclature of A1-w-N1 waterbridges.py’, to assign a name to each base-mediated water bridge present in ‘waterbridges_Base.txt’. This will generate a file ‘classified_A1-w-B1_waterbridges’.


# Waterbridge: A python code for detection of water bridges in RNA–protein complexes

Water plays an important role in the assembly, stability, and function of RNA–protein complexes. Water forms directional interactions, and its dynamic clusters contribute to molecular recognition. To gain a deeper understanding, we developed a graph theory-based classification scheme for water-mediated amino acid-ribonucleotide motifs, categorizing them into triplets, quartets, and quintet bridging topologies, with further sub-topologies. This categorization not only enhances insights into biomolecular dynamics but also informs the rational design of RNA–protein complexes, providing a framework for potential applications in bioinformatics and therapeutics.

<!--- BADGES: START --->
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
![Static Badge](https://img.shields.io/badge/build-MIT-brightgreen?style=flat&logo=gitbook&logoColor=black&logoSize=auto&label=License&labelColor=%23a9a9a9&color=brightgreen&link=https%3A%2F%2Fgithub.com%2FRamanCompChem%2Fwaterbridges%2Fblob%2Fmain%2FLICENSE&link=https%3A%2F%2Fgithub.com%2FRamanCompChem%2Fwaterbridges%2Fblob%2Fmain%2FLICENSE)
![Static Badge](https://img.shields.io/badge/build-v0.1.3-brightgreen?style=flat&logo=pypi&logoColor=blue&logoSize=auto&label=PyPi&labelColor=%23a9a9a9&color=brightgreen&link=https%3A%2F%2Fpypi.org%2Fproject%2Fwaterbridge%2F&link=https%3A%2F%2Fpypi.org%2Fproject%2Fwaterbridge%2F)
![Static Badge](https://img.shields.io/badge/build-v0.1.3-brightgreen?style=flat&logo=appveyor&logoColor=violet&logoSize=auto&label=Installable%20via%20pip&color=orange&link=https%3A%2F%2Fgithub.com%2FRamanCompChem%2Fwaterbridges%2Fblob%2Fmain%2FLICENSE&link=https%3A%2F%2Fgithub.com%2FRamanCompChem%2Fwaterbridges%2Fblob%2Fmain%2FLICENSE)
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
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Benchmark](#benchmark)
- [Example Files](#example-files)
- [Tests](#tests)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Installation

1. On Windows machine:
      Save the repository to your local machine and python scripts can be run they are.
   
2. On Linux machine
      Directly via pip
      ```bash
      pip install waterbridge
      ```

## How to run

After installation, you can use the `waterbridge` command in your terminal. Here’s how to run the different parts of the tool:

### Command-Line Usage

```bash
waterbridge -1
