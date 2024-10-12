############################################################################################################################
#                                       Nomenclature of waterbridges                                                       #
############################################################################################################################
import math
import os

to_update = {':A': ':rA', ':U': ':rU',':G': ':rG',':C': ':rC',}

def group_list_prefix(old_list):
    grouped = {}
    for item in old_list:
        prefix = item[:9]
        if prefix not in grouped:
            grouped[prefix] = []
        grouped[prefix].append(item)
    new_list = list(grouped.values())
    return new_list

# User-defined function to convert PDB record of a water molecule from HBPLUS format (B0026-HOH) 
# to original PDB format (HOH B 26).  
def convert_water_to_pdb(W_id):
    W_id_a = W_id[0]
    W_id_b = W_id[1:5]
    W_id = W_id_a + str(int(W_id_b))
    
    W_id = list(W_id)
    while len(W_id) <5:
        W_id.insert(1, ' ')
    W_id = ''.join(W_id)
    w_pdb = 'O   HOH ' + W_id
    return w_pdb

# User-defined function to convert PDB record of nucleotide to HBPLUS format to original PDB format.  
def convert_nt_to_pdb(idd):
    id_a = idd[0]
    id_b = idd[1:5]
    id_c = idd[8]
    id_atom = idd[10:13]
    chainresi = id_a + str(int(id_b))
    chainresi = list(chainresi)
    while len(chainresi) <5:
        chainresi.insert(1, ' ')
    chainresi.insert(0, ' ' )
    chainresi.insert(0, id_c)
    chainresi.insert(0, id_atom)
    while len(chainresi) <11:
        chainresi.insert(1, ' ')
    nt_pdb = ''.join(chainresi)
    return nt_pdb

# User-defined function to decide interacting edge of nucleotide.
def decision(nt1, nt2, hoh):
    global nt1_xyz, nt2_xyz, hoh_xyz 
    for cords in cord_data:
        if hoh in cords and 'HETATM' in cords: 
            hoh_xyz = [float(cords[31:38].strip()), float(cords[39:46].strip()), float(cords[47:55].strip())]
        if nt1 in cords and 'ATOM' in cords: 
            nt1_xyz = [float(cords[31:38].strip()), float(cords[39:46].strip()), float(cords[47:55].strip())]
        if nt2 in cords and 'ATOM' in cords: 
            nt2_xyz = [float(cords[31:38].strip()), float(cords[39:46].strip()), float(cords[47:55].strip())]
    if (math.dist(hoh_xyz, nt1_xyz)) < (math.dist(hoh_xyz, nt2_xyz)): near_atom = nt1[0:3].strip()
    else: near_atom = nt2[0:3].strip()
    return near_atom
    
# User-defined function to sort ribonucleotides in nomenclature:
def sort_nucleos(coded_ribonucleotides_list):
    WC_list = [nu for nu in coded_ribonucleotides_list if '(WC)' in nu] 
    HG_list = [nu for nu in coded_ribonucleotides_list if '(HG)' in nu]
    SG_list = [nu for nu in coded_ribonucleotides_list if '(SG)' in nu]
    PH_list = [nu for nu in coded_ribonucleotides_list if '(Ph)' in nu]
    RB_list = [nu for nu in coded_ribonucleotides_list if '(Rb)' in nu]
    
    WC_list.sort() 
    HG_list.sort()
    SG_list.sort()
    PH_list.sort()
    RB_list.sort()
    
    final_parts = []
    
    # Add parts in the specified order
    if WC_list: final_parts.append('|'.join(WC_list))
    if HG_list: final_parts.append('|'.join(HG_list))
    if SG_list: final_parts.append('|'.join(SG_list))
    if PH_list: final_parts.append('|'.join(PH_list))
    if RB_list: final_parts.append('|'.join(RB_list))
    
    final_string = '|'.join(final_parts)
    return final_string

# User-defined function to sort amino acids in nomenclature:
def sort_amino_acids(coded_amino_acids_list):
    # part 1: categorize amino acids based on their suffixes
    main_chain_list = []
    main_chain_side_list = []
    side_chain_list = []
    
    for acid in coded_amino_acids_list:
        if '(m)' in acid: main_chain_list.append(acid)
        elif '(ms)' in acid: main_chain_side_list.append(acid)
        elif '(s)' in acid: side_chain_list.append(acid)
            
    # sort each category in reverse order as amino acids will be written on left side of :w:
    main_chain_list.sort(reverse=True)
    main_chain_side_list.sort(reverse=True)
    side_chain_list.sort(reverse=True)
    
    # exclude empty lists
    final_parts = []
    if side_chain_list: final_parts.append('|'.join(side_chain_list))
    if main_chain_side_list: final_parts.append('|'.join(main_chain_side_list))
    if main_chain_list: final_parts.append('|'.join(main_chain_list))
    
    final_string = '|'.join(final_parts)
    return final_string

# List of nucleotides and amino acids.
NT = ['-  A', '-  U', '-  G', '-  C']
AA = ['-GLY', '-ALA', '-SER', '-THR', '-CYS', '-VAL', '-LEU', '-ILE', '-MET', '-PRO', '-PHE', '-TYR', '-TRP', '-ASP', '-GLU', '-ASN', '-GLN', '-HIS', '-LYS', '-ARG']

# Combination of atoms to decide interacting portion of amino acid
MC = ['O  ', 'N  ', 'OXT', ['N  ', 'O  '], ['N  ', 'OXT'], ['O  ', 'OXT'], ['N  ', 'O  ', 'OXT']] # Main chain
mc_atoms = ['O', 'N', 'OXT']
# Combination of atoms to decide interacting portion of amino acid
WC_edge = [['N1', 'N6'],  ['N1', 'N2'], ['N1', 'O6'],  ['N3', 'N4'], ['N3', 'O2'], ['N3', 'O4']] # Watson-Crick edge
HG_edge = [['N6', 'N7'],  ['N7', 'O6'], ['N7', 'O5\'']] # Hoogsteen edge
SG_edge = [['N3', 'N9'], ['N3', 'O2\''], ['N1', 'O2'], ['N2', 'N3'], ['O2', 'O2\''], ['O2\'', 'O3\''], ['N3', 'O4\''], ['N3', "O2\'", "O3\'"]] # Sugar edge
rib1 = [['O4\''], ['O3\'']]
rib2 = [['O3\'', 'O4\''], ['O4\'', 'O5\''],['O2\'', 'O3\'']]
phos1 = [['OP1'], ['OP2'], ['O5\'']]
phos2 = [['OP1', 'OP2'], ['O5\'', 'OP1'], ['O5\'', 'OP2']]

def nucleotide_three_atoms(ntt):
    nt_atoms = []
    for nt in ntt:
        nt_a = nt[-3:].strip()  # interacting atoms
        nt_resi = str(int(nt[1:5])) # residue number
        nt_nt = 'r' + nt[8] # nucleotide
        nt_F = nt_nt + nt_resi
        nt_atoms.append(nt_a)
    nt_atoms.sort()
    if nt_atoms in SG_edge: z_n1 = '(SG)'
    else: z_n1 = 'ERROR'
    final = nt_F + z_n1
    return final

def nucleotide_two_atoms(ntt):
    nt_atoms = []
    for nt in ntt:
        nt_a = nt[-3:].strip()  # interacting atoms
        nt_resi = str(int(nt[1:5])) # residue number
        nt_nt = 'r' + nt[8] # nucleotide
        nt_F = nt_nt + nt_resi
        nt_atoms.append(nt_a)
    nt_atoms.sort()
    if nt_atoms in WC_edge: z_n1 = '(WC)'
    elif nt_atoms in HG_edge: z_n1 = '(HG)'
    elif nt_atoms in SG_edge: z_n1 = '(SG)'
    elif nt_atoms in phos2: z_n1 = '(Ph)'
    elif nt_atoms in rib2: z_n1 = '(Rb)'
    else: z_n1 = 'ERROR'
    final = nt_F + z_n1
    return final

def nucleotide_one_atom(ntt):
    nt = ntt[0] + ' ' # nnt is a list of single entry, so for using that entry as a string 
    nt_a = nt[-3:]  # interacting atoms
    nt_resi = str(int(nt[1:5])) # residue number
    nt_nt = 'r' + nt[8] # nucleotide
    nt_F = nt_nt + nt_resi
    w_pdb =  convert_water_to_pdb(x)    
    if '-  A' in nt:
        if ' N6 ' in nt:
            new_n_h = nt.replace('N6', 'C5')
            new_n_w = nt.replace('N6', 'N1')
            nt_pdb_h = convert_nt_to_pdb(new_n_h)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_h, nt_pdb_w, w_pdb)
            if nearest_atom == 'C5': z_n1 = '(HG)'
            else: z_n1 = '(WC)'
        elif ' N7 ' in nt: z_n1 = '(HG)'
        elif ' N1 ' in nt: z_n1 = '(WC)'
        elif ' N3 ' in nt: z_n1 = '(SG)'
        elif ' N9 ' in nt: z_n1 = '(SG)'
        elif ' O2\'' in nt: z_n1 = '(SG)'
        elif ' O4\'' in nt: z_n1 = '(Rb)'
        elif ' O3\'' in nt: z_n1 = '(Rb)'
        elif ' OP1' in nt: z_n1 = '(Ph)'
        elif ' OP2' in nt: z_n1 = '(Ph)'
        elif ' O5\'' in nt: z_n1 = '(Ph)'
        else: z_n1 = 'ERR'

    elif '-  G' in nt:
        if ' O6 ' in nt:
            new_n_h = nt.replace('O6', 'C5')
            new_n_w = nt.replace('O6', 'N1')
            nt_pdb_h = convert_nt_to_pdb(new_n_h)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_h, nt_pdb_w, w_pdb) #nt_pdb_h will decide edge
            if nearest_atom == 'C5': z_n1 = '(HG)'
            else: z_n1 = '(WC)'
        elif ' N2 ' in nt:
            new_n_s = nt.replace('N2', 'N3')
            new_n_w = nt.replace('N2', 'N1')
            nt_pdb_s = convert_nt_to_pdb(new_n_s)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_s, nt_pdb_w, w_pdb) #nt_pdb_s will decide edge
            if nearest_atom == 'N3': z_n1 = '(SG)'
            else: z_n1 = '(WC)'
        elif ' N7 ' in nt: z_n1 = '(HG)'
        elif ' N1 ' in nt: z_n1 = '(WC)'
        elif ' N3 ' in nt: z_n1 = '(SG)'
        elif ' N9 ' in nt: z_n1 = '(SG)'
        elif ' O2\'' in nt: z_n1 = '(SG)'
        elif ' O4\'' in nt: z_n1 = '(Rb)'
        elif ' O3\'' in nt: z_n1 = '(Rb)'
        elif ' OP1' in nt: z_n1 = '(Ph)'
        elif ' OP2' in nt: z_n1 = '(Ph)'
        elif ' O5\'' in nt: z_n1 = '(Ph)'
        else: z_n1 = 'ERR'

    elif '-  U' in nt:
        if ' O4 ' in nt:
            new_n_h = nt.replace('O4', 'C5')
            new_n_w = nt.replace('O4', 'N3')
            nt_pdb_h = convert_nt_to_pdb(new_n_h)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_h, nt_pdb_w, w_pdb) #nt_pdb_s will decide edge
            if nearest_atom == 'C5': z_n1 = '(HG)'
            else: z_n1 = '(WC)' 
        elif ' O2 ' in nt:
            new_n_s = nt.replace('O2', 'N1')
            new_n_w = nt.replace('O2', 'N3')
            nt_pdb_s = convert_nt_to_pdb(new_n_s)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_s, nt_pdb_w, w_pdb) #nt_pdb_s will decide edge
            if nearest_atom == 'N1': z_n1 = '(SG)'
            else: z_n1 = '(WC)' 
        elif ' N1 ' in nt: z_n1 = '(SG)'
        elif ' N3 ' in nt: z_n1 = '(WC)'
        elif ' O2\'' in nt: z_n1 = '(SG)'
        elif ' O4\'' in nt: z_n1 = '(Rb)'
        elif ' O3\'' in nt: z_n1 = '(Rb)'
        elif ' OP1' in nt: z_n1 = '(Ph)'
        elif ' OP2' in nt: z_n1 = '(Ph)'
        elif ' O5\'' in nt: z_n1 = '(Ph)'
        else: z_n1 = 'ERR'
            
    elif '-  C' in nt:
        if ' O2 ' in nt:
            new_n_s = nt.replace('O2', 'N1')
            new_n_w = nt.replace('O2', 'N3')
            nt_pdb_s = convert_nt_to_pdb(new_n_s)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_s, nt_pdb_w, w_pdb) #nt_pdb_s will decide edge
            if nearest_atom == 'N1': z_n1 = '(SG)'
            else: z_n1 = '(WC)' 
        elif ' N4 ' in nt:
            new_n_h = nt.replace('N4', 'C5')
            new_n_w = nt.replace('N4', 'N3')
            nt_pdb_h = convert_nt_to_pdb(new_n_h)
            nt_pdb_w = convert_nt_to_pdb(new_n_w)
            nearest_atom = decision(nt_pdb_h, nt_pdb_w, w_pdb) #nt_pdb_s will decide edge
            if nearest_atom == 'C5': z_n1 = '(HG)'
            else: z_n1 = '(WC)' 
        elif ' O2\' ' in nt: z_n1 = '(SG)'
        elif ' N1 ' in nt: z_n1 = '(SG)'
        elif ' N3 ' in nt: z_n1 = '(WC)'
        elif ' O4\'' in nt: z_n1 = '(Rb)'
        elif ' O3\'' in nt: z_n1 = '(Rb)'
        elif ' OP1' in nt: z_n1 = '(Ph)'
        elif ' OP2' in nt: z_n1 = '(Ph)'
        elif ' O5\'' in nt: z_n1 = '(Ph)'
        else: z_n1 = 'ERR'
    final = nt_F + z_n1
    return final
        


files_from_previous = []
path = os.getcwd()
directory = os.listdir(path)  
for fis in directory: 
    if '.txt' in fis and '-w-' in fis and 'waterbridges' not in fis :
        files_from_previous.append(fis)
files_from_previous.sort()

# Read previously generated file containing water bridges.
for previousfile in files_from_previous:
    file = open(previousfile, 'r')
    outfile = 'named_'+ previousfile[:2] + 'w' + previousfile[5:7] + '_waterbridges.txt'
    data = file.readlines()
    for lines in data:
        print(lines)
        if 'pdb' in lines: 
            #print(lines)
            # Reading PDB files to extract coordinates of atoms.
            cord_file = open(lines.split('.')[0].split(' ')[2] + '.ent')
            cord_data = cord_file.readlines()
            print(lines.strip(), file=open(outfile, "a"))
        if 'pdb' not in lines and 'HOH' in lines:
            AA_info = [] # List of amino acids.
            NT_info = [] # List of nucleotides.
            # Determine whether the water bridge is cyclic or acyclic.
            if "NIL" in lines: zc = ''
            else: zc = 'cyc-'
            # Extracting amino acids and nucleotides.
            x = lines.split('[')[0]
            for i in range(16, len(x), 23):
                y = lines[i:i+14].strip()
                y1 = y[5:9]
                if y1 in NT: 
                    NT_info.append(y) 
                elif y1 in AA: 
                    AA_info.append(y)
            AA_info = group_list_prefix(AA_info)
            NT_info = group_list_prefix(NT_info)
            #print(AA_info, NT_info)
            
            code_nt_ALL = [] # contains all coded nucleotides
            for ntt in NT_info:
                # Case 1: When 3 atoms of a nucleotide are forming hydrogen bonds with water.
                if len(ntt) == 3:
                    z_n = nucleotide_three_atoms(ntt)
                    #print(z_n)
                # Case 2: When 2 atoms of a nucleotide are forming hydrogen bonds with water.
                elif len(ntt) == 2:
                    z_n = nucleotide_two_atoms(ntt)
                    #print(z_n)
                # Case 3: When only 1 atom of nucleotide is forming a hydrogen bond with water.
                elif len(ntt) == 1:
                    z_n = nucleotide_one_atom(ntt)
                    #print(z_n)
                code_nt_ALL.append(z_n)
        
            #print(code_nt_ALL)
            right_nt = sort_nucleos(code_nt_ALL)
            #print(right_nt)
            
            # Deciding interacting portion of amino acid.
            code_aa_ALL = [] # contains all coded nucleotides
            for aat in AA_info:
                aa_id_ms = []
                #print(aat)
                for aa in aat:
                    aa_resi = str(int(aa[1:5])) # residue number
                    aa_aa = aa [6:9] # amino acid
                    aa_atom = aa[10:].strip()  # interacting amino acid atom
                    if aa_atom in mc_atoms:
                        z_a = 'm'
                    else: z_a = 's'
                    if z_a not in aa_id_ms: aa_id_ms.append(z_a) 
                aa_id_ms.sort()
                joined_mc_sc = ''.join(aa_id_ms)
                code_aa = aa_aa + aa_resi + '(' + joined_mc_sc + ')'
                code_aa_ALL.append(code_aa)
            #print(code_aa_ALL)
            left_aa = sort_amino_acids(code_aa_ALL)
            #print(left_aa)
            
            name =  zc + left_aa + ':w:' + right_nt
            interaction = lines.strip() + ' —> ' + name
            print(interaction)
            print(interaction.strip(), file=open(outfile, "a"))
    file.close()


                


            
        
                
                 
            
            
            
            
            

