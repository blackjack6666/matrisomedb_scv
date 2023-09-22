"""
apply [PTM] to multiple peptide lines, in a {},separated by line break. e.g {PEPTIDE1\nPEPTIDE2\nPEPTIDE3\n}[PTM1]

"""

import os
import re
import traceback

import pymol
import numpy as np
from collections import defaultdict
import sys
from pymol2glmol import *
from glob import glob
import ahocorasick
import commons
import sys
import json
import time
import sqlite3
from gen_3d import show_cov_3d_v2
from bs4 import BeautifulSoup as Soup


# aa_str = ''.join([aa for aa in commons.aa_mass_table])

def automaton_trie(peptide_list):
    A = ahocorasick.Automaton()
    for idx, peptide in enumerate(peptide_list):
        A.add_word(peptide, (idx, peptide))
    A.make_automaton()
    return A


def automaton_matching(A, seq_line):
    result = []
    for end_idx, (insert_order, original_value) in A.iter(seq_line):
        start_idx = end_idx - len(original_value) + 1
        result.append((start_idx, end_idx, original_value))
        assert seq_line[start_idx:start_idx + len(original_value)] == original_value
    return result


def fasta_reader(fasta_file_path):
    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def fasta_reader2(fasta_file_path):
    # use gene_uniprot as key
    gene_protein_seq_dict = {}
    with open(fasta_file_path, 'r') as f_o:
        file_split = f_o.read().split('\n>')

    for each in file_split:
        first_line, seq = each.split('\n')[0], ''.join(each.split('\n')[1:])
        uniprot_id = first_line.split('|')[1]
        gene = first_line.split('GN=')[1].split(' ')[0] if 'GN=' in first_line else 'N/A'
        gene_protein_seq_dict[gene + '_' + uniprot_id] = seq
    return gene_protein_seq_dict


def fasta_reader3(fasta_path: str):
    protein_dict = {}
    with open(fasta_path, 'r') as f_o:
        file_split = f_o.read().split('\n>')

    for each in file_split:
        first_line, seq = each.split('\n')[0], ''.join(each.split('\n')[1:])
        uniprot_id = first_line.split('|')[1]
        gene = first_line.split('GN=')[1].split(' ')[0] if 'GN=' in first_line else 'N/A'
        des = ' '.join(first_line.split(' ')[1:]).split(' OS=')[0]
        protein_dict[uniprot_id] = (seq, gene, des)
    return protein_dict


def peptide_counting(peptide_tsv_file):
    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list


def my_replace(match_obj):
    match_obj = match_obj.group()
    matched_aa = match_obj[0]
    if matched_aa != 'n':
        return matched_aa  # gives back the first element of matched object as string
    else:
        # if first match is n, then n acetylation, get rid of n
        return ''

def pdb_file_reader(pdb_file_list: list):
    """
    reads a pdb file into protein sequence
    :param pdb_file:
    :return:
    """
    aa_dict = {'ALA': 'A',
               'ARG': 'R',
               'ASN': 'N',
               'ASP': 'D',
               'ASX': 'B',
               'CYS': 'C',
               'GLU': 'E',
               'GLN': 'Q',
               'GLX': 'Z',
               'GLY': 'G',
               'HIS': 'H',
               'ILE': 'I',
               'LEU': 'L',
               'LYS': 'K',
               'MET': 'M',
               'PHE': 'F',
               'PRO': 'P',
               'SER': 'S',
               'THR': 'T',
               'TRP': 'W',
               'TYR': 'Y',
               'VAL': 'V'}

    aa_reg_str = '|'.join([key for key in aa_dict])

    import re
    import os

    pdb_protein_seq_dict = {}

    for pdb_file in pdb_file_list:
        with open(pdb_file, 'r') as f_o:
            f_split = f_o.read().split('\nATOM')[1:]
            pos_aa_list = [(int(re.search('\d+(?=\s+[+-]?\d+\.)', each).group()),
                            re.search(aa_reg_str, each).group(0)) for each in f_split]
            protein_seq = ''
            for i in range(len(pos_aa_list) - 1):
                if pos_aa_list[i + 1][0] == pos_aa_list[i][0]:
                    continue
                else:
                    protein_seq += aa_dict[pos_aa_list[i][1]]

        # add last aa
        protein_seq += aa_dict[pos_aa_list[-1][1]]
        pdb_protein_seq_dict[os.path.split(pdb_file)[-1]] = (protein_seq, 'User upload')
        print(protein_seq)
    return pdb_protein_seq_dict


def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, pdbs, regex_dict=None, psm_group_dict = {}):
    """
    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """
    from collections import Counter
    from collections import defaultdict
    from heapq import heappush, heappop

    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    h = []
    regex_pat = '\w{1}\[\d+\.?\d+\]'  # universal ptm pattern
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)
    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line_dict = {each: commons.zero_line_for_seq(seq_line)
                           for each in regex_dict} if regex_dict else False

    separtor_pos_array = commons.separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1] + 1] += len(peptide_psm_dict[matched_pep])
        # ptm assign might need to be optimized
        if ptm_index_line_dict:  # if ptm enabled
            # print(regex_dict)

            for psm in peptide_psm_dict[matched_pep]:
                # first label group PTM
                if psm in psm_group_dict:
                    group_ptm = psm_group_dict[psm]
                    ptm_index_line_dict[group_ptm][tp[0]:tp[1] + 1] +=1

                # then label each PTM inside each PSM
                for ptm in ptm_index_line_dict:
                    print (ptm)
                    new_psm = re.sub('n?\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')).replace('*',
                                                                                                             ptm.replace(
                                                                                                                 '\\',
                                                                                                                 ''))
                    ptm_mod = re.findall(ptm, new_psm)

                    if ptm_mod:

                        for ele in ptm_mod:
                            num_of_mod = len(
                                re.findall(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm))
                            PTM_index = [m.start() for m in
                                         re.finditer(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'),
                                                     new_psm)]
                            PTM_index_clean = [ind - num * (len(ele) - 1) for ind, num in
                                               zip(PTM_index, range(num_of_mod))]
                            for indx in PTM_index_clean:
                                ptm_index_line_dict[ptm][tp[0] + indx] += 1


    time_start = time.time()
    for i in range(len(separtor_pos_array) - 1):
        zero_line_slice = zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
        if len(zero_line_slice) > 0:
            percentage_cov = np.count_nonzero(zero_line_slice) / len(zero_line_slice) * 100
        else:
            percentage_cov = 0.0
        if percentage_cov != 0.0:
            heappush(h, (percentage_cov, (id_list[i], protein_dict[id_list[i]]), zero_line_slice.tolist(),
                         any(id_list[i] in pdb for pdb in pdbs)))
        id_freq_array_dict[id_list[i]] = zero_line_slice.tolist()
        if ptm_index_line_dict:
            id_ptm_idx_dict[id_list[i]] = {ptm: np.array(
                np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[0]).tolist()
                                           for ptm in ptm_index_line_dict}

    print(time.time() - time_start)

    return id_freq_array_dict, id_ptm_idx_dict, [heappop(h) for i in range(len(h))][::-1]


def pq_gen(psm_list, protein_dict, pdbs):
    from heapq import heappush, heappop
    info_dict= {}
    h = []
    regex_pat = '\w{1}\[\d+\.?\d+\]'  # universal ptm pattern
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)
    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)

    separtor_pos_array = commons.separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1] + 1] += len(peptide_psm_dict[matched_pep])

    time_start = time.time()
    heap_id_list = []
    for i in range(len(separtor_pos_array) - 1):
        zero_line_slice = zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
        if len(zero_line_slice) > 0:
            percentage_cov = np.count_nonzero(zero_line_slice) / len(zero_line_slice) * 100
        else:
            percentage_cov = 0.0
        if percentage_cov != 0.0:
            # heappush(h, (percentage_cov, (id_list[i], protein_dict[id_list[i]]), zero_line_slice.tolist(),
            #              any(id_list[i] in pdb for pdb in pdbs)))
            # heap_id_list.append(id_list[i])

            info_dict[id_list[i]] = [(percentage_cov, (id_list[i], protein_dict[id_list[i]]), zero_line_slice.tolist(),
                         any(id_list[i] in pdb for pdb in pdbs))]

    # pq_list = [heappop(h) for i in range(len(h))]

    return info_dict


def modified_peptide_from_psm(psm_path):
    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]', line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list


def fetch_all_proteins_from_db(cur):
    cur.execute("SELECT protein FROM pdbstr")
    records = cur.fetchall()
    return records


def create_table(db_name):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    cur.execute('''CREATE TABLE IF NOT EXISTS results (
                    job_number text PRIMARY KEY,
                    pq text NOT NULL,
                    id_ptm_idx_dict text NOT NULL,
                    regex_dict text NOT NULL,
                    background_color text NOT NULL,
                    pdb_dest text NOT NULL) ''')
    con.commit()


def insert_to_table(db_name, job_number, pq, id_ptm_idx_dict, regex_dict, background_color, pdb_dest):
    con = sqlite3.connect(db_name)
    cur = con.cursor()
    cur.execute("INSERT INTO results VALUES (?, ?, ?, ?, ?, ?)",
                (job_number,
                 json.dumps(pq),
                 json.dumps(id_ptm_idx_dict),
                 json.dumps(regex_dict),
                 background_color,
                 pdb_dest)
                )
    con.commit()


def get_ptms(psms):
    ptms = set()
    res = {}
    for p in psms:
        reg = re.search('\w{1}\[\d+\.?\d+\]', p)
        if reg is not None:
            ptms.add(reg.group(0))
    for ptm in ptms:
        res[ptm] = [0, 255, 0]
    return res


def get_ptms_2(psm_str:str):

    # check if {} in string, if in, extract what's the [ptm] after each {}
    if '}' in psm_str and '{' in psm_str:
        psm_split = {each.split('}')[1].split('\n')[0]:each.split('}')[0].split('\n') for each in psm_str.split('{')[1:]}
        reg_dict_update = {group:[0,255,0] for group in psm_split} # assign groups with RGB
        psms = re.sub('\{|\}\[\w+\]','',psm_str).split('\n')
        psm_group_dict = {psm: group for group in psm_split for psm in psm_split[group]}
    else:
        reg_dict_update = {}
        psms = psm_str.split('\n')
        psm_group_dict = {}
    # print (psm_split)
    ptms = set()
    res = {}
    for p in psms:
        reg = re.search('\w{1}\[\d+\.?\d+\]', p)
        if reg is not None:
            ptms.add(reg.group(0))
    for ptm in ptms:
        res[ptm] = [0, 255, 0]
    res.update(reg_dict_update)
    return res, psms, psm_group_dict


if __name__ == '__main__':

    # PSM str with {}[PTM] as a group
    PSMs_str = '{MAGLM[100]TIVTSLLF\nAHHIIP[200]TGSV}[group1]\n{VSKRIPENRVVSYQ[300]LSSRSTCL\nTTKKGQQFCGDPKQEWVQRYMKNLDAKQ}[group2]\nRARAVAVK[400]GPVQRY'

    # get regex dict, PSM list, psm_group dict
    regex_dict,PSMs, psm_group_dict = get_ptms_2(PSMs_str)
    psm_group_dict = {psm:psm_group_dict[psm].replace('[','\[').replace(']','\]').replace('.','\.') for psm in psm_group_dict}
    print (psm_group_dict)
    print (PSMs)
    protein_dict = fasta_reader3('F:/matrisomedb2.0/mat.fasta')
    pdbs = ['xxx']

    # regex_dict = get_ptms(PSMs)
    new_regex_set = {each.replace('[', '\[').replace(']', '\]').replace('.', '\.') for each in regex_dict}
    regex_dict = {val:[0, int(255/(ind+1)), int(128/(ind+1))] for ind,val in enumerate(new_regex_set)}


    id_freq_array_dict, id_ptm_idx_dict, pq = freq_ptm_index_gen_batch_v2(PSMs, protein_dict, pdbs,
                                                                        regex_dict=regex_dict,psm_group_dict=psm_group_dict)

    print(regex_dict)
    print(id_ptm_idx_dict['O00175'])

    # generate 3d cov
    pdb_file = 'F:/matrisomedb2.0/mat_alphafold_pdbs/AF-' + 'O00175' + '-F1-model_v1.pdb'
    output_path = 'multi_psm_test.html'


    dump_dict = show_cov_3d_v2('O00175', 'AAA', pdb_file, np.array(id_freq_array_dict['O00175']).astype(np.int32),
                               id_ptm_idx_dict=id_ptm_idx_dict['O00175'],
                               regex_color_dict=regex_dict,
                               background_color=str(16777215))

    percentage_cov = np.count_nonzero(id_freq_array_dict['O00175'])/len(id_freq_array_dict['O00175'])*100
    print ('cov:',percentage_cov)
    pq_dump = [(percentage_cov, ('O00175', protein_dict['O00175']), id_freq_array_dict['O00175'],
                         True)]
    view_file = open('./view.html')
    soup = Soup(view_file.read(), 'lxml')
    js_data = soup.find(id='data')
    js_data.string = 'let id_ptm_idx_dict = ' + json.dumps({'O00175':id_ptm_idx_dict['O00175']}) + ';' + \
                     'let regex_dict = ' + json.dumps(regex_dict) + ';' + \
                     'let background_color = ' + str(16777215) + ';' + \
                     'let pq = ' + json.dumps(pq_dump) + ';'
    pdb_str = soup.find(id='glmol01_src')
    pdb_str.string = dump_dict['pbdstr']
    ret = soup.find(id='glmol01_rep')
    ret.string = dump_dict['ret']
    output_file = open(output_path, "w")
    # output_file = open('P11276_test.html','w')
    output_file.write(str(soup))
    output_file.close()


