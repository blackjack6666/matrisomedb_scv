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


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string, regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map single protein seq
    :param peptide_list:
    :param protein_seq_string:
    :param regex_pat:
    :return:
    """

    freq_array = np.zeros(len(protein_seq_string))
    PTM_sites_counting = defaultdict(int)
    PTM_loc_list = []

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation

    for pep, new_pep, PTM in zip(peptide_list, new_pep_list, PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:

            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) - 1
            # print (start_pos,end_pos,new_pep)
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:
                    PTM_index = pep.find(ele)
                    # PTM_site = pep[PTM_index] # single amino acid
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos + PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)

    return freq_array, PTM_loc_list, PTM_sites_counting


def freq_ptm_index_gen_batch(psm_list, protein_dict, regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map large psm list on whole proteome
    :param psm_list: psm list with ptms
    :param protein_dict: protein sequence dictionary
    :param regex_pat: regex pattern
    :return:
    """
    from collections import Counter
    from collections import defaultdict
    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}

    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line = commons.zero_line_for_seq(seq_line)
    separtor_pos_array = commons.separator_pos(seq_line)

    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)

    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1] + 1] += len(peptide_psm_dict[matched_pep])
        # print (tp[0],tp[1],matched_pep, 'aho')
        for psm in peptide_psm_dict[matched_pep]:
            psm_mod = re.findall(regex_pat, psm)
            if psm_mod:  # if psm has mod
                for ele in psm_mod:
                    ptm_idx = psm.find(ele)
                    ptm_index_line[tp[0] + ptm_idx] += 1

    for i in range(len(separtor_pos_array) - 1):
        id_freq_array_dict[id_list[i]] = zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]].tolist()
        id_ptm_idx_dict[id_list[i]] = np.nonzero(ptm_index_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[
                                          0] + 1

    return id_freq_array_dict, id_ptm_idx_dict


# def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, regex_dict=None):
#     """
#
#     :param psm_list:
#     :param protein_dict:
#     :param regex_dict: {regex:HEX color}
#     :return:
#     """
#
#     from collections import Counter
#     from collections import defaultdict
#     id_freq_array_dict = {}
#     ptm_site_counting = defaultdict(int)
#     id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
#
#     regex_pat = '\w{1}\[\d+\.?\d+\]' # universal ptm pattern
#     peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
#     for each in psm_list:
#
#         each_reg_sub = re.sub(regex_pat, my_replace, each)
#         peptide_psm_dict[each_reg_sub].append(each)
#
#     # aho mapping
#     id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
#     seq_line = commons.creat_total_seq_line(seq_list, sep="|")
#     zero_line = commons.zero_line_for_seq(seq_line)
#     ptm_index_line_dict = {each:commons.zero_line_for_seq(seq_line) for each in regex_dict} if regex_dict else False
#     separtor_pos_array = commons.separator_pos(seq_line)
#
#     aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
#     for tp in aho_result:
#         matched_pep = tp[2]  # without ptm site
#         zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
#         # ptm assign might need to be optimized
#         if ptm_index_line_dict:  # if ptm enabled
#             for psm in peptide_psm_dict[matched_pep]:
#                 for ptm in regex_dict:
#                     # only keep one ptm in psm if there are multiple for correct index finding
#                     new_psm = re.sub('(?!'+ptm[1:]+')\[\d+\.?\d+\]','',psm)
#
#                     ptm_mod = re.findall(ptm, new_psm)
#                     print(ptm_mod)
#                     if ptm_mod:
#                         for ele in ptm_mod:
#                             print (new_psm)
#                             ptm_idx = new_psm.find(ele)
#                             print(matched_pep, tp[0], ele, ptm_idx)
#                             ptm_index_line_dict[ptm][tp[0] + ptm_idx] += 1
#
#     for i in range(len(separtor_pos_array)-1):
#         zl = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]].tolist()
#         # id_freq_array_dict[id_list[i]] = []
#         for elem in zl:
#             if elem != 0:
#                 id_freq_array_dict[id_list[i]] = zl
#                 break
#
#         # id_freq_array_dict.pop(id_list[i])
#         # id_freq_array_dict[id_list[i]] = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]].tolist()
#         if ptm_index_line_dict:
#             id_ptm_idx_dict[id_list[i]]= {ptm:np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1
#                                           for ptm in ptm_index_line_dict}
#     print (id_ptm_idx_dict)
#     return id_freq_array_dict, id_ptm_idx_dict

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


def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, pdbs, regex_dict=None):
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
                # print (psm)
                for ptm in ptm_index_line_dict:

                    #                     print(ptm)
                    # only keep one ptm in psm if there are multiple for correct index finding
                    #                    new_psm = re.sub('(?!'+ptm[1:]+')\[\d+\.?\d+\]','',psm)
                    #                    new_psm = re.sub('['+aa_str.replace(ptm[0],'')+']\[\d+\.?\d+\]',my_replace,psm)
                    new_psm = re.sub('n?\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')).replace('*',
                                                                                                           ptm.replace(
                                                                                                               '\\',
                                                                                                               ''))
                    ptm_mod = re.findall(ptm, new_psm)
                    # print (ptm_mod)
                    #                     print(ptm_mod)
                    if ptm_mod:
                        #                        for ele in ptm_mod:
                        #                            print (new_psm)
                        #                            ptm_idx = new_psm.find(ele)
                        #                            print(matched_pep, tp[0], ele, ptm_idx)
                        #                            ptm_index_line_dict[ptm][tp[0] + ptm_idx] += 1
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
            id_ptm_idx_dict[id_list[i]] = {ptm.replace('\[','[').replace('\]',']').replace('\.','.'): np.array(
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


if __name__ == '__main__':
    from glob import glob
    import pickle
    # global_protein_psm_dict = json.load(open('F:/matrisomedb2.0/global_protein_psm.dict_fromsample.json', 'r'))
    # sample_prot_psm_dict = json.load(
    #     open('F:/matrisomedb2.0/sample_protein_psm_dict_3.json', 'r'))  # sample name has illgal charcters
    # count =0
    # print (len(sample_prot_psm_dict))
    # for sample in sample_prot_psm_dict:
    #     for char in '\?%*:|"<>/':
    #         if char in sample:
    #             print("Illegal character %s found in %s" % (char, sample))


    # protein_dict = fasta_reader3(sys.argv[3])
    protein_dict = fasta_reader3('F:/matrisomedb2.0/mat.fasta')

    # peptide_file = open(sys.argv[1], "r")
    # peptide_file = open('./peptide_files/P11276_1080D.txt')
    # protein = sys.argv[1].split('/')[-1].split('_')[0]
    # protein = 'P11276'
    # for f in os.listdir(sys.argv[2]):
    #     if f.split('-')[1] == protein:
    #         pdbs.append(os.path.join(sys.argv[2], f))
    # pdbs = ['./pdbs/AF-P11276-F1-model_v1.pdb']
    # print(pdbs)
    # protein_dict = pdb_file_reader(pdbs)

    # print(protein_dict)
    # protein_dict = fasta_reader3('./fastas/uniprot-proteome_UP000000589_sp_only_mouse.fasta')

    # PSMs = peptide_file.read().splitlines()

    # print(protein)
    # print(PSMs)
    # PSMs = pickle.load(open(r'F:\matrisomedb2.0/all_psm.p','rb'))
    # print ('PSMs loaded' )
    # regex_dict = get_ptms(PSMs)
    regex_set = {'M\\[15\\.9949\\]','P\\[15\\.9949\\]','K\\[15\\.9949\\]','n\\[42\\.0106\\]',
              'R\\[0\\.9840\\]','T\\[79\\.9663\\]','S\\[79\\.9663\\]','Y\\[79\\.9663\\]'}
    regex_dict = {each:[0, 255, 0] for each in regex_set}
    # new_regex_set = {each.replace('[', '\[').replace(']', '\]').replace('.', '\.') for each in regex_dict}

    # id_freq_array_dict, id_ptm_idx_dict, pq = freq_ptm_index_gen_batch_v2(PSMs, protein_dict, pdbs,
    #                                                                       regex_dict=new_regex_set)
    # base_path = 'F:/matrisomedb2.0/3dcov_htmls/'
    # pq_dict = pq_gen(PSMs,protein_dict,pdbs)
    # count = 0
    # for each in pq_dict:
    #     count+=1
    #     if count<10:
    #         print (each, pq_dict[each])
    # pickle.dump(pq_dict,open('F:/matrisomedb2.0/pq.data','wb'),protocol=5)

    # load data
    pq_dict = pickle.load(open('F:/matrisomedb2.0/pq.data', 'rb'))
    id_freq_array_dict = pickle.load(open('F:/matrisomedb2.0/glob_prot_freq_dict.p','rb'))
    id_ptm_idx_dict = pickle.load(open('F:/matrisomedb2.0/glob_prot_ptm_ind_dict.p','rb'))
    # converting ndarray to list in id_ptm_idx_dict for json feeding
    clean_id_ptm_idx_dict = {}
    for prot in id_ptm_idx_dict:
        dict1 = {}
        ptm_idx_dict = id_ptm_idx_dict[prot]
        for ptm in ptm_idx_dict:
            index_list = ptm_idx_dict[ptm].tolist()
            dict1[ptm] = index_list
        clean_id_ptm_idx_dict[prot] = dict1

    # for loop to generate HTMLs starts here!

    for protein in id_ptm_idx_dict:
    # protein = 'O15230'
        if os.path.exists('F:/matrisomedb2.0/mat_alphafold_pdbs/AF-' + protein + '-F1-model_v1.pdb'):
            time_start = time.time()
            pdb_file = 'F:/matrisomedb2.0/mat_alphafold_pdbs/AF-' + protein + '-F1-model_v1.pdb'
            output_path = 'F:/matrisomedb2.0/3dcov_htmls/' + protein + '_3dcov.html'
            ptm_idx_dict = clean_id_ptm_idx_dict[protein]
            regex_dict = {each: regex_dict[each] for each in regex_dict if each in ptm_idx_dict}
            json_id_ptm_dict_dict = {protein: clean_id_ptm_idx_dict[protein]}
            # pdb_file = r'D:\data\alphafold_pdb\UP000000589_10090_MOUSE/AF-Q99MQ5-F1-model_v1.pdb'

            try:
                dump_dict = show_cov_3d_v2(protein, 'AAA', pdb_file, np.array(id_freq_array_dict[protein]).astype(np.int32),
                                           id_ptm_idx_dict=clean_id_ptm_idx_dict[protein],
                                           regex_color_dict=regex_dict,
                                           background_color=str(16777215))

                view_file = open('./view.html')
                soup = Soup(view_file.read(), 'lxml')
                js_data = soup.find(id='data')
                js_data.string = 'let id_ptm_idx_dict = ' + json.dumps(json_id_ptm_dict_dict) + ';' + \
                                 'let regex_dict = ' + json.dumps(regex_dict) + ';' + \
                                 'let background_color = ' + str(16777215) + ';' + \
                                 'let pq = ' + json.dumps(pq_dict[protein]) + ';'
                pdb_str = soup.find(id='glmol01_src')
                pdb_str.string = dump_dict['pbdstr']
                ret = soup.find(id='glmol01_rep')
                ret.string = dump_dict['ret']
                output_file = open(output_path, "w")
                # output_file = open('P11276_test.html','w')
                output_file.write(str(soup))
                output_file.close()
                print(f'total time{time.time() - time_start}')
            except KeyError:
                continue
