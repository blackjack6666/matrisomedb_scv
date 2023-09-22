def freq_ptm_index_gen_batch_v2(psm_dict, protein_dict, pdbs, regex_dict=None):
    """
    :param psm_dict: psms, output from js, {"unlabeled":[psm1, psm2], "group1":[...], "group2":[...]}
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

    psm_list = set([psm for v in psm_dict.values() for psm in v])
    psm_group_dict = {psm:group for group in psm_dict if group != "unlabeled" for psm in psm_dict[group]}

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
                    # print (ptm)
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
            id_ptm_idx_dict[id_list[i]] = {ptm.replace('\[','[').replace('\]',']').replace('\.','.'): np.array(
                np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[0]).tolist()
                                           for ptm in ptm_index_line_dict}

    print(time.time() - time_start)

    return id_freq_array_dict, id_ptm_idx_dict, [heappop(h) for i in range(len(h))][::-1]