# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from collections import OrderedDict
import copy
import re
from sqlalchemy.exc import DBAPIError

import SequenceDataORM_updt as sqd


def gd_line_to_dict(gd_line):
    if gd_line[0] == '#':
        raise ValueError("Comment lines aren't mutations or evidence for one")
    entries = gd_line.strip('\n').split('\t')
    gd_dict = OrderedDict()
    gd_dict['type'] = entries[0]
    gd_dict['id'] = int(entries[1])
    try:
        gd_dict['parent_id'] = int(entries[2])
    except ValueError:
        gd_dict['parent_id'] = entries[2]
    # types of genome diff mutations are SNP, SUB, DEL, INS, MOB, AMP, CON, INV
    # types of genome diff evidence are RA, MC, JC, UN
    # this doesn't handle validation types
    type_ = gd_dict['type']
    if type_ == 'SNP':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['new_seq'] = entries[5]
        next_idx = 6
    elif type_ == 'SUB':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        gd_dict['new_seq'] = entries[6]
        next_idx = 7
    elif type_ == 'DEL':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        next_idx = 6
    elif type_ == 'INS':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = entries[4]
        gd_dict['new_seq'] = entries[5]
        next_idx = 6
    elif type_ == 'MOB':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['repeat_name'] = entries[5]
        gd_dict['strand'] = entries[6]
        gd_dict['duplication_size'] = int(entries[7])
        next_idx = 8
    elif type_ == 'AMP':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        gd_dict['new_copy_number'] = int(entries[6])
        next_idx = 7
    elif type_ == 'CON':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        gd_dict['region'] = entries[6]
        next_idx = 7
    elif type_ == 'INV':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['size'] = int(entries[5])
        next_idx = 6
    elif type_ == 'RA':
        gd_dict['seq_id'] = entries[3]
        gd_dict['position'] = int(entries[4])
        gd_dict['insert_position'] = int(entries[5])
        gd_dict['ref_base'] = entries[6]
        gd_dict['new_base'] = entries[7]
        next_idx = 8
    elif type_ == 'MC':
        gd_dict['seq_id'] = entries[3]
        gd_dict['start'] = int(entries[4])
        gd_dict['end'] = int(entries[5])
        gd_dict['start_range'] = int(entries[6])
        gd_dict['end_range'] = int(entries[7])
        next_idx = 8
    elif type_ == 'JC':
        gd_dict['side_1_seq_id'] = entries[3]
        gd_dict['side_1_position'] = int(entries[4])
        gd_dict['side_1_strand'] = entries[5]
        gd_dict['side_2_seq_id'] = entries[6]
        gd_dict['side_2_position'] = int(entries[7])
        gd_dict['side_2_strand'] = entries[8]
        gd_dict['overlap'] = int(entries[9])
        next_idx = 10
    elif type_ == 'UN':
        gd_dict['seq_id'] = entries[3]
        gd_dict['start'] = int(entries[4])
        gd_dict['end'] = int(entries[5])
        next_idx = 6
    else:
        raise RuntimeError('unknown .gd entry type')
    for entry in entries[next_idx:]:
        try:
            key, value = entry.split('=')
            gd_dict[key] = value
        except ValueError:
            print('invalid key-value entry in .gd line')
    return gd_dict


def RA_dict_to_SNP_dicts(sample, gd_dict):
    # first unpack what we need to look up or make the entry in the genes table
    gene = {}
    gene['name'] = gd_dict.get('gene_name')
    if gene['name'] is None:
        return None
    gene['description'] = gd_dict.get('gene_product')
    gene['intergenic'] = (gd_dict['snp_type'] == 'intergenic')
    # next unpack what we need to look up or make the SNP_mutations entry
    snp_mutation = {}
    snp_mutation['chr_position'] = gd_dict.get('position')
    snp_mutation['ref_base'] = gd_dict.get('ref_base')
    snp_mutation['new_base'] = gd_dict.get('new_base')
    if (type(snp_mutation['chr_position']) is not int) or \
       (snp_mutation['ref_base'] is None) or \
       (snp_mutation['new_base'] is None):
        return None
    snp_mutation['gene'] = gene['name']
    snp_mutation['gene_pos'] = gd_dict.get('gene_position')
    snp_mutation['ref_aa'] = gd_dict.get('aa_ref_seq')
    snp_mutation['new_aa'] = gd_dict.get('aa_new_seq')
    if gene['intergenic'] and (snp_mutation['gene_pos'] is not None):
        nums = re.findall('[-+]?[0-9]+', snp_mutation['gene_pos'])
        if len(nums) == 2:
            snp_mutation['intergenic_left'] = int(nums[0])
            snp_mutation['intergenic_right'] = int(nums[1])
        else:
            snp_mutation['intergenic_left'] = None
            snp_mutation['intergenic_right'] = None
    else:
        snp_mutation['intergenic_left'] = None
        snp_mutation['intergenic_right'] = None
    # finally unpack what we need to look up or make the SNP_evidence entry
    snp_evidence = {}
    snp_evidence['sample'] = sample
    snp_evidence['chr_position'] = snp_mutation['chr_position']
    snp_evidence['ref_base'] = snp_mutation['ref_base']
    snp_evidence['new_base'] = snp_mutation['new_base']
    snp_evidence['frequency'] = gd_dict.get('frequency')
    snp_evidence['bias_e_value'] = gd_dict.get('bias_e_value')
    snp_evidence['bias_p_value'] = gd_dict.get('bias_p_value')
    snp_evidence['consensus_score'] = gd_dict.get('consensus_score')
    snp_evidence['fisher_strand_p_value'] = \
        gd_dict.get('fisher_strand_p_value')
    snp_evidence['ref_cov'] = gd_dict.get('ref_cov')
    snp_evidence['new_cov'] = gd_dict.get('new_cov')
    snp_evidence['total_cov'] = gd_dict.get('total_cov')
    return gene, snp_mutation, snp_evidence


def commit_SNP_dicts(gene, snp_mutation, snp_evidence, session):
    # check if the gene is already in the database. If not, add it.
    gene_matches = session.query(sqd.Gene).filter(sqd.Gene.name ==
                                                  gene['name'])
    if gene_matches.count() == 0:
        session.add(sqd.Gene(**gene))
    # check if the snp_mutation is already in the database. If not, add it.
    mu_query = session.query(sqd.SNP_Mutation)
    mu_matches = mu_query.filter(sqd.SNP_Mutation.chr_position ==
                                 snp_mutation['chr_position'])
    mu_matches = mu_matches.filter(sqd.SNP_Mutation.ref_base ==
                                   snp_mutation['ref_base'])
    mu_matches = mu_matches.filter(sqd.SNP_Mutation.new_base ==
                                   snp_mutation['new_base'])
    if mu_matches.count() == 0:
        session.add(sqd.SNP_Mutation(**snp_mutation))
    # check if the snp_evidence is already in the database. If not, add it.
    ev_query = session.query(sqd.SNP_Evidence)
    ev_matches = ev_query.filter(sqd.SNP_Evidence.sample ==
                                 snp_evidence['sample'])
    ev_matches = ev_matches.filter(sqd.SNP_Evidence.chr_position ==
                                   snp_evidence['chr_position'])
    ev_matches = ev_matches.filter(sqd.SNP_Evidence.ref_base ==
                                   snp_evidence['ref_base'])
    ev_matches = ev_matches.filter(sqd.SNP_Evidence.new_base ==
                                   snp_evidence['new_base'])
    if ev_matches.count() == 0:
        session.add(sqd.SNP_Evidence(**snp_evidence))
    # try committing the changes to the database. If it fails, rollback.
    try:
        session.commit()
    except DBAPIError:
        session.rollback()


def gdMC_dict_to_MC_dict(sample, gd_dict):
    mc = {}
    mc['sample'] = sample
    mc['start'] = gd_dict['start']
    mc['end'] = gd_dict['end']
    mc['start_range'] = gd_dict['start_range']
    mc['end_range'] = gd_dict['end_range']
    mc['left_inside_cov'] = gd_dict.get('left_inside_cov')
    mc['left_outside_cov'] = gd_dict.get('left_outside_cov')
    mc['right_inside_cov'] = gd_dict.get('right_inside_cov')
    mc['right_outside_cov'] = gd_dict.get('right_outside_cov')
    mc['gene'] = gd_dict.get('gene_name')
    return mc


def add_MC_dict(mc, session):
    mc_matches = (session.query(sqd.MC_Evidence)
                         .filter(sqd.MC_Evidence.sample == mc['sample'])
                         .filter(sqd.MC_Evidence.start == mc['start'])
                         .filter(sqd.MC_Evidence.end == mc['end']))
    if mc_matches.count() == 0:
        session.add(sqd.MC_Evidence(**mc))


def gdMOB_dict_to_MOB_dicts(sample, gd_dict):
    mob = {}
    mob['chr_position'] = gd_dict['position']
    mob['repeat_name'] = gd_dict['repeat_name']
    if gd_dict['strand'] == '1':
        mob['strand'] = True
    elif gd_dict['strand'] == '-1':
        mob['strand'] = False
    else:
        raise ValueError('strand must be either "1" or "-1"')
    mob['duplication_size'] = int(gd_dict['duplication_size'])
    mob_mut = copy.deepcopy(mob)
    mob_ev = copy.deepcopy(mob)
    mob_mut['gene'] = gd_dict['gene_name']
    mob_mut['repeat_size'] = int(gd_dict['repeat_size'])
    mob_ev['sample'] = sample
    mob_ev['frequency'] = float(gd_dict['frequency'])
    return mob_mut, mob_ev


def add_MOB_dicts(mob_mut, mob_ev, session):
    mob_mut_matches = (session.query(sqd.MOB_Mutation)
                              .filter(sqd.MOB_Mutation.chr_position
                                      == mob_mut['chr_position'])
                              .filter(sqd.MOB_Mutation.repeat_name
                                      == mob_mut['repeat_name'])
                              .filter(sqd.MOB_Mutation.strand
                                      == mob_mut['strand'])
                              .filter(sqd.MOB_Mutation.duplication_size
                                      == mob_mut['duplication_size']))
    if mob_mut_matches.count() == 0:
        session.add(sqd.MOB_Mutation(**mob_mut))
    mob_ev_matches = (session.query(sqd.MOB_Evidence)
                             .filter(sqd.MOB_Evidence.chr_position
                                     == mob_ev['chr_position'])
                             .filter(sqd.MOB_Evidence.repeat_name
                                     == mob_ev['repeat_name'])
                             .filter(sqd.MOB_Evidence.strand
                                     == mob_ev['strand'])
                             .filter(sqd.MOB_Evidence.duplication_size
                                     == mob_ev['duplication_size'])
                             .filter(sqd.MOB_Evidence.sample
                                     == mob_ev['sample']))
    if mob_ev_matches.count() == 0:
        session.add(sqd.MOB_Evidence(**mob_ev))


def gdDEL_dict_to_del_dicts(sample, gd_dict):
    del_mut, del_ev = {}, {}
    del_mut['chr_position'] = gd_dict['position']
    del_mut['size'] = gd_dict['size']
    del_mut['gene'] = gd_dict['gene_name']
    del_ev['sample'] = sample
    del_ev['chr_position'] = gd_dict['position']
    del_ev['size'] = gd_dict['size']
    del_ev['frequency'] = float(gd_dict['frequency'])
    return del_mut, del_ev


def check_merge_RAdel_dicts(*del_dicts):
    if len(del_dicts) <= 1:
        return True
    del_list = list(copy.deepcopy(del_dicts))
    del_list.sort(key=lambda x: x['position'])
    curr = del_list[0]['position']
    new_cov = del_list[0]['new_cov']
    for item in del_dicts[1:]:
        if item['position'] == curr + 1 and item['new_cov'] == new_cov:
            curr = item['position']
        else:
            return False
    return True


def merge_RAdel_dicts_to_del_dicts(sample, *del_dicts):
    if check_merge_RAdel_dicts(*del_dicts):
        del_list = list(copy.deepcopy(del_dicts))
        del_list.sort(key=lambda x: x['position'])
        merged_del_mut = {}
        merged_del_mut['chr_position'] = del_list[0]['position']
        merged_del_mut['size'] = len(del_list)
        merged_del_mut['gene'] = del_list[0]['gene_name']
        merged_del_ev = {}
        merged_del_ev['sample'] = sample
        merged_del_ev['chr_position'] = merged_del_mut['chr_position']
        merged_del_ev['size'] = merged_del_mut['size']
        merged_del_ev['frequency'] = sum(float(del_dict['frequency'])
                                         for del_dict
                                         in del_list)/len(del_list)
    return merged_del_mut, merged_del_ev


def add_DEL_dicts(del_mut, del_ev, session):
    del_mut_matches = (session.query(sqd.DEL_Mutation)
                              .filter(sqd.DEL_Mutation.chr_position
                                      == del_mut['chr_position'])
                              .filter(sqd.DEL_Mutation.size
                                      == del_mut['size']))
    if del_mut_matches.count() == 0:
        session.add(sqd.DEL_Mutation(**del_mut))
    del_ev_matches = (session.query(sqd.DEL_Evidence)
                             .filter(sqd.DEL_Evidence.sample
                                     == del_ev['sample'])
                             .filter(sqd.DEL_Evidence.chr_position
                                     == del_ev['chr_position'])
                             .filter(sqd.DEL_Evidence.size
                                     == del_ev['size']))
    if del_ev_matches.count() == 0:
        session.add(sqd.DEL_Evidence(**del_ev))


def gdINS_dict_to_ins_dicts(sample, gd_dict):
    ins_mut, ins_ev = {}, {}
    ins_mut['chr_position'] = gd_dict['position']
    ins_mut['new_seq'] = gd_dict['new_seq']
    ins_mut['gene'] = gd_dict['gene_name']
    ins_ev['sample'] = sample
    ins_ev['chr_position'] = gd_dict['position']
    ins_ev['new_seq'] = gd_dict['new_seq']
    ins_ev['frequency'] = float(gd_dict['frequency'])
    return ins_mut, ins_ev


def check_merge_RAins_dicts(*ins_dicts):
    if len(ins_dicts) <= 1:
        return True
    ins_list = list(copy.deepcopy(ins_dicts))
    ins_list.sort(key=lambda x: x['position'])
    curr = ins_list[0]['position']
    new_cov = ins_list[0]['new_cov']
    for item in ins_dicts[1:]:
        if item['position'] == curr and item['new_cov'] == new_cov:
            pass
        else:
            return False
    return True


def merge_RAins_dicts_to_ins_dicts(sample, *ins_dicts):
    if check_merge_RAins_dicts(*ins_dicts):
        ins_list = list(copy.deepcopy(ins_dicts))
        ins_list.sort(key=lambda x: x['position'])
        merged_ins_mut = {}
        merged_ins_mut['chr_position'] = ins_list[0]['position']
        merged_ins_mut['new_seq'] = ''.join(ins_dict['new_base']
                                            for ins_dict in ins_list)
        merged_ins_mut['gene'] = ins_list[0]['gene_name']
        merged_ins_ev = {}
        merged_ins_ev['sample'] = sample
        merged_ins_ev['chr_position'] = merged_ins_mut['chr_position']
        merged_ins_ev['new_seq'] = merged_ins_mut['new_seq']
        merged_ins_ev['frequency'] = sum(float(ins_dict['frequency'])
                                         for ins_dict
                                         in ins_list)/len(ins_list)
    return merged_ins_mut, merged_ins_ev


def add_INS_dicts(ins_mut, ins_ev, session):
    ins_mut_matches = (session.query(sqd.INS_Mutation)
                              .filter(sqd.INS_Mutation.chr_position
                                      == ins_mut['chr_position'])
                              .filter(sqd.INS_Mutation.new_seq
                                      == ins_mut['new_seq']))
    if ins_mut_matches.count() == 0:
        session.add(sqd.INS_Mutation(**ins_mut))
    ins_ev_matches = (session.query(sqd.INS_Evidence)
                             .filter(sqd.INS_Evidence.sample
                                     == ins_ev['sample'])
                             .filter(sqd.INS_Evidence.chr_position
                                     == ins_ev['chr_position'])
                             .filter(sqd.INS_Evidence.new_seq
                                     == ins_ev['new_seq']))
    if ins_ev_matches.count() == 0:
        session.add(sqd.INS_Evidence(**ins_ev))
