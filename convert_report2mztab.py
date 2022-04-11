#!/usr/bin/env python
from doctest import OutputChecker
from pickle import APPEND
from sqlite3 import DatabaseError
import pandas as pd
import click
import os
import re
import numpy as np
import time
import sys
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@click.command("convert_report2mztab")
@click.option("--mztab", "-mz")
@click.option("--diann_report", "-r")
@click.option("--unimod_csv", "-u")
@click.option("--pg_matrix", "-pg")
@click.option("--pr_matrix", "-pr")
@click.option("--unique_matrix", "-un")
@click.option("--openms", "-o")
@click.option("--exp_design", "-e")
@click.option("--fasta", "-f")
@click.pass_context

def convert_report2mztab(ctx, mztab, diann_report, unimod_csv, pg_matrix, pr_matrix, unique_matrix, openms, exp_design, fasta):
    start = time.time()
    report = pd.read_csv(diann_report, sep = "\t", header = 0, dtype = 'str')  
    pg = pd.read_csv(pg_matrix, sep = '\t', header = 0, dtype = 'str')
    pr = pd.read_csv(pr_matrix, sep = '\t', header = 0, dtype = 'str')  
    unimod_data = pd.read_csv(unimod_csv, sep = ",", header = 0, dtype = 'str')  
    unique = pd.read_csv(unique_matrix, sep = '\t', header = 0, dtype = 'str')  
    opms = pd.read_csv(openms, sep = '\t', header = 0, dtype = 'str') 
    with open(exp_design, 'r') as f:
        data = f.readlines()
        empty_row = data.index('\n')
        f_table = [i.replace("\n", '').split("\t") for i in data[1:empty_row]]
        f_header = data[0].replace("\n", "").split("\t")
        f_table = pd.DataFrame(f_table, columns=f_header)
        f_table.loc[:,"run"] = f_table.apply(lambda x: os.path.basename(x["Spectra_Filepath"].split(".")[-2]), axis=1)
        f_table.loc[:, "ms_run"] = f_table.apply(lambda x: x["Fraction_Group"], axis=1)
        f_table.loc[:, "study_variable"] = f_table.apply(lambda x: x["Sample"], axis=1)
        s_table = [i.replace("\n", '').split("\t") for i in data[empty_row + 1:]][1:]
        s_header = data[empty_row + 1].replace("\n", "").split("\t")
        s_DataFrame = pd.DataFrame(s_table, columns=s_header)
    index_ref = f_table
    index_ref.loc[:, 'ms_run'] = index_ref.loc[:, 'ms_run'].astype('int')
    index_ref.loc[:, 'study_variable'] = index_ref.loc[:, 'study_variable'].astype('int')
    report.loc[:, 'ms_run'] = report.apply(lambda x: get_ms_run(x['Run'], index_ref), axis=1, result_type="expand")
    report.loc[:, 'ms_run'] = report.loc[:, 'ms_run'].astype('str')
    #charge = sys.argv[2]
    #missed_cleavages = sys.argv[3]
    (MTD, database) = mztab_MTD(index_ref, opms, report, unimod_data, fasta)
    RPH = mztab_RPH(report, pg, unique, index_ref, database)
    PEH = mztab_PEH(report, pr, unique, unimod_data, index_ref, database)
    PSH = mztab_PSH(unique, unimod_data, report)
    MTD.loc['', :] = ''
    RPH.loc[len(RPH) + 1, :] = ''
    PEH.loc[len(PEH) + 1, :] = ''
    with open(mztab, 'w', newline = '') as f:
        MTD.to_csv(f, mode='f+', index = False, header = False)
        RPH.to_csv(f, mode='f+', index = False, header = True)
        PEH.to_csv(f, mode='f+', index = False, header = True)
        PSH.to_csv(f, mode='f+', index = False, header = True)
    end = time.time()
    print("Total run time:{t}".format(t=end - start))

def mztab_MTD(index_ref, opms, report, unimod_data, fasta):
    start = time.time()
    FragmentMassTolerance = opms['FragmentMassTolerance'].values[0]
    FragmentMassToleranceUnit = opms['FragmentMassToleranceUnit'].values[0]
    PrecursorMassTolerance = opms['PrecursorMassTolerance'].values[0]
    PrecursorMassToleranceUnit = opms['PrecursorMassToleranceUnit'].values[0]
    Enzyme = opms['Enzyme'].values[0]
    FixedModifications = opms['FixedModifications'].values[0]
    VariableModifications = opms['VariableModifications'].values[0]
    fixed_mods = mod_inf(FixedModifications, unimod_data)
    variable_mods = mod_inf(VariableModifications, unimod_data)
    out_mztab_MTD = pd.DataFrame()
    out_mztab_MTD.loc[1, 'mzTab-version'] = '1.0.0'
    out_mztab_MTD.loc[1, 'mzTab-mode'] = 'Summary'
    out_mztab_MTD.loc[1, 'mzTab-type'] = 'Quantification'
    out_mztab_MTD.loc[1, 'title'] = 'ConsensusMap export from OpenMS'
    out_mztab_MTD.loc[1, 'description'] = 'OpenMS export from consensusXML'
    out_mztab_MTD.loc[1, 'protein_search_engine_score[1]'] = '[MS, MS:null, DIANN:Protein.Q.Value, ]'
    out_mztab_MTD.loc[1, 'peptide_search_engine_score[1]'] = '[MS, MS:null, DIANN:Q.Value, ]'
    out_mztab_MTD.loc[1, 'psm_search_engine_score[1]'] = '[, , DIANN null, ]'
    out_mztab_MTD.loc[1, 'software[1]'] = '[MS, MS:null, DIANN, Release (v1.8)]'
    out_mztab_MTD.loc[1, 'software[1]-setting[1]'] = fasta
    out_mztab_MTD.loc[1, 'software[1]-setting[2]'] = 'db_version:null'
    out_mztab_MTD.loc[1, 'software[1]-setting[3]'] = 'fragment_mass_tolerance:' + FragmentMassTolerance
    out_mztab_MTD.loc[1, 'software[1]-setting[4]'] = 'fragment_mass_tolerance_unit:' + FragmentMassToleranceUnit
    out_mztab_MTD.loc[1, 'software[1]-setting[5]'] = 'precursor_mass_tolerance:' + PrecursorMassTolerance
    out_mztab_MTD.loc[1, 'software[1]-setting[6]'] = 'precursor_mass_tolerance_unit:' + PrecursorMassToleranceUnit
    out_mztab_MTD.loc[1, 'software[1]-setting[7]'] = 'enzyme:' + Enzyme
    out_mztab_MTD.loc[1, 'software[1]-setting[8]'] = 'enzyme_term_specificity:full'
    #out_mztab_MTD.loc[1, 'software[1]-setting[9]'] = 'charges:' + charge
    #out_mztab_MTD.loc[1, 'software[1]-setting[10]'] = 'missed_cleavages:' + missed_cleavages
    out_mztab_MTD.loc[1, 'software[1]-setting[11]'] = 'fixed_modifications:' + FixedModifications
    out_mztab_MTD.loc[1, 'software[1]-setting[12]'] = 'variable_modifications:' + VariableModifications
    for i in range(1, len(fixed_mods) // 2 + 1):
        out_mztab_MTD.loc[1, 'fixed_mod[' + str(i) + ']'] = fixed_mods[2 * i - 2]
        out_mztab_MTD.loc[1, 'fixed_mod[' + str(i) + ']-site'] = fixed_mods[2 * i - 1]
        out_mztab_MTD.loc[1, 'fixed_mod[' + str(i) + ']-position'] = 'Anywhere'
    for i in range(1, len(variable_mods) // 2 + 1):
        out_mztab_MTD.loc[1, 'variable_mods[' + str(i) + ']'] = variable_mods[2 * i - 2]
        out_mztab_MTD.loc[1, 'variable_mods[' + str(i) + ']-site'] = variable_mods[2 * i - 1]
        out_mztab_MTD.loc[1, 'variable_mods[' + str(i) + ']-position'] = 'Anywhere'
    out_mztab_MTD.loc[1, 'quantification_method'] = '[MS, MS:1001834, LC-MS label-free quantitation analysis, ]'
    out_mztab_MTD.loc[1, 'protein-quantification_unit'] = '[, , Abundance, ]'
    out_mztab_MTD.loc[1, 'peptide-quantification_unit'] = '[, , Abundance, ]'
    for i in range(1, max(index_ref['ms_run']) + 1):
        out_mztab_MTD.loc[1, 'ms_run[' + str(i) + ']-format'] = '[MS, MS:1000584, mzML file, ]'
        out_mztab_MTD.loc[1, 'ms_run[' + str(i) + ']-location'] = 'file://' + index_ref[index_ref['ms_run'] == i]['Spectra_Filepath'].values[0]
        out_mztab_MTD.loc[1, 'ms_run[' + str(i) + ']-id_format'] = '[MS, MS:1000777, spectrum identifier nativeID format, ]'
        out_mztab_MTD.loc[1, 'assay[' + str(i) + ']-quantification_reagent'] = '[MS, MS:1002038, unlabeled sample, ]'
        out_mztab_MTD.loc[1, 'assay[' + str(i) + ']-ms_run_ref'] = 'ms_run[' + str(i) + ']'
    for i in range(1, max(index_ref['study_variable']) + 1):
        study_variable = []
        for j in list(index_ref[index_ref['study_variable'] == i]['ms_run'].values):
            study_variable.append('assay[' + str(j) + ']')
        out_mztab_MTD.loc[1, 'study_variable[' + str(i) + ']-assay_refs'] = ','.join(study_variable)
        out_mztab_MTD.loc[1, 'study_variable[' + str(i) + ']-description'] = 'no description given'   
    out_mztab_MTD.loc[2, :] = 'MTD'
    col = list(out_mztab_MTD.columns)
    row = list(out_mztab_MTD.index)
    out_mztab_MTD_T = pd.DataFrame(out_mztab_MTD.values.T, index = col, columns = row)
    out_mztab_MTD_T.columns = ['inf', 'index']
    out_mztab_MTD_T.insert(0, 'title', out_mztab_MTD_T.index)
    index = out_mztab_MTD_T.loc[:, 'index']
    out_mztab_MTD_T.drop(labels = ['index'], axis = 1, inplace = True)
    out_mztab_MTD_T.insert(0, 'index', index)
    database = re.split('\\\\|/', fasta)[-1].split('.')[-2]
    print('<--      MTD has been processed !        -->')
    end = time.time()
    print("MTD run time:{t}".format(t=end - start))
    return out_mztab_MTD_T, database

def mztab_RPH(report, pg, unique, index_ref, database):
    start = time.time()
    file = list(pg.columns[5:])
    col = {}
    for i in file:
        col[i] = 'protein_abundance_assay[' + str(index_ref[index_ref['run'] == i.split("/")[-1].split(".")[-2]]['ms_run'].values[0]) + ']'
    pg = pg.rename(columns = col)
    out_mztab_RPH = pd.DataFrame()
    out_mztab_RPH = pg.drop(['Protein.Group', 'Protein.Names'], axis = 1)
    out_mztab_RPH = out_mztab_RPH.rename(columns = {'Protein.Ids':'accession', 'First.Protein.Description':'description'})
    out_mztab_RPH.loc[:, 'PRH'] = 'PRT'
    index = out_mztab_RPH.loc[:, 'PRH']
    out_mztab_RPH.drop(labels = ['PRH'], axis = 1, inplace = True)
    out_mztab_RPH.insert(0, 'RPH', index)
    out_mztab_RPH.loc[:, 'database'] = database
    null_col = ['taxid', 'species', 'database_version', 'search_engine', 'protein_coverage', 'opt_global_Posterior_Probability_score', 'opt_global_nr_found_peptides', 'opt_global_cv_PRIDE:0000303_decoy_hit']
    for i in null_col:
        out_mztab_RPH.loc[:, i] = 'null'
    out_mztab_RPH.loc[:, 'ambiguity_members'] = out_mztab_RPH.loc[:, 'accession']
    out_mztab_RPH[["modifiedSequence", "best_search_engine_score[1]"]] = out_mztab_RPH.apply(lambda x: RPHfrom_report(report, x["accession"]), axis=1, result_type="expand")
    out_mztab_RPH.loc[:, 'modifications'] = out_mztab_RPH.apply(lambda x: find_modification(x['modifiedSequence']), axis = 1, result_type = 'expand') 
    max_assay = max(index_ref['ms_run'])
    max_study_variable = max(index_ref['study_variable'])
    ## Change the type of each abundance_assay
    for i in range(1, max_assay + 1):
        out_mztab_RPH.loc[:, 'protein_abundance_assay[' + str(i) +']'] = out_mztab_RPH.loc[:, 'protein_abundance_assay[' + str(i) +']'].astype("float")
    for i in range(1, max_study_variable + 1):
        out_mztab_RPH.loc[:, 'protein_abundance_study_variable[' + str(i) + ']'] = 0
        assay_num = max(index_ref[index_ref['study_variable'] == i]['ms_run'].values)
        for j in list(index_ref[index_ref['study_variable'] == i]['ms_run'].values):
            out_mztab_RPH.loc[:, 'protein_abundance_study_variable[' + str(i) + ']'] += out_mztab_RPH.loc[:, 'protein_abundance_assay[' + str(j) +']']
        out_mztab_RPH.loc[:, 'protein_abundance_study_variable[' + str(i) + ']'] = out_mztab_RPH.loc[:, 'protein_abundance_study_variable[' + str(i) + ']'] / assay_num
        out_mztab_RPH[['protein_abundance_stdev_study_variable[' + str(i) + ']', 'protein_abundance_std_error_study_variable[' + str(i) + ']']] = out_mztab_RPH.apply(lambda x:('null', 'null'),axis = 1, result_type = 'expand')
    out_mztab_RPH.loc[:, "opt_global_result_type"] = out_mztab_RPH.apply(lambda x: classify_protein(x["Genes"], unique["Genes"], 'single_protein', 'indistinguishable_protein_group'), axis = 1, result_type = 'expand')
    out_mztab_RPH = out_mztab_RPH.drop(['Genes', 'modifiedSequence'], axis = 1)
    out_mztab_RPH.fillna('null', inplace=True)  
    out_mztab_RPH.to_csv('./mztab_RPH.csv', sep=',', index=False)
    print('<--      RPH has been processed !        -->')
    end = time.time()
    print("RPH run time:{t}".format(t=end - start))
    return out_mztab_RPH

def mztab_PEH(report, pr, unique, unimod_data, index_ref, database):
    start = time.time()
    file = list(pr.columns[10:])
    col = {}
    for i in file:
        col[i] = 'protein_abundance_assay[' + str(index_ref[index_ref['run'] == i.split("/")[-1].split(".")[-2]]['ms_run'].values[0]) + ']'
    pr = pr.rename(columns = col)
    out_mztab_PEH = pd.DataFrame()
    out_mztab_PEH = pr.drop(['Protein.Group', 'Protein.Names', 'First.Protein.Description', 'Proteotypic'], axis = 1)
    out_mztab_PEH = out_mztab_PEH.rename(columns = {'Stripped.Sequence':'sequence', 'Protein.Ids':'accession', 
                                                    'Modified.Sequence':'opt_global_cv_MS:1000889_peptidoform_sequence', 'Precursor.Charge':'charge'})
    out_mztab_PEH.loc[:, 'PEH'] = 'PEP'
    index = out_mztab_PEH.loc[:, 'PEH']
    out_mztab_PEH.drop(labels = ['PEH'], axis = 1, inplace = True)
    out_mztab_PEH.insert(0, 'PEH', index)
    out_mztab_PEH.loc[:, 'database'] = database
    out_mztab_PEH.loc[:, 'modifications'] = out_mztab_PEH.apply(lambda x: find_modification(x['opt_global_cv_MS:1000889_peptidoform_sequence']), axis = 1, result_type = 'expand')
    out_mztab_PEH.loc[:, 'opt_global_cv_MS:1000889_peptidoform_sequence'] = out_mztab_PEH.apply(lambda x: convert_modification(x["opt_global_cv_MS:1000889_peptidoform_sequence"], unimod_data), axis=1)
    out_mztab_PEH.loc[:, 'unique'] = out_mztab_PEH.apply(lambda x: classify_protein(x['Genes'], unique["Genes"], '1', '0'), axis=1, result_type="expand")
    null_col = ['database_version', 'search_engine', 'retention_time_window']
    for i in null_col:
        out_mztab_PEH.loc[:, i] = 'null'
    out_mztab_PEH.loc[:, 'opt_global_cv_MS:1002217_decoy_peptide'] = '0'
    # average value of each study_variable
    max_assay = max(index_ref['ms_run'])
    max_study_variable = max(index_ref['study_variable'])
    for i in range(1, max_assay + 1):                                                                                                                 
        out_mztab_PEH.loc[:, 'search_engine_score[1]_ms_run[' + str(i) + ']'] = out_mztab_PEH.apply(lambda x: find_each_params(report, x['Precursor.Id'], index_ref, i, "Q.Value", 0), axis = 1, result_type = 'expand')
    for i in range(1, max_study_variable + 1):                                                                                                                 
        out_mztab_PEH.loc[:, 'peptide_abundance_study_variable[' + str(i) + ']'] = out_mztab_PEH.apply(lambda x: find_each_params(report, x['Precursor.Id'], index_ref, i, "Precursor.Normalised", 1), axis = 1, result_type = 'expand')
        out_mztab_PEH.loc[:, 'opt_global_mass_to_charge_study_variable[' + str(i) + ']'] = out_mztab_PEH.apply(lambda x: find_each_params(report, x['Precursor.Id'], index_ref, i, "Precursor.Mz", 1), axis = 1, result_type = 'expand')
        out_mztab_PEH.loc[:, 'opt_global_retention_time_study_variable[' + str(i) + ']'] = out_mztab_PEH.apply(lambda x: find_each_params(report, x['Precursor.Id'], index_ref, i, "RT", 1), axis = 1, result_type = 'expand')
        out_mztab_PEH[['peptide_abundance_stdev_study_variable[' + str(i) + ']', 'peptide_abundance_std_error_study_variable[' + str(i) + ']']] = out_mztab_PEH.apply(lambda x:('null', 'null'), axis = 1, result_type = 'expand')  
    out_mztab_PEH[["best_search_engine_score[1]", "retention_time", "opt_global_q-value", "opt_global_SpecEValue_score", "mass_to_charge"]] = out_mztab_PEH.apply(lambda x: 
                    PEHfrom_report(report, x["Precursor.Id"], ["Q.Value", "RT", "Global.Q.Value", "Lib.Q.Value", "Precursor.Mz"]), axis=1, result_type="expand")
    out_mztab_PEH[['opt_global_feature_id', 'spectra_ref']] = out_mztab_PEH.apply(lambda x:('null', 'null'),axis = 1, result_type = 'expand')
    ##opt_global_mass_to_charge_study_variable[n], opt_global_retention_time_study_variable[n]
    out_mztab_PEH = out_mztab_PEH.drop(['Precursor.Id', 'Genes'], axis = 1)
    out_mztab_PEH.fillna('null', inplace = True)  
    out_mztab_PEH.to_csv('./mztab_PEH.csv', sep=',', index=False)
    print('<--      PEH has been processed !        -->')
    end = time.time()
    print("PEH run time:{t}".format(t=end - start))
    return out_mztab_PEH

def mztab_PSH(unique, unimod_data, report, database):
    start = time.time()
    out_mztab_PSH = pd.DataFrame()
    out_mztab_PSH = report[['Stripped.Sequence', 'Protein.Ids', 'Genes', 'Q.Value', 'RT', 'Precursor.Charge', 'Precursor.Mz', 'Modified.Sequence', 'PEP', 'Global.Q.Value', 'Global.Q.Value']]
    out_mztab_PSH.columns = ['sequence', 'accession', 'Genes', 'search_engine_score[1]', 'retention_time', 'charge', 'exp_mass_to_charge', 'opt_global_cv_MS:1000889_peptidoform_sequence', 
                            'opt_global_SpecEValue_score', 'opt_global_q-value', 'opt_global_q-value_score']
    out_mztab_PSH.loc[:, 'spectra_ref'] = report.apply(lambda x: 'ms_run[' + x['ms_run'] + ']:index=' + x['Lib.Index'], axis=1, result_type="expand")
    out_mztab_PSH.loc[:, 'opt_global_spectrum_reference'] = report.apply(lambda x: 'index=' + x['Lib.Index'], axis=1, result_type="expand")
    out_mztab_PSH.loc[:, 'opt_global_cv_MS:1002217_decoy_peptide'] = '0'
    out_mztab_PSH.loc[:, 'PSH'] = 'PSM'
    index = out_mztab_PSH.loc[:, 'PSH']
    out_mztab_PSH.drop(labels = ['PSH'], axis = 1, inplace = True)
    out_mztab_PSH.insert(0, 'PSH', index)
    out_mztab_PSH.loc[:, 'PSM_ID'] = out_mztab_PSH.index
    out_mztab_PSH.loc[:, 'unique'] = out_mztab_PSH.apply(lambda x: classify_protein(x['Genes'], unique["Genes"], '1', '0'), axis=1, result_type="expand")
    out_mztab_PSH.loc[:, 'database'] = database
    null_col = ['database_version', 'search_engine', 'calc_mass_to_charge', 'pre', 'post', 'start', 'end', 'opt_global_feature_id', 'opt_global_map_index']
    for i in null_col:
        out_mztab_PSH.loc[:, i] = 'null'
    out_mztab_PSH.loc[:, 'modifications'] = out_mztab_PSH.apply(lambda x: find_modification(x['opt_global_cv_MS:1000889_peptidoform_sequence']), axis = 1, result_type = 'expand')
    out_mztab_PSH.loc[:, 'opt_global_cv_MS:1000889_peptidoform_sequence'] = out_mztab_PSH.apply(lambda x: convert_modification(x["opt_global_cv_MS:1000889_peptidoform_sequence"], unimod_data), axis=1, result_type="expand")
    out_mztab_PSH = out_mztab_PSH.drop(['Genes'], axis = 1)
    out_mztab_PSH.fillna('null', inplace = True)  
    #out_mztab_PSH.to_csv('./mztab_PSH.csv', sep=',', index=False)
    print('<--      PSH has been processed !        -->')
    end = time.time()
    print("PSH run time:{t}".format(t=end - start))
    return out_mztab_PSH

def classify_protein(target, unique, t, f):
    if any(unique == target):
        return t
    else:
        return f

def find_each_params(report, target, index_ref, run, col, average):
    if average == 1:
        target_file = index_ref[index_ref['study_variable'] == run]['run'].values[0]
        result = report[report['Precursor.Id'] == target]
        match = result[result['Run'] == target_file]
        return match[col].mean() if match[col].values.size > 0 else np.nan
    elif average == 0:
        target_file = index_ref[index_ref['ms_run'] == run]['run'].values[0]
        result = report[report['Precursor.Id'] == target]
        match = result[result['Run'] == target_file]
        return match[col].values[0] if match[col].values.size > 0 else np.nan

def get_ms_run(target, index_ref):
    match = index_ref[index_ref['run'] == target]
    ms_run = match['ms_run'].values[0]
    return str(ms_run)

def convert_modification(peptide, unimod_data):
    pattern = re.compile(r"\((.*?)\)")
    origianl_mods = re.findall(pattern, peptide)
    for mod in set(origianl_mods):
        name = unimod_data[unimod_data["id"] == mod]["name"].values[0]
        peptide = peptide.replace(mod, name)
    if peptide.startswith("("):
        peptide = peptide + "."
    return peptide

def RPHfrom_report(report, target):
    match = report[report['Protein.Ids'] == target]
    modSeq = match['Modified.Sequence'].values[0] if match['Modified.Sequence'].values.size > 0 else np.nan
    score = min(match['Protein.Q.Value'].values)
    return modSeq, score

def PEHfrom_report(report, target, col):
    for i in col:
        report.loc[:, i] = report.loc[:, i].astype('float')
    match = report[report['Precursor.Id'] == target]
    search_score = min(match['Q.Value'].values) if match['Q.Value'].values.size > 0 else np.nan
    time = match['RT'].mean() if match['RT'].values.size > 0 else np.nan 
    q_score = match['Global.Q.Value'].values[0] if match['Global.Q.Value'].values.size > 0 else np.nan
    spec_e = match['Lib.Q.Value'].values[0] if match['Lib.Q.Value'].values.size > 0 else np.nan
    mz = match['Precursor.Mz'].mean() if match['Precursor.Mz'].values.size > 0 else np.nan
    return search_score, time, q_score, spec_e, mz

def find_modification(peptide):
    pattern = re.compile(r"\((.*?)\)")
    position = [j.start() for j in re.finditer('\(', peptide)]
    original_mods = re.findall(pattern, peptide)
    for i in range(0,len(original_mods)):
        original_mods[i] = str(position[i]) + '-' + original_mods[i] 
    original_mods = ','.join(str(i) for i in original_mods)
    return original_mods if len(original_mods) > 0 else 'null'

def mod_inf(modification, unimod):
    mods = []
    origianl_mods = re.split('\,|\, ', modification)
    if origianl_mods:
        for i in origianl_mods:
            mod = '[UNIMOD, ' + unimod[unimod['name'] == i.split(' ')[0]]['id'].values[0] + ', ' + i.split(' ')[0] + ', ]'
            pattern = re.compile(r"\((.*?)\)")
            site = re.findall(pattern, i)
            mods.append(mod)
            mods.append(site)
    else: return "null"
    return mods

cli.add_command(convert_report2mztab)
if __name__ == "__main__":
    cli()
