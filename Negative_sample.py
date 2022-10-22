import EntryModule as em
import peptide_assembly as pept

def create_negative_sample_tables(path_to_table_ABL, explorable_drug):
    all_drugs = em.ABL_table_analyzer.get_unique_values_in_column(path_to_table_ABL, ' DRUG_NAME')
