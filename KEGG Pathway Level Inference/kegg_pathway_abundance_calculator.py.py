import pandas as pd
import os


if __name__ == '__main__':

    ko_abundance_file = 'KO_pred_metagenome_unstrat.tsv.gz'
    ko_df = pd.read_csv(ko_abundance_file, sep='\t', compression='gzip')

    mapping_ko_kegg_df = pd.read_csv('../KEGG db/KO2KEGG.csv')

    selected_col = [v for v in ko_df.columns if v != 'function']
    final_col = selected_col.copy()
    final_col.insert(0, 'Pathway')
    pathways_df = pd.DataFrame(columns=final_col)

    check_ko = []

    for kegg_pathway in mapping_ko_kegg_df['Unnamed: 0'].values:
        ko_numbers = mapping_ko_kegg_df.loc[mapping_ko_kegg_df['Unnamed: 0'] == kegg_pathway].values.tolist()[0]
        ko_numbers = [v for v in ko_numbers if str(v) != 'nan']
        check_ko.extend(ko_numbers)
        selected_co = ko_df.loc[ko_df['function'].isin(ko_numbers)]
        values_list = [kegg_pathway]
        values_list.extend(selected_co[selected_col].sum(axis=0).values.tolist())
        pathways_df.loc[pathways_df.shape[0] + 1, final_col] = values_list

    check_ko = set(check_ko)
    total_ko_functions = ko_df['function'].values.tolist()
    missed_ko_numbers = [v for v in total_ko_functions if v not in check_ko]
    if missed_ko_numbers:
        print('missed KO numbers:')
        print('#: ' + str(len(missed_ko_numbers)))
        print(missed_ko_numbers)

    output_path = os.path.dirname(ko_abundance_file)
    pathways_df.to_csv(output_path + '/KEGG_path_abun_unstrat.csv', index=False)
