import pandas as pd
from bioservices import KEGG


if __name__ == '__main__':

    df = pd.read_csv('KO_pred_metagenome_unstrat.tsv.gz',
                     sep='\t', 
                     compression='gzip')

    # Initialize KEGG service
    kegg = KEGG()

    pathways_dict = {}
    cnt = 0
    # Fetch pathways for each KO
    for ko in df['function'].values:
        cnt += 1
        if cnt % 1000 == 0:
            print(cnt)
        pathways = kegg.link("pathway", ko)
        try:
            consider_pathways = pathways.split('\n')
            selected_pathways = [v.split('path:')[1] for v in consider_pathways if v.count('map') == 0 and v != '']

            for path in selected_pathways:
                if path in pathways_dict.keys():
                    pathways_dict[path].append(ko)
                else:
                    pathways_dict[path] = [ko]
        except:
            print(ko + ' not found')

    df = pd.DataFrame.from_dict(pathways_dict, orient='index')
    df.to_csv('KO2KEGG.csv')

