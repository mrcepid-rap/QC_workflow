#!/usr/bin/env python3

import csv
import gzip

import numpy as np
import pandas as pd


def check_dupes(df: pd.DataFrame, current_code: str, age: float) -> bool:

    found_df = df[df['ICD10'].str.startswith(current_code)]
    if len(found_df) == 1:
        return False
    else:
        deduped = found_df[['eid', 'ICD10']].drop_duplicates()
        if len(deduped) == 1:  # Multiple identical codes... Keep the one with age
            if pd.isna(age):
                return True
            else:
                return False
        else:
            return True


def recover_time(df: pd.DataFrame, current_code: str) -> float:

    found_table = df[df['ICD10'] == current_code[0:3]]

    if len(found_table) == 1:
        found = search_table[search_table['ICD10'] == current_code[0:3]].iloc[0].to_dict()
        return found['age_at_incidence']
    else:
        return -1.0


search_table = pd.read_csv(gzip.open('./search_table.tsv.gz', 'rb'), sep="\t")
grouped_search_table = search_table.groupby(by='eid')

counter = 1
with open('processed_icd10.tsv', 'w') as outfile:
    output_writer = csv.DictWriter(outfile, delimiter="\t",
                                   fieldnames=['eid', 'ICD10', 'current_age', 'age_at_incidence'],
                                   extrasaction='ignore')
    output_writer.writeheader()

    for group in grouped_search_table:

        current_df = group[1]

        # Get the age of this participant
        current_age = current_df[pd.isna(current_df['current_age']) == False]['current_age'].unique()
        if len(current_age) == 0:
            current_age = np.NaN
        else:
            current_age = current_df[pd.isna(current_df['current_age']) == False]['current_age'].unique()[0]

        current_df['DROP'] = current_df.apply(lambda x: check_dupes(current_df, x['ICD10'], x['age_at_incidence']), axis=1)
        deduped_df = current_df[current_df['DROP'] == False].copy()

        # Now try and recover date of incidence...
        search_table = current_df[pd.isna(current_df['current_age']) == False]
        deduped_df['age_at_incidence'] = deduped_df['ICD10'].apply(lambda x: recover_time(search_table, x))

        # And add in age of incidence...
        deduped_df['current_age'] = current_age

        for record in deduped_df.to_dict(orient='records'):
            output_writer.writerow(record)

        # If I ever want to record progress...
        if counter % 1000 == 0:
            pass
            # print(f'{counter} individuals processed')

        counter += 1

    outfile.close()

print('processed_icd10.tsv\n')
