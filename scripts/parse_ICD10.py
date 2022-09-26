#!/usr/bin/env python3

import re
import sys
import csv
import math
import argparse
import datetime
from dataclasses import dataclass
from pathlib import Path

from dateutil.relativedelta import relativedelta
from typing import Dict, TypedDict, List, Tuple

import numpy as np
import pandas as pd


# This is a hard-coded max date that should only be changed on death database/first occurrence
# being updated (whatever is earlier)
MAX_DATE = datetime.date.fromisoformat('2021-11-12')


def file_type(file_string: str) -> Path:

    current_file = Path(file_string)
    if not current_file.exists():
        print(f'The provided file - {file_string} – does not appear to exist.\n\nExiting...\n')
        sys.exit(1)

    return current_file


@dataclass
class Args:

    eid: str
    codes: str
    tree: Path
    participants: Path
    conditions: Path
    set_na: bool
    is_tte: bool
    output: Path


def parse_args() -> Args:

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--eid', help="participant eid to search codes for.", type=str, dest='eid', required=False)
    parser.add_argument('-d', '--codes', help="ICD-10 code(s) to extract phenotype(s) for.", type=str, dest='codes', required=False, nargs='*')
    parser.add_argument('-i', '--icd_tree', help="Path to the ICD-10 tree.", type=file_type, dest='tree', required=True)
    parser.add_argument('-p', '--participants', help="Path to vital statistics for study participants to include.", type=file_type, dest='participants', required=True)
    parser.add_argument('-c', '--conditions', help="Path to processed participant conditions.", type=file_type, dest='conditions', required=True)
    parser.add_argument('-n', '--na', help="Set individuals with codes higher on the ICD10 tree to 'NA' [False]", action='store_true', dest='set_na')
    parser.add_argument('-t', '--tte', help="Generate a time to event phenotype for the requested code [False]", action='store_true', dest='is_tte')
    parser.add_argument('-o', '--output', help="Name of the output file", type=Path, dest='output', default=None, required=False)

    opts = Args(**vars(parser.parse_args()))

    if opts.eid is None and opts.codes is None:
        print('At least eid (-e/--eid) or code (-c/--code) must be provided.\n\nExiting...\n')
        sys.exit(1)
    elif opts.eid is not None and opts.codes is not None:
        print('Both eid and code cannot be provided.\n\nExiting...\n')
        sys.exit(1)

    if opts.is_tte and opts.set_na is False:
        print('WARNING: Building a time-to-event phenotype without setting individuals with ICD-10 conditions higher in'
              ' the ICD-10 hierarchy to "NA". This may result in erroneous phenotype definitions!')

    return opts


class ICDInfo(TypedDict):
    description: dict
    long_code: str
    parent: int
    children: List[str]


def load_icd_tree(icd_tree_file: Path) -> Dict[str, ICDInfo]:

    icd_tree = {}

    with icd_tree_file.open('r') as icd_tree_reader:
        icd_reader = csv.DictReader(icd_tree_reader, delimiter="\t")
        for code in icd_reader:
            chapter_match = re.match('(Chapter [XVI]*)\\s([\\S\\s]*)', code['meaning'])
            name_match = re.match('(\\S+)\\s{1}([\\S\\s]*)', code['meaning'])

            long_code = name_match.group(1) if chapter_match is None else chapter_match.group(1)
            description = name_match.group(2) if chapter_match is None else chapter_match.group(2)
            icd_tree[int(code['node_id'])] = {'short_code': long_code if 'Block' in code['coding'] else code['coding'],
                                              'long_code': long_code,
                                              'description': description,
                                              'parent': int(code['parent_id']),
                                              'children': []}
        icd_tree_reader.close()

    for key, value in icd_tree.items():
        if value['parent'] != 0:
            icd_tree[value['parent']]['children'].append(value['short_code'])
            icd_tree[key]['parent'] = icd_tree[value['parent']]['short_code']
        else:
            icd_tree[key]['parent'] = None

    final_tree = {}
    for key, value in icd_tree.items():
        code = value['short_code']
        del value['short_code']
        final_tree[code] = value

    return final_tree


def load_individual_codes(individual_code_file: Path) -> pd.DataFrame:

    codes = pd.read_csv(individual_code_file, sep="\t", dtype={'eid': np.character})
    codes = codes.set_index('eid')

    return codes


class ParticipantInfo(TypedDict):
    date_of_death: datetime.date
    date_of_birth: datetime.date
    is_deceased: bool
    age: int


def load_participants(participants_file: Path) -> Dict[str, ParticipantInfo]:

    participants_dict = {}

    with participants_file.open('r') as participants_reader:
        participants_csv = csv.DictReader(participants_reader, delimiter="\t")
        for participant in participants_csv:

            dob = datetime.date.fromisoformat(participant['dob'])
            if participant['dod'] == 'NA':
                dod = None
                age = relativedelta(MAX_DATE, dob).years  # This is a hardcoded max date!
                is_deceased = False
            else:
                dod = datetime.date.fromisoformat(participant['dod'])
                age = relativedelta(dod, dob).years
                is_deceased = True

            participants_dict[participant['eid']] = {'date_of_death': dod,
                                                     'date_of_birth': dob,
                                                     'age': age,
                                                     'is_deceased': is_deceased}

        return participants_dict


# This is a simple query to find all ICD-10 codes for a given individual and give long names to them
def get_individual_codes(eid: str, codes: pd.DataFrame, tree: Dict[str, ICDInfo], participants: Dict[str, ParticipantInfo]):

    vital_stats = participants[eid]
    print(f'Participant id: {eid}')
    print(f'Participant DoB: {vital_stats["date_of_birth"]}')
    print(f'Participant DoD: {vital_stats["date_of_death"]}')
    print(f'Participant is: {"Alive" if vital_stats["is_deceased"] is False else "Deceased"}')
    print(f'Participant age: {vital_stats["age"]}\n')
    print(f'Found Codes:')

    found_codes = codes.loc[eid]

    for row in found_codes.iterrows():
        current_code = row[1]['ICD10']
        age_at_incidence = float('NaN') if row[1]['age_at_incidence'] == -1.0 else row[1]['age_at_incidence']
        code_meaning = tree[current_code]['description']
        print(f'ICD10 Code: {current_code}')
        print(f'\tMeaning: {code_meaning}')
        print(f'\tAge at incidence: {age_at_incidence: 0.2f}')


def find_children(current_code: str, tree: Dict[str, ICDInfo], relevant_codes=None) -> List[str]:

    if relevant_codes is None:
        relevant_codes = [current_code]
    else:
        relevant_codes.append(current_code)

    for child in tree[current_code]['children']:
        relevant_codes = find_children(current_code=child, tree=tree, relevant_codes=relevant_codes)

    return relevant_codes


def find_parents(current_code: str, tree: Dict[str, ICDInfo], relevant_codes=None) -> List[str]:

    parent = tree[current_code]['parent']
    if relevant_codes is None:
        relevant_codes = [parent]
    elif parent is not None:
        relevant_codes.append(parent)

    if parent is not None:
        relevant_codes = find_parents(parent, tree, relevant_codes)

    return relevant_codes


def get_individuals(code: str, codes: pd.DataFrame, tree: Dict[str, ICDInfo]) -> Tuple[pd.DataFrame, pd.DataFrame]:

    if code in tree:

        # First walk down and up the ICD tree (in that order) to find codes relevant to the code of interest
        relevant_codes = find_children(code, tree)
        parents = find_parents(code, tree)

        # Query the listing of codes for:
        # 1. Queried code and all codes 'below' it in on the ICD-10 tree
        # 2. All 'parents' above the queried code in case we want to set to NA individuals
        found_individuals = codes[codes['ICD10'].isin(relevant_codes)]
        found_parents = codes[codes['ICD10'].isin(parents)]

        return found_individuals, found_parents

    else:
        print(f'Code {code} not found in the ICD-10 tree. Please try again...\nExiting\n\n')
        sys.exit(1)


def add_information(found_individuals: pd.DataFrame, found_parents: pd.DataFrame, participant: str,
                    set_na: bool, is_tte: bool):

    # First just set the code regardless of the 'set_na' flag
    if participant in found_individuals.index:
        if is_tte:
            # Define rules to get age_at_incidence if searching a code that can return multiple results (e.g. Blocks)
            time_info = found_individuals.loc[participant]
            if type(time_info) == pd.Series:  # Will always be a series if only one code found
                age_at_incidence = time_info['age_at_incidence']
            else:
                time_info = time_info[time_info['age_at_incidence'] >= 0.0]
                if len(time_info) == 0:
                    age_at_incidence = -1.0
                else:
                    age_at_incidence = time_info['age_at_incidence'].min()

            if set_na and participant in found_parents.index:
                return pd.NA
            else:
                if age_at_incidence == -1.0:
                    return pd.NA
                else:
                    return math.floor(age_at_incidence)
        else:
            if set_na and participant in found_parents.index:
                return pd.NA
            else:
                return 1
    else:
        if is_tte:
            if set_na and participant in found_parents.index:
                return pd.NA
            else:
                return 0
        else:
            if set_na and participant in found_parents.index:
                return pd.NA
            else:
                return 0


# First load required data
args = parse_args()
icd_tree = load_icd_tree(args.tree)
participants = load_participants(args.participants)
individual_codes = load_individual_codes(args.conditions)

# icd_tree = load_icd_tree('../coding19.tsv')
# participants = load_participants('../vital_stats.tsv')
# individual_codes = load_individual_codes('../processed_icd10.tsv')

# Then do requested processing. This is safe without an else because we check input parameters at launch, above
if args.eid:
    get_individual_codes(eid=args.eid, codes=individual_codes, tree=icd_tree, participants=participants)
elif args.codes:
    # Set up a pandas DataFrame to hold necessary information:
    participant_frame = pd.DataFrame.from_dict(data=participants, orient='index')
    participant_frame['IID'] = participant_frame.index
    participant_frame['FID'] = participant_frame.index
    participant_frame = participant_frame.rename(columns={'age': 'exposure_years'})
    if args.is_tte:
        participant_frame = participant_frame[['IID', 'FID', 'exposure_years']]
    else:
        participant_frame = participant_frame[['IID', 'FID']]

    # And iterate through all codes given to the 'codes' input param
    for code in args.codes:
        individuals_with_code, individuals_with_parent = get_individuals(code=code,
                                                                         codes=individual_codes,
                                                                         tree=icd_tree)

        participant_frame[code] = participant_frame.apply(lambda x: add_information(
            found_individuals=individuals_with_code,
            found_parents=individuals_with_parent,
            participant=x['IID'],
            set_na=args.set_na,
            is_tte=args.is_tte), axis=1)

        code_counts = participant_frame[code].value_counts(dropna=False)
        total_individuals = code_counts.sum()
        non_zero_counts = code_counts.iloc[code_counts.index != 0].sum()
        if pd.NA in code_counts.index:
            na_counts = code_counts.loc[pd.NA]
        else:
            na_counts = 0

        print(f'{"Total individuals parsed for code " + code:{40}}: {total_individuals}')
        print(f'{"Total individuals with non-zero codes":{40}}: {non_zero_counts}')
        print(f'{"Total individuals with NA codes":{40}}: {na_counts} ({(na_counts / non_zero_counts)*100:0.2f}%)')

    # Write output file
    if args.output is not None:
        pheno_file = args.output
    elif len(args.codes) == 1:
        pheno_file = Path(f'{args.codes[0]}.pheno')
    else:
        pheno_file = Path(f'phenotype.pheno')
    participant_frame.to_csv(pheno_file, sep="\t", na_rep='NA',index=False)
