import sys

import pandas as pd

from bioscreen._imports import *
import requests
from pathlib import Path

# when done use the stable address given here
#   https://string-db.org/api/json/version
APIURL = 'https://string-db.org/api/'
APPNAME = 'JTs STRINGdb interface'

# note: identifiers are 9606.ENSP0000###
def get_string_ids(genes:Collection[str], species='9606') \
        -> pd.DataFrame:

    queryurl = APIURL+'json/get_string_ids'
    params = dict(identifiers='\r'.join(genes),
                  species=species,
                  echo_query=1,
                  limit=1)
    res = requests.post(queryurl, data=params)
    if res.ok:
        return pd.DataFrame(res.json())
    else:
        logging.warning(f'Query failed with reason {res.reason}')


# def get_network_table(
#         stringids:Collection[str],
#         species='9606',
#         required_score=0.4,
#         caller_identity=APPNAME,
# ) -> pd.DataFrame:
#     """
#
#     See: https://string-db.org/help/api/#getting-the-string-network-interactions"""
#
#     #[output-format]/network?identifiers=[your_identifiers]&[optional_parameters]
#     queryurl = APIURL+'network'
#
#     params = dict(
#         identifiers='\r'.join(stringids),
#         species=species,
#
#     )



import gzip

def pair_score_no_textmining(fullscores_fn, out_picklefn, min_score=99):
    # assume it's a test string
    if type(fullscores_fn) is list:
        opener = lambda x: iter(x)
    elif fullscores_fn.endswith('.gz'):
        opener = lambda fn: gzip.open(fn, mode='rt')
    else:
        opener = open
    prior = 0.041

    def compute_prior_away(score, prior):

        if score < prior: score = prior
        score_no_prior = (score - prior) / (1 - prior)

        return score_no_prior

    scores_no_txt = {}
    failures = 0
    f = opener(fullscores_fn, )
    next(f)
    for line in f:
        line = line
        l = line.split()
        ## load the line

        try:
            (protein1, protein2,
             neighborhood, neighborhood_transferred,
             fusion, cooccurrence,
             homology,
             coexpression, coexpression_transferred,
             experiments, experiments_transferred,
             database, database_transferred,
             textmining, textmining_transferred,
             initial_combined) = l
        except:
            failures += 1
            continue

        ## divide by 1000

        neighborhood = float(neighborhood) / 1000
        neighborhood_transferred = float(neighborhood_transferred) / 1000
        fusion = float(fusion) / 1000
        cooccurrence = float(cooccurrence) / 1000
        #homology = float(homology) / 1000 # Missing on purpose
        coexpression = float(coexpression) / 1000
        coexpression_transferred = float(coexpression_transferred) / 1000
        experiments = float(experiments) / 1000
        experiments_transferred = float(experiments_transferred) / 1000
        database = float(database) / 1000
        database_transferred = float(database_transferred) / 1000
        #textmining = float(textmining) / 1000
        #textmining_transferred = float(textmining_transferred) / 1000


        ## compute prior away

        neighborhood_prior_corrected = compute_prior_away(neighborhood, prior)
        neighborhood_transferred_prior_corrected = compute_prior_away(neighborhood_transferred, prior)
        fusion_prior_corrected = compute_prior_away(fusion, prior)
        cooccurrence_prior_corrected = compute_prior_away(cooccurrence, prior)
        coexpression_prior_corrected = compute_prior_away(coexpression, prior)
        coexpression_transferred_prior_corrected = compute_prior_away(coexpression_transferred, prior)
        experiments_prior_corrected = compute_prior_away(experiments, prior)
        experiments_transferred_prior_corrected = compute_prior_away(experiments_transferred, prior)
        database_prior_corrected = compute_prior_away(database, prior)
        database_transferred_prior_corrected = compute_prior_away(database_transferred, prior)
        #textmining_prior_corrected = compute_prior_away(textmining, prior)
        #textmining_transferred_prior_corrected = compute_prior_away(textmining_transferred, prior)

        ## then, combine the direct and transferred scores for each category:

        neighborhood_both_prior_corrected = 1.0 - (1.0 - neighborhood_prior_corrected) * (
                    1.0 - neighborhood_transferred_prior_corrected)
        coexpression_both_prior_corrected = 1.0 - (1.0 - coexpression_prior_corrected) * (
                    1.0 - coexpression_transferred_prior_corrected)
        experiments_both_prior_corrected = 1.0 - (1.0 - experiments_prior_corrected) * (
                    1.0 - experiments_transferred_prior_corrected)
        database_both_prior_corrected = 1.0 - (1.0 - database_prior_corrected) * (
                    1.0 - database_transferred_prior_corrected)
        #textmining_both_prior_corrected = 1.0 - (1.0 - textmining_prior_corrected) * (
         #           1.0 - textmining_transferred_prior_corrected)

        ## next, do the 1 - multiplication:

        combined_score_one_minus = (
                (1.0 - neighborhood_both_prior_corrected) *
                (1.0 - fusion_prior_corrected) *
                (1.0 - cooccurrence_prior_corrected) *
                (1.0 - coexpression_both_prior_corrected) *
                (1.0 - experiments_both_prior_corrected) *
                (1.0 - database_both_prior_corrected)
                #(1.0 - textmining_both_prior_corrected) *
        )

        ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*

        combined_score = (1.0 - combined_score_one_minus)  ## 1- conversion
        combined_score *= (1.0 - prior)  ## scale down
        combined_score += prior  ## and add prior.

        ## round

        combined_score = int(combined_score * 1000)
        if combined_score > min_score:
            scores_no_txt[(protein1, protein2)] = combined_score
    print(failures)
    import pickle
    with open(out_picklefn, 'wb') as f:
        pickle.dump(scores_no_txt, f)

def _test():
    tabl = '''protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
    9606.ENSP00000000233 9606.ENSP00000356607 0 0 0 0 0 0 45 0 134 0 0 0 81 173
    9606.ENSP00000000233 9606.ENSP00000427567 0 0 0 0 0 0 0 0 128 0 0 0 70 154'''.split('\n')

    pair_score_no_textmining('/mnt/m/tmp/string.txt', '/mnt/m/tmp/fff.pickle')



testset = """9606.ENSP00000265986
9606.ENSP00000265038
9606.ENSP00000329797
9606.ENSP00000417610
9606.ENSP00000351035
9606.ENSP00000345344
9606.ENSP00000263934
9606.ENSP00000376475
9606.ENSP00000393566
9606.ENSP00000273261
9606.ENSP00000248437
9606.ENSP00000333926
9606.ENSP00000478887
9606.ENSP00000425824
9606.ENSP00000359520
9606.ENSP00000358799
9606.ENSP00000377262
9606.ENSP00000215071
9606.ENSP00000201031
9606.ENSP00000402060""".split('\n')

if __name__ == '__main__':
    if sys.argv[1] == 'scorenotext':
        pair_score_no_textmining(sys.argv[2], sys.argv[3])