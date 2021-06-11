import numpy as np
from db_connect import SQLConnection

with SQLConnection('okcum', 'Herres69!') as conn:
    query_results = conn.query('select fragment_id, quality_reads '
                               'from fragment')

    for row in query_results:
        fragment_id, read_scores = row
        score_list = []
        for score in read_scores:
            adjusted_score = ord(score)
            score_list.append(adjusted_score)
        fragment_score = np.median(score_list)
        print(fragment_score)
        conn.query('update fragment set quality_fragment=%s where '
                   'fragment_id = %s', (fragment_score,
                                        fragment_id))
