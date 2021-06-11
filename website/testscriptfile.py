
'''
def connect():
    conn = mysql.connector.connect(host='bio-inf.han.nl', user='maoqq',
                                   db='jds', password='Aa514784!',
                                   auth_plugin='mysql_native_password')

    cursor = conn.cursor()
    cursor.execute('SET autocommit = ON')

    return conn, cursor


def close(conn):
    conn.commit()
    conn.close()


def common_lineage(tax1: str, tax2: str) -> bool:
    """Checks if the taxa have one or more taxa in common in their
    lineage.

    Args:
        tax1: taxID
        tax2: taxID

    Returns:
        Bool
    """
    # Creating NCBITaxa instance containing needed methods
    ncbi = NCBITaxa()
    if ncbi._common_lineage((tax1, tax2)):
        return True
    return False
'''
from ete3 import NCBITaxa
import mysql.connector

import ete3
from db_connect import SQLConnection

ncbi = ete3.NCBITaxa()

with SQLConnection('okcum', 'Herres69!') as conn:
    taxa = conn.query('select tax_id from taxonomy')

prokaryote_count = 0
fungi_count = 0

for taxon in taxa:
    lineage = ncbi.get_lineage(taxon[0])

    if 2 in lineage or 2157 in lineage:
        prokaryote_count += 1
    elif 4751 in lineage:
        fungi_count += 1

print(fungi_count)
print(prokaryote_count)
with open('output_results.txt', 'a+') as f:
    f.write(f'Fungi count: {fungi_count}\n')
    f.write(f'Prokaryota count: {prokaryote_count}\n')