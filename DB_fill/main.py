from ete3 import NCBITaxa
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Record import Alignment
from db_connect import SQLConnection
import json
import csv
import threading
import time
from datetime import datetime
from typing import Sequence, Iterator, TextIO
from io import StringIO


def connect(credentials_file: TextIO) -> SQLConnection:
    """Returns a SQLConnection object, for interacting with the
    HAN Bio-Informatics database.

    Args:
        credentials_file: JSON file containing login information.

    Returns:
        SQLConnection object, wrapping sql.connector
    """
    creds = json.load(credentials_file)
    db_user = creds['han_user']
    db_pass = creds['han_pass']
    return SQLConnection(db_user, db_pass)


def get_lineage(organism_name: str) -> list:
    """Returns the lineage of the input organism_name.

    Args:
        organism_name: Scientific name of an organism_name.

    Returns:
        2D list containing taxID's and names of the higher taxa.
    """
    # Creating NCBITaxa instance containing needed methods
    ncbi = NCBITaxa()

    # Retrieving tax_id of organism
    try:
        tax_id = ncbi.get_name_translator(
            [organism_name])[organism_name][0]
        # Retrieving the whole lineage of the organism
        lineage = ncbi.get_lineage(tax_id)
    finally:
        ncbi.db.close()
    # Returning the lineage of the organism
    return lineage[1:]  # Excluding the first element of list 'Root'


def title_parse(title: str) -> tuple:
    """Splits the title of an alignment into the needed subunits

    Args:
        title: Hit def of an alignment.

    Returns:
        The name of the protein and organism of the hit.
    """
    protein_name, start_organism_name = title.split('[', 1)

    # Setting default value for organism_name
    organism_name = None
    nested_count = 0
    for pos, char in enumerate(start_organism_name):
        if char == '[':
            nested_count += 1
        elif char == ']':
            nested_count -= 1

        if nested_count == -1:
            organism_name = start_organism_name[:pos]
            break

    return protein_name.strip(), organism_name


def extract_results(alignment: Alignment) -> dict:
    """Extracts the blast_results of a BLAST alignment for one alignment
    object.

    Args:
        alignment: Object containing information about one alignment
        made by BLAST.

    Returns:
        Dictionary containing all relevant fragment_data from alignment.
    """

    title = alignment.hit_def
    protein_name, organism_name = title_parse(title)

    results = {
        'score': alignment.hsps[0].score,
        'expect': alignment.hsps[0].expect,
        'identities': alignment.hsps[0].identities,
        'gaps': alignment.hsps[0].gaps,
        'positives': alignment.hsps[0].positives,
        'accession': alignment.accession,
        'protein_name': protein_name,
        'organism_name': organism_name,
        'title': alignment.hit_def
    }
    return results


def update_taxonomy(current_tax: set, lineage: Sequence,
                    connection: SQLConnection) -> None:
    """Adds all inputted taxa that missing in the database, to the
    taxonomy table.

    Args:
        current_tax: A set of all taxa currently in the database.
        lineage: A sequence of taxa.
        connection: Connection to the database.

    Returns:
        None
    """
    # Creating NCBITaxa instance containing needed methods
    ncbi = NCBITaxa()
    # pos is used to track the position of the current element
    for pos, tax_id in enumerate(lineage):
        if tax_id not in current_tax:
            parent_pos = pos-1

            # Checking if a parent exists
            if parent_pos >= 0:
                parent = lineage[pos-1]
            else:
                parent = None
            name = ncbi.get_taxid_translator([tax_id])[tax_id]
            query = 'insert into taxonomy values(%s, %s, %s)'
            connection.query(query, (tax_id, name, parent))
    ncbi.db.close()


def check_lineage(organism_name: str, connection: SQLConnection) -> str:
    """Assures that all required taxa for an organism are present
    in the taxonomy table.

    Args:
        organism_name: Scientific organism_name of an organism.
        connection: Connection to the database.

    Returns:
        The taxon of the input organism.
    """
    present_taxa = connection.query('select tax_id from taxonomy')

    # Creating a set of taxa
    tax_set = {row[0] for row in present_taxa}
    lineage = get_lineage(organism_name)
    # Iterating over the lineage, starting at species
    for tax_id in reversed(lineage):
        # Checking if the taxon is not already in the database
        if tax_id not in tax_set:
            update_taxonomy(tax_set, lineage, connection)
            break
    # Returning the taxon of the organism
    return lineage[-1]


def blast_results_save(fragment_data: Sequence, blast_results: StringIO,
                       connection: SQLConnection) -> None:
    """Saves the fragments and associated results to the database.

    Args:
        fragment_data: List containing a header,
                       DNA sequence and quality scores.
        blast_results: Results of BLAST query.
        connection: Connection to the database.

    Returns:
        None
    """
    # Parsing the BLAST results
    parsed_blast_results = NCBIXML.read(blast_results)

    # Saving the fragment data to the database
    query = 'insert ignore into ' \
            'fragment(header, seq, quality) ' \
            'values(%s, %s, %s)'
    connection.query(query, fragment_data)

    # Iterating over alignments from blast
    for alignment in parsed_blast_results.alignments:
        blast_results = extract_results(alignment)
        tax_id = None

        try:
            check_lineage(blast_results['organism_name'], connection)
            tax_id = get_lineage(blast_results['organism_name'])[-1]
        except KeyError:
            print(f'An error occurred while looking up the data of '
                  f'{fragment_data[0]} as '
                  f'{blast_results["organism_name"]},'
                  f'with title {blast_results["title"]}')

        # Saving the protein data to the database
        query = 'insert ignore into protein values(%s, %s, %s, %s)'
        connection.query(query, (blast_results['accession'], tax_id,
                                 None, blast_results['protein_name']))

        # Saving the alignment data to the database
        query = 'insert into blast_result(' \
                'fragment_id, accession_number, score, expect,' \
                'identities, positives, gaps) ' \
                'values((select fragment_id  as id from fragment ' \
                'where header = %s), %s, %s, %s, %s, %s, %s)'
        connection.query(query, (fragment_data[0],
                                 blast_results['accession'],
                                 blast_results['score'],
                                 blast_results['expect'],
                                 blast_results['identities'],
                                 blast_results['positives'],
                                 blast_results['gaps']))


def processing(fragment_data: list) -> None:
    """Handles the processing of sequence fragments, which consists of
    using blast on the read and adding annotation to the blast_results.
    These will then get saved to the database.

    Args:
        fragment_data: List containing a header,
                       DNA sequence and quality scores.

    Returns:
        None
    """
    print(f'Sending blast request at {datetime.now()} for: '
          f'{fragment_data[0]}')

    # Filter to only search in sequences of prokaryotes and fungi.
    entrez = 'Procaryotae[Organism] OR Fungi[Organism]'

    result = NCBIWWW.qblast(program="blastx", database="nr",
                            sequence=fragment_data[1], ncbi_gi=True,
                            hitlist_size=5, format_type='XML',
                            matrix_name='BLOSUM62', gapcosts='11 1',
                            word_size=6, entrez_query=entrez,
                            expect=10 ^ -5)

    print(f'BLAST complete at {datetime.now()} for: {fragment_data[0]}')
    with open('credentials.json') as credentials_file:
        with connect(credentials_file) as han_conn:
            blast_results_save(fragment_data, result, han_conn)


def fragment_reader(csv_file: TextIO) -> Iterator[list]:
    """Returns the reads stored in the inputted csv file.

    Args:
        csv_file: Name of file storing fragment data.

    Returns:
        Generator that return reads in pairs.
    """
    csv_reader = csv.reader(csv_file, delimiter=';')
    for line in csv_reader:
        read1 = line[0:3]
        read2 = line[3:]
        yield read1, read2


def create_threads(allow_duplicates: bool) -> None:
    """Creates the threads for the program to run in.

    Args:
        allow_duplicates: Indicates if fragments that are already in
        the database are used or ignored.

    Returns:
        None
    """

    if allow_duplicates:
        with open('dataset_small.csv') as csv_file:
            for reads in fragment_reader(csv_file):
                for read in reads:
                    threading.Thread(target=processing,
                                     args=[read]).start()
                    time.sleep(0.6)

    else:
        with open('credentials.json') as credentials_file:
            with connect(credentials_file) as han_conn:
                query = 'select header from fragment'
                present_taxa = han_conn.query(query)
                fragment_set = {row[0] for row in present_taxa}

        with open('dataset_small.csv') as csv_file:
            for reads in fragment_reader(csv_file):
                for read in reads:
                    if read[0] in fragment_set:
                        continue
                    threading.Thread(target=processing,
                                     args=[read]).start()
                    time.sleep(0.6)


if __name__ == '__main__':
    create_threads(True)
