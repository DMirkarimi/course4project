from typing import Tuple, List
import matplotlib.pyplot as plt
import mysql.connector
from mysql.connector import MySQLConnection
from mysql.connector.cursor import MySQLCursor
from ete3 import NCBITaxa
from flask import Flask, render_template, request

QueryResults = Tuple[Tuple[str, ...], ...]
Header = Tuple[str, ...]
app = Flask(__name__)


def connect() -> Tuple[MySQLConnection, MySQLCursor]:
    """connecteert aan de databse

    Returns: conn-object-een connectie aan de de database
             cursor-object-een cursor te zoeken in de database

    """
    conn = mysql.connector.connect(host='bio-inf.han.nl', user='maoqq',
                                   db='jds', password='Aa514784!',
                                   auth_plugin='mysql_native_password')

    cursor = conn.cursor()
    cursor.execute('SET autocommit = ON')

    return conn, cursor


def close(conn: MySQLConnection) -> None:
    """
    closes the connection with the database when not needed anymore
    Args:
        conn: conn-object-a connection to the database
    """
    conn.commit()
    conn.close()


def get_results(cursor: MySQLCursor, query: str,
                search_word: str = None) -> QueryResults:
    """runs the given query with the given searchword and returns the
    results
    Args:
        cursor-object-an object to search with in the databse
        query-str-str of query needed
        search_word-str- str of the search word

    Returns:
        query_results-tuple- tuple containing the results of the query

    """
    if search_word:
        cursor.execute(query, (search_word,))
    else:
        cursor.execute(query)
    query_results = cursor.fetchall()
    return query_results


def search(cursor: MySQLCursor, search_word: str) -> \
        Tuple[QueryResults, Header]:
    """determines what the type of search is by search_type and than
    calls the correct function according to the search type
    Args:
        cursor-object-an object to search with in the databse
        search_word-str- str of the search word

    Returns:
        results-tuple- tuple containing the results of the query
        headers-list-list of headers for the tables on the site
    """
    search_type = request.values.get('search_type')
    search_query = ''
    header = ''
    if search_type == 'protein_name':
        search_word = '%' + search_word + '%'
        search_query, header = search_protein()
    elif search_type == 'read_id':
        search_word = search_word + '%'
        search_query, header = search_fragment()
    elif search_type == 'tax_id':
        search_query, header = search_tax(search_word)
        search_word = None
    elif search_type == 'accession_code':
        search_query, header = search_accession()
    results = get_results(cursor, search_query, search_word)
    return results, header


def search_protein() -> Tuple[str, Header]:
    """
    returns the correct search query for the type of search
    Returns:
        search_query-str-str of query needed
        headers-list-list of headers for the tables on the site

    """
    search_query = 'select * from blast_result ' \
            'join protein using(accession_number)' \
            'where protein_name like %s;'
    header = ('Accession number', 'blast id', 'fragment id', 'score',
              'e-value', 'identity', 'positives', 'gaps', 'tax_id',
              'protein name')
    return search_query, header


def search_fragment() -> Tuple[str, Header]:
    """ returns the correct search query for the type of search
    Returns:
        search_query-str-str of query needed
        headers-list-list of headers for the tables on the site
    """
    search_query = 'select header, seq, quality_fragment, ' \
                   'accession_number, score, expect, ' \
                   'identities, positives, gaps ' \
                   'from fragment f join blast_result br' \
                   ' on f.fragment_id = br.fragment_id' \
                   ' where Right(header, 6) like %s;'

    header = ('fragment id', 'sequence', 'quality score',
              'accession number', 'score', 'e-value', 'identity',
              'positives', 'gaps')
    return search_query, header


def search_tax(search_word: str) -> Tuple[str, Header]:
    """
    returns the correct search query for the type of search
    Args:
        search_word-str- str of the search word

    Returns:
        search_query-str-str of query needed
        headers-list-list of headers for the tables on the site
    """
    taxa_list = NCBITaxa().get_descendant_taxa(
        search_word, intermediate_nodes=True)

    taxa_string = f'({",".join([str(taxon) for taxon in taxa_list])})'

    search_query = 'select accession_number, ' \
                   'protein_name, p.tax_id, tax_name ' \
                   'from protein  p natural join taxonomy t' \
                   ' where p.tax_id in {};'.format(taxa_string)
    header = ('accession number', 'protein name',
              'tax id', 'tax name')
    return search_query, header


def search_accession() -> Tuple[str, Header]:
    """
      returns the correct search query for the type of search
    Returns:
        search_query-str-str of query needed
        headers-list-list of headers for the tables on the site
    """
    search_query = 'select * from blast_result b ' \
            'join protein p using(accession_number)' \
            'where b.accession_number like %s;'

    header = ('Accession number', 'blast id', 'fragment id', 'score',
              'e-value', 'identity', 'positives', 'gaps',
              'tax_id', 'protein name')
    return search_query, header


def get_common_lineage(tax1: str, tax2: str) -> bool:
    """Checks if the taxa have one or more taxa in common in their
    lineage.
    Args:
        tax1: search term taxID
        tax2: taxID to be checked
    Returns:
        Bool
    """
    # Creating NCBITaxa instance containing needed methods
    ncbi = NCBITaxa()
    common_lineage = ncbi._common_lineage((tax1, tax2))
    if common_lineage:
        return True
    return False


def search_common(results: QueryResults,
                  search_word: str) -> QueryResults:
    """looks through
    Args:
        results:
        search_word:

    Returns:

    """
    filtered_results = []
    for row in results:
        if get_common_lineage(row[1], search_word):
            filtered_results.append(row)
    return tuple(filtered_results)


def default_search(cursor: MySQLCursor) -> Tuple[
                                        QueryResults,
                                        Header]:
    """passes the query for the default page to the get_resulsts
    function and parses it.
    Args:
        cursor-object-an object to search with in the databse

    Returns:
        search_query-str-str of query needed
        headers-list-list of headers for the tables on the site

    """

    search_query = 'select * from blast_result ' \
                'join protein using(accession_number); '
    cursor.execute(search_query)
    query_results = cursor.fetchall()
    header = ('Accession number', 'blast id', 'fragment id', 'score',
              'e-value', 'identity', 'positives', 'gaps', 'tax_id',
              'protein name')
    return query_results, header


def search_info(cursor: MySQLCursor, search_word_tax: str) -> None:
    """runs query with search_word and parses the results to other
    functions.
    Args:
        cursor-object-an object to search with in the databse
        search_word_tax-str or int- searcht term for query
    """
    query = 'select tax_name, count(p.tax_id)  counts from protein p ' \
            'join taxonomy t on p.tax_id = t.tax_id ' \
            'where tax_name like %s ' \
            'group by tax_name order by counts desc;'
    query_total = 'select count(*) from protein;'
    search_word = '%' + search_word_tax + '%'
    my_results = get_results(cursor, query, search_word)
    total = get_results(cursor, query_total)[0][0]
    labels, counts = parse_list(my_results, total)
    pie_chart(labels, counts)


def parse_list(results, total) -> Tuple[List[str],
                                        List[str]]:
    """parses the results to make two different lists. A value list
    and a labels list.
    Args:
        results-tuple- tuple containing the results of the query
        total-tuple- tuple of number of total found

    Returns:
        labels-list-list of labels of counted product
        counts-list-list of counts of counted product
    """
    labels = []
    counts = []
    for item in results:
        labels.append(item[0])
        count_percentage = str(round(item[1]/total*100, 1))
        counts.append(count_percentage)
    return labels, counts


def pie_chart(labels: List[str], counts: List[str]) -> None:
    """Takes a list of values and a list of labels and makes a
    pie chart, then saves the figure.
    Args:
        labels-list-list of labels of counted product
        counts-list-list of counts of counted product

    Returns:

    """
    plt.close()
    plt.pie(counts, labels=counts, autopct='%1.2f%%')
    plt.legend(labels, loc='upper center', bbox_to_anchor=(0.5, 0.05),
          ncol=3, fontsize='x-small')
    plt.suptitle('percentage of count of results', fontsize=14)
    plt.title('inner ring is percentage of search, outer ring is percentage of all db',
                 y=1, fontsize=8)
    plt.savefig('static/images/plaatje.png')


@app.route('/', methods=['GET', 'POST'])
def homepage() -> str:
    return render_template('home.html')


@app.route('/input', methods=['GET', 'POST'])
def inputpage() -> str:
    search_word = request.values.get('zoekwoord')
    conn, cursor = connect()
    if search_word:
        results, header = search(cursor, search_word)
        close(conn)
    else:
        results, header = default_search(cursor)
        close(conn)
    return render_template('input.html', info=results, header=header)


@app.route('/infopage', methods=['GET', 'POST'])
def info() -> str:
    search_word_tax = request.values.get('tax_zoek_woord')
    if search_word_tax:
        conn, cursor = connect()
        search_info(cursor, search_word_tax)
        return render_template('infopage.html',
                               img_url='/static/images/plaatje.png')
    return render_template(
        'infopage.html', img_url='/static/images/insert_image_here.png')


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
