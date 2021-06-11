from main import connect
from Bio import Entrez
from db_connect import SQLConnection

def get_accesions(conn):
    query = 'select accesion_number from protein'
    return conn.query(query)


def get_info(accession_number):
    Entrez.email = "d.mirkarimi1@student.han.nl"

    with Entrez.efetch(db="protein", id="GFE78581", rettype="gp",
                       retmode="xml") as handle:
        record = Entrez.read(handle)[0]
    print(record)


def main():
    with open('credentials.json') as credentials_file:
        with connect(credentials_file) as han_conn:
            get_info(get_accesions(han_conn))


if __name__ == '__main__':
    main()
