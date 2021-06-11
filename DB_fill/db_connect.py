from typing import Sequence
import mysql.connector


class SQLConnection:
    def __init__(self, username: str, password: str = None):
        self.params = {'host': 'bio-inf.han.nl', 'user': username,
                       'db': 'jds', 'password': password,
                       'auth_plugin': 'mysql_native_password'}
        self.conn = mysql.connector.connect(**self.params)
        self.cursor = self.conn.cursor()

        self.cursor.execute('SET autocommit = ON')

    def get_cursor(self):
        return self.cursor

    def query(self, query: str, args: Sequence = None) -> tuple:
        if args:
            self.cursor.execute(query, args)
        else:
            self.cursor.execute(query)

        try:
            return self.cursor.fetchall()
        except mysql.connector.errors.InterfaceError:
            return (None,)

    def close(self):
        self.conn.commit()
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
