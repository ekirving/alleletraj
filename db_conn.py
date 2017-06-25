#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mysql.connector


class db_conn:

    db_config = {
      'user':     'root',
      'password': '',
      'host':     '127.0.0.1',
      'database': 'allele_trajectory'
    }

    def __init__(self):
        # connect to the db
        self.cnx = mysql.connector.connect(**self.db_config)
        self.cursor = self.cnx.cursor(dictionary=True)

    def __del__(self):
        # close the connection
        self.cursor.close()
        self.cnx.close()

    def __get_records(self, table, conds = None):
        """
        Helper function for fetching records
        """
        sql = "SELECT * FROM {table} ".format(table=table)

        if conds:
            sql += "WHERE {conds}".format(conds="AND".join(conds))

        self.cursor.execute(sql)

    def get_records(self, table, conds=None):
        """
        Get all matching records
        """
        self.__get_records(table, conds)

        return {item['id']: item for item in self.cursor}

    def get_record(self, table, conds = None):
        """
        Get a single record
        """
        self.__get_records(table, conds)

        return self.cursor.fetchone()

    def save_record(self, table, record):
        """
        Insert a new record
        """
        sql = "INSERT INTO {table} (`{keys}`) VALUES ('{values}') ".format(
            table=table, keys='`,`'.join(record.keys()), values="','".join(record.values()))

        self.cursor.execute(sql)
        self.cnx.commit()

