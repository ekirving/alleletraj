#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mysql.connector


class db_conn:
    """
    Class for handling all the db connectivity.
    """

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

    def __get_records(self, table, conds=None):
        """
        Helper function for fetching records
        """
        sql = "SELECT * FROM {table} ".format(table=table)

        if conds:
            sub = [u"`{}`='{}'".format(key, value) for key, value in conds.iteritems()]
            sql += "WHERE {conds}".format(conds=" AND ".join(sub))

        self.cursor.execute(sql)

    def get_records(self, table, conds=None, key='id'):
        """
        Get all matching records
        """
        self.__get_records(table, conds)

        return {item[key]: item for item in self.cursor}

    def get_record(self, table, conds=None):
        """
        Get a single record
        """
        self.__get_records(table, conds)

        return self.cursor.fetchone()

    def exists_record(self, table, conds=None):
        """
        Get a single record
        """
        self.__get_records(table, conds)

        return self.cursor.fetchone() is not None

    def save_record(self, table, record):
        """
        Insert a new record
        """
        data = {
            'table':  table,
            'fields': u"`,`".join(record.keys()),
            'values': u"','".join(record.values()),
            'update': u",".join([u"`{}`='{}'".format(key, value) for key, value in record.iteritems() if key != 'id'])
        }

        sql = u"INSERT INTO {table} (`{fields}`) VALUES ('{values}') ON DUPLICATE KEY UPDATE {update}".format(**data)

        try:
            self.cursor.execute(sql)
            self.cnx.commit()

        except Exception as e:
            # dump the record before throwing the exception
            print "ERROR: db_conn.save_record()"

            for key, value in record.iteritems():
                print key, ": ", value

            raise e
