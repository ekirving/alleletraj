#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mysql.connector
import itertools

from pprint import pprint
from collections import OrderedDict


# the maximum number of rows to insert in a single operation
MAX_INSERT_SIZE = 50000


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

    def __format_data(self, params):

        data = {}

        for key, value in params.iteritems():
            new_key = u"`{}`".format(key)
            new_val = u"'{}'".format(self.cnx.converter.escape(value)) if value is not None and value != '' else 'NULL'

            data[new_key] = new_val

        return data

    def __get_records(self, table, conds=None):
        """
        Helper function for fetching records
        """
        sql = "SELECT * FROM {table} ".format(table=table)

        if conds:
            sub = [u"{}={}".format(key, value) for key, value in self.__format_data(conds).iteritems()]
            sql += u"WHERE {conds}".format(conds=u" AND ".join(sub))

        self.cursor.execute(sql)

    def get_records(self, table, conds=None, key='id'):
        """
        Get all matching records
        """
        self.__get_records(table, conds)

        return OrderedDict((item[key], item) for item in self.cursor)

    def get_records_sql(self, sql, key='id'):
        """
        Get all matching records for a given SQL query
        """
        self.cursor.execute(sql)

        if key is None:
            return list(self.cursor)
        else:
            return OrderedDict((item[key], item) for item in self.cursor)

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
        record = self.__format_data(record)

        data = {
            'table':  table,
            'fields': u", ".join(record.keys()),
            'values': u", ".join(record.values()),
            'update': u", ".join([u"{}={}".format(key, value) for key, value in record.iteritems() if key != 'id'])
        }

        sql = u"INSERT INTO {table} ({fields}) VALUES ({values}) " \
              u"ON DUPLICATE KEY UPDATE {update}".format(**data)

        try:
            self.cursor.execute(sql)
            self.cnx.commit()

        except Exception as e:
            # dump the record before throwing the exception
            print "ERROR: db_conn.save_record()"
            pprint(record)
            raise e

    def save_records(self, table, fields, records):
        """
        Batch insert new records
        """

        try:
            while True:
                # throws a StopIteration exception when we're done
                first = records.next()

                # split bulk inserts into chunks so we don't exceed the max_allowed_packet size in the DB
                data = {
                    'table': table,
                    'fields': u", ".join("`{}`".format(field) for field in fields),
                    'values': u", ".join(
                        "('" + "','".join(record) + "')" for record in
                        itertools.islice(itertools.chain([first], records), MAX_INSERT_SIZE))
                }

                sql = u"INSERT INTO {table} ({fields}) " \
                      u"VALUES ({values})".format(**data)

                self.cursor.execute(sql)
                self.cnx.commit()

        except StopIteration:
            # we're done
            pass

        except Exception as e:
            # dump the record before throwing the exception
            print "ERROR: db_conn.save_records()"
            # pprint(records)
            # pprint(sql)
            raise e
