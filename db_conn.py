#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

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

    def __format_conditions(self, conds):

        sub = [u"{}={}".format(key, value) for key, value in self.__format_data(conds).iteritems()]

        return u"WHERE {conds}".format(conds=u" AND ".join(sub))

    def __get_records(self, table, conds=None, sort=None):
        """
        Helper function for fetching records
        """
        sql = "SELECT * FROM {table} ".format(table=table)

        if conds:
            sql += self.__format_conditions(conds)

        if sort:
            sql += "ORDER BY {sort} ".format(sort=sort)

        self.cursor.execute(sql)

    def __delete_records(self, table, conds=None):
        """
        Helper function for deleting records
        """
        sql = "DELETE FROM {table} ".format(table=table)

        if conds:
            sql += self.__format_conditions(conds)

        return self.cursor.execute(sql)

    def __count_records(self, table, conds=None):
        """
        Helper function for counting records
        """
        sql = "SELECT COUNT(*) FROM {table} ".format(table=table)

        if conds:
            sql += self.__format_conditions(conds)

        self.cursor.execute(sql)

    def get_records(self, table, conds=None, sort=None, key='id'):
        """
        Get all matching records
        """
        self.__get_records(table, conds, sort)

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

    def delete_records(self, table, conds=None):
        """
        Delete all matching records
        """
        return self.__delete_records(table, conds)

    def count_records(self, table, conds=None):
        """
        Count all matching records
        """
        self.__count_records(table, conds)

        return self.cursor.fetchone()[u'COUNT(*)']

    def save_record(self, table, record, insert=None):
        """
        Insert/update a record
        """
        formatted = self.__format_data(record)

        data = {
            'table': u"`{}`".format(table),
            'fields': u", ".join(formatted.keys()),
            'values': u", ".join(formatted.values()),
            'update': u", ".join([u"{}={}".format(key, value) for key, value in formatted.iteritems() if key != u'`id`'])
        }

        if 'id' in record and not insert:
            data['id'] = record['id']

            # update an existing record
            sql = u"UPDATE {table} SET {update} " \
                  u"WHERE `id` = {id}".format(**data)
        else:
            # insert new record
            sql = u"INSERT INTO {table} ({fields}) VALUES ({values}) " \
                  u"ON DUPLICATE KEY UPDATE {update}".format(**data)

        try:
            self.cursor.execute(sql)
            self.cnx.commit()

            return self.cursor.lastrowid

        except Exception as e:
            # dump the record before throwing the exception
            print("ERROR: db_conn.save_record()")
            pprint(sql)
            raise e

    def save_records(self, table, fields, records):
        """
        Batch insert new records
        """

        if type(records) is list:
            # convert to an iterable
            records = itertools.chain.from_iterable([records])

        try:
            while True:
                # throws a StopIteration exception when we're done
                first = records.next()

                # split bulk inserts into chunks so we don't exceed the max_allowed_packet size in the DB
                data = {
                    'table': table,
                    'fields': u", ".join("`{}`".format(field) for field in fields),
                    'values': u", ".join(
                        "('" + "','".join(map(str, record)) + "')" for record in
                        itertools.islice(itertools.chain([first], records), MAX_INSERT_SIZE))
                }

                sql = u"INSERT INTO {table} ({fields}) " \
                      u"VALUES {values}".format(**data)

                self.cursor.execute(sql)
                self.cnx.commit()

        except StopIteration:
            # we're done
            pass

        except Exception as e:
            # dump the record before throwing the exception
            print("ERROR: db_conn.save_records()")
            pprint(list(records))
            # pprint(sql)
            raise e

    def execute_sql(self, sql):
        """
        Execute a given SQL query
        """
        self.cursor.execute(sql)
        self.cnx.commit()