#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import itertools
import mysql.connector

from collections import OrderedDict
from datetime import timedelta
from pprint import pprint
from time import time

# import my custom modules
from pipeline_consts import REF_ASSEMBLY
from pipeline_qtls import QTLDB_RELEASE


class DBConn:
    """
    Class for handling all the db connectivity.
    """

    db_config = {
      'user':     'root',
      'password': '',
      'host':     '127.0.0.1'
    }

    # the maximum number of rows to insert in a single operation
    max_insert_size = 50000

    # the maximum number of conditions in a single query
    max_query_size = 5000

    def __init__(self, species):
        # set the database name
        self.db_config['database'] = DBConn.__get_name(species)

        # connect to the db
        self.cnx = mysql.connector.connect(**self.db_config)
        self.cursor = self.cnx.cursor(dictionary=True)

    def __del__(self):
        # close the connection
        self.cursor.close()
        self.cnx.close()

    @staticmethod
    def __get_name(species):
        """
        Embed the reference assembly and the QTLdb release number into the database name.
        """
        return 'alleletraj_{}_{}_{}'.format(species, REF_ASSEMBLY[species], QTLDB_RELEASE).lower()

    def __format_data(self, params):

        data = {}

        for key in params:
            new_key = u"`{}`".format(key)
            new_val = u"'{}'".format(self.cnx.converter.escape(params[key])) \
                if params[key] is not None and params[key] != '' else 'NULL'

            data[new_key] = new_val

        return data

    def __format_conditions(self, conds):
        conds = self.__format_data(conds)
        sub = [u"{}={}".format(key, conds[key]) for key in conds]

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
        sql = u"DELETE FROM {table} ".format(table=table)

        if conds:
            sql += self.__format_conditions(conds)

        return self.cursor.execute(sql)

    def __count_records(self, table, conds=None):
        """
        Helper function for counting records
        """
        sql = u"SELECT COUNT(*) FROM {table} ".format(table=table)

        if conds:
            sql += self.__format_conditions(conds)

        self.cursor.execute(sql)

    @staticmethod
    def create_database(species):
        """
        Create an empty database
        """
        cnx = mysql.connector.connect(**DBConn.db_config)
        cursor = cnx.cursor()
        name = DBConn.__get_name(species)

        cursor.execute(u"CREATE DATABASE `{}`".format(name))

        return name

    def execute_file(self, sql_file):
        """
        Execute multiple SQL queries from a file.
        """
        with open(sql_file, 'r') as fin:
            self.cursor.execute(fin.read().decode('utf-8'))
            while True:
                # keep executing until all queries are done
                if not self.cursor.nextset():
                    break

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
        try:
            self.cursor.execute(sql)

        except Exception as e:
            # dump the record before throwing the exception
            print("ERROR: db_conn.save_records()")
            pprint(sql)
            raise e

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
            'update': u", ".join([u"{}={}".format(key, formatted[key]) for key in formatted
                                  if key != u'`id`'])
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
                        itertools.islice(itertools.chain([first], records), self.max_insert_size))
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
            pprint(sql)
            raise e

    def execute_sql(self, sql):
        """
        Execute a given SQL query, and return the time it took.
        """
        start = time()

        self.cursor.execute(sql)
        self.cnx.commit()

        return timedelta(seconds=time() - start)
