#!/usr/bin/env python
# -*- coding: utf-8 -*-

# third party modules
import luigi

# local modules
from alleletraj import utils
from alleletraj.db.conn import Database


class ExternalDatabaseSQL(utils.PipelineExternalTask):
    """
    External task dependency for SQL table definitions.
    """

    def output(self):
        return luigi.LocalTarget('data/alleletraj_database.sql')


class CreateDatabase(utils.MySQLTask):
    """
    Create and new db and add all the empty tables.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ExternalDatabaseSQL()

    def queries(self):
        sql_file = self.input()

        # create the empty db
        Database.create_database(self.species)

        # load the CREATE TABLE sql file
        self.dbc.execute_file(sql_file.path)


if __name__ == '__main__':
    luigi.run()
