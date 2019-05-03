#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

# import my custom modules
from alleletraj.utils import PipelineTask, PipelineExternalTask
from alleletraj.database.api import Database


class ExternalDatabaseSQL(PipelineExternalTask):
    """
    External task dependency for SQL table definitions.
    """
    def output(self):
        return luigi.LocalTarget('alleletraj_database.sql')


class CreateDatabase(PipelineTask):
    """
    Create and new database and add all the empty tables.

    :type species: str
    """
    species = luigi.Parameter()

    def requires(self):
        return ExternalDatabaseSQL()

    def output(self):
        return luigi.LocalTarget('db/{}-{}.log'.format(self.basename, self.classname))

    def run(self):
        sql_file = self.input()

        # create the empty database
        name = Database.create_database(self.species)

        # open a connection to the new db
        dbc = self.db_conn()

        # load the CREATE TABLE sql file
        dbc.execute_file(sql_file.path)

        with self.output().open('w') as fout:
            fout.write('INFO: Created database `{}`'.format(name))


if __name__ == '__main__':
    luigi.run()
