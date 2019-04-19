#!/usr/bin/env python
# -*- coding: utf-8 -*-

import luigi

# import my custom modules
from pipeline_utils import PipelineTask
from dbconn import DBConn


class CreateDatabase(PipelineTask):
    """
    Create and new database and add all the empty tables.

    :type species: str
    """
    species = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget('db/{}-database.log'.format(self.species))

    def run(self):
        # create the empty database
        name = DBConn.create_database(self.species)

        # open a connection to the new db
        dbc = self.db_conn()

        # load the CREATE TABLE sql file
        dbc.execute_file('alleletraj_database.sql')

        with self.output().open('w') as fout:
            fout.write('INFO: Created database `{}`'.format(name))


if __name__ == '__main__':
    luigi.run()
