#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import urllib
import urllib2
import xmltodict

QTLDB_API_URL = 'http://www.animalgenome.org/cgi-bin/QTLdb/API'
QTLDB_CHUNK_SIZE = 10

QTLDB_TRAIT_CLASSES = [
    'Meat and Carcass',
    'Health',
    'Exterior',
    'Production',
    'Reproduction'
]

# TODO remove me
import pprint

class qtldb_api:

    def __execute(self, request):
        """
        Execute an API request.
        """
        xml = urllib2.urlopen(request).read()

        # fix entity issue
        xml = re.sub(r'\s&\s', '&amp;', xml)

        print request

        return xmltodict.parse(xml, xml_attribs=True)

    def __iinfo(self, species):
        """
        Fetch a summary record for a species.
        """
        # compose the request
        request = "{url}/iinfo?s={s}&q=info".format(url=QTLDB_API_URL, s=species)

        return self.__execute(request)

    def __iquery(self, species, query, type):
        """
        Fetch records by free-text query.
        """
        request = "{url}/iquery?s={s}&q={q}&h={h}".format(url=QTLDB_API_URL, s=species, q=urllib.quote_plus(query), h=type)

        return self.__execute(request)

    def __ifetch(self, species, query):
        """
        Fetch QTL records by ID list.
        """
        request = "{url}/ifetch?s={s}&q={q}".format(url=QTLDB_API_URL, s=species, q=urllib.quote_plus(query))

        return self.__execute(request)

    def get_qtls(self, species, ids):

        # fetch all the QTL records in chunks
        for i in range(0, len(ids), QTLDB_CHUNK_SIZE):

            chunk = ids[i:i + QTLDB_CHUNK_SIZE]

            # make a comma separated list
            query = ','.join(chunk)

            # fetch the results
            data = self.__ifetch(species, query)

            for record in data['EFETCHresults']['QTL']:
                yield record

    def get_traits(self, species):

        traits = []

        for tclass in QTLDB_TRAIT_CLASSES:

            # get all the traits
            data = self.__iquery(species, tclass, 'traits')

            if 'trait' in data['EqueryResults']:
                for trait in data['EqueryResults']['trait']:
                    yield trait

