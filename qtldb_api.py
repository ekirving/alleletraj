#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import urllib
import urllib2
import xmltodict
from time import sleep

from collections import OrderedDict
from xml.sax.saxutils import escape

QTLDB_API_URL = 'http://www.animalgenome.org/cgi-bin/QTLdb/API'
QTLDB_CHUNK_SIZE = 10

QTL_FILES = {
    'pig':    'data/AnimalQTLdb/rel35/cM/pig.txt',
    'cattle': 'data/AnimalQTLdb/rel35/cM/cattle.txt',
    'horse':  'data/AnimalQTLdb/rel35/cM/horse.txt',
}

# wait 10 seconds beforey retyring a failed request
QTLDB_WAIT_TIME = 10

class qtldb_api:

    @staticmethod
    def __execute(request):
        """
        Execute an API request.
        """
        try:
            xml = urllib2.urlopen(request).read().decode("utf-8", 'ignore')

        except urllib2.URLError as e:
            # wait a while and try again
            sleep(QTLDB_WAIT_TIME)

            xml = urllib2.urlopen(request).read().decode("utf-8", 'ignore')


        result = None

        loops = 0

        while not result:
            try:
                result = xmltodict.parse(xml, xml_attribs=True)

            except xmltodict.expat.ExpatError as err:
                # find the problem character and encode it
                line, column = [int(s) for s in re.split('[ ,]', str(err)) if s.isdigit()]

                lines = xml.splitlines(True)
                problem = list(lines[line-1])

                # handle entity encoding issues
                problem[column - 1] = escape(problem[column - 1])

                # if ord(problem[column]) > 160:
                #     # handle utf8 encoding issues
                #     problem[column] = '&#%s;' % ord(problem[column])

                lines[line-1] = ''.join(problem)
                xml = ''.join(lines)

                if loops > 100:
                    # don't keep looping forever!
                    raise err

            loops += 1

        return result

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
        request = "{url}/iquery?s={s}&q={q}&h={h}".format(url=QTLDB_API_URL, s=species,
                                                          q=urllib.quote_plus(query), h=type)

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
            query = ','.join(map(str, chunk))

            # fetch the results
            data = self.__ifetch(species, query)

            for record in data['EFETCHresults']['QTL']:
                yield record

    def get_trait_type(self, species, trait_id, name):

        # get all trait
        data = self.__iquery(species, name, 'traits')

        if 'trait' in data['EqueryResults']:
            traits = data['EqueryResults']['trait']
            if type(traits) is not list:
                traits = [traits]

            for trait in traits:
                if isinstance(trait, OrderedDict):
                    if trait['@traitID'] == trait_id:
                        return trait['traitType']
        return None

    def get_publication(self, species, pubmedid):

        # get all trait
        data = self.__iquery(species, pubmedid, 'publications')

        data['EqueryResults'].pop('@xmlns')

        return data['EqueryResults']