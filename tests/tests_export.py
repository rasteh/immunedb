import os
import requests
import unittest

import immunedb.common.config as config
from immunedb.exporting.clone_export import CloneExport
from immunedb.exporting.sequence_export import SequenceExport


def request(endpoint, data={}):
    return requests.post('http://localhost:8891' + endpoint, data=data)


class ExportTest(unittest.TestCase):
    def check(self, expected_path, url, data={}):
        path = 'tests/data/export/' + expected_path
        response = request(url, data).text
        if os.getenv('GENERATE'):
            with open(path, 'w+') as fh:
                fh.write(response)
        else:
            with open(path) as fh:
                expected = fh.read()
                if expected != response:
                    with open('expected.log', 'w+') as fh:
                        fh.write(expected)
                    with open('response.log', 'w+') as fh:
                        fh.write(response)
                assert expected == response

    def test_sequences(self):
        for ext in ('fasta', 'fastq', 'clip'):
            self.check(
                'sequences.' + ext,
                '/export/sequences/sample/T2', data={
                    'format': ext,
                    'fields': ','.join(sorted(
                        SequenceExport.allowed_fields.keys()))
                })

    def test_clones(self):
        self.check(
            'clones.csv',
            '/export/clones/sample/T2', data={
                'fields': ','.join(sorted(CloneExport.allowed_fields.keys()))
            })

    def test_mutations(self):
        self.check(
            'mutations.csv',
            '/export/mutations/sample/T2', data={
                'fields': ','.join(sorted(CloneExport.allowed_fields.keys())),
            })
