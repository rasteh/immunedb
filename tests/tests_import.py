from regression import NamespaceMimic, RegressionTest

from sldb.importing.delimited import DEFAULT_MAPPINGS, run_import

class TestImport(RegressionTest):
    def testAll(self):
        self.imgt_import()

    def imgt_import(self):
        ns = NamespaceMimic(
            input_file='tests/data/imgt/input/'
               '2_IMGT-gapped-nt-sequences.txt',
            v_germlines='tests/imgt_human_v.fasta',
            j_germlines='tests/imgt_human_j.fasta',
            upstream_of_cdr3=31,
            anchor_len=18,
            min_anchor_len=12,
            study_name='Test',
            sample_name='Sample 1',
            subject='Subject 1',
            date='2016-01-01',
            ties=True,
            tissue='PBL',
            subset='MN',
            disease='Something bad',
            lab='Sanger',
            experimenter='Bob',
            ig_class='IgE',
            v_primer='Leader',
            j_primer='Trailer',
            unpaired=False,
        )
        for k, v in DEFAULT_MAPPINGS.iteritems():
            setattr(ns, k, v)

        run_import(self.session, ns)
        self.session.commit()
        populate_regression()
