from collections import OrderedDict

from celery.exceptions import TimeoutError
import dnautils
import json
import os
import traceback

from Bio import SeqIO

from sqlalchemy.sql import exists

import sldb.common.config as config
from sldb.common.log import logger
import sldb.common.modification_log as mod_log
from sldb.identification import add_as_noresult, add_as_sequence
from sldb.common.models import (DuplicateSequence, NoResult, Sample, Sequence,
                                Study, Subject)
from sldb.identification import AlignmentException
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines
from sldb.identification.tasks import align_sequence, align_to_vties
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups

def get_result(task, max_tries=None, timeout=.5):
    tries = 0
    while True:
        if max_tries is not None and tries > max_tries:
            raise TimeoutError()

        try:
            return task.get(timeout=timeout)
        except TimeoutError:
            tries += 1


class SampleMetadata(object):
    def __init__(self, specific_config, global_config=None):
        self._specific = specific_config
        self._global = global_config

    def get(self, key, require=True):
        if key in self._specific:

            return self._specific[key]
        elif self._global is not None and key in self._global:
            return self._global[key]
        if require:
            raise Exception(('Could not find metadata for key {}'.format(key)))


def setup_sample(session, meta):
    study, new = funcs.get_or_create(
        session, Study, name=meta.get('study_name'))

    if new:
        logger.info('Created new study "%s"', study.name)
        session.commit()

    name = meta.get('sample_name')
    sample, new = funcs.get_or_create(session, Sample, name=name, study=study)
    if new:
        sample.date = meta.get('date')
        logger.info('Created new sample "%s"', sample.name)
        for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                    'ig_class', 'v_primer', 'j_primer'):
            setattr(sample, key, meta.get(key, require=False))
        subject, new = funcs.get_or_create(
            session, Subject, study=study,
            identifier=meta.get('subject'))
        sample.subject = subject
        session.commit()

    return study, sample


def process_sample(path, session, sample, meta):
    unique_seqs = OrderedDict()
    ftype = 'fasta' if path.endswith('.fasta') else 'fastq'

    print 'Collapsing identical sequences'
    for i, record in enumerate(SeqIO.parse(path, ftype)):
        if i > 100:
            break
        seq = str(record.seq)
        if seq not in unique_seqs:
            quality = funcs.ord_to_quality(
                record.letter_annotations.get('phred_quality')
            )
            unique_seqs[seq] = {
                'ids': [],
                'quality': quality,
                'align_result': align_sequence.delay(seq, quality)
            }
        unique_seqs[seq]['ids'].append(record.description)

    mutations = []
    lengths = []
    print 'Waiting for workers to finish'
    for seq in unique_seqs.keys():
        info = unique_seqs[seq]
        try:
            info['align_result'] = get_result(info['align_result'])
            lengths.append(info['align_result']['v']['length'])
            mutations.append(info['align_result']['v']['mutation_fraction'])
        except AlignmentException as e:
            #add_as_noresult(session, sample, info['ids'], seq, info['quality'])
            del unique_seqs[seq]

    avg_mut = sum(mutations) / float(len(mutations))
    avg_len = sum(lengths) / float(len(lengths))
    for seq, info in unique_seqs.iteritems():
        info['vties_result'] = align_to_vties.delay(
            info['align_result']['v'],
            info['align_result']['j'],
            seq,
            info['quality'],
            avg_len,
            avg_mut
        )

    for seq, info in unique_seqs.iteritems():
        try:
            info['vties_result'] = get_result(info['vties_result'])
            print info['vties_result']
        except AlignmentException as e:
            pass

def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)
    tasks = []

    sample_names = set([])
    fail = False
    for directory in args.sample_dirs:
        # If metadata is not specified, assume it is "metadata.json" in the
        # directory
        if args.metadata is None:
            meta_fn = os.path.join(directory, 'metadata.json')
        else:
            meta_fn = args.metadata

        # Verify the metadata file exists
        if not os.path.isfile(meta_fn):
            print 'Metadata file not found.'
            return

        with open(meta_fn) as fh:
            metadata = json.load(fh)

        # Create the tasks for each file
        for fn in sorted(metadata.keys()):
            if fn == 'all':
                continue
            meta = SampleMetadata(
                metadata[fn],
                metadata['all'] if 'all' in metadata else None)
            if session.query(Sample).filter(
                    Sample.name == meta.get('sample_name'),
                    exists().where(
                        Sequence.sample_id == Sample.id
                    )).first() is not None:
                print 'Sample {} already exists. {}'.format(
                    meta.get('sample_name'), 'Skipping.' if
                    args.warn_existing else 'Cannot continue.'
                )
                fail = True
            elif meta.get('sample_name') in sample_names:
                print ('Sample {} exists more than once in metadata. Cannot '
                       'continue.').format(meta.get('sample_name'))
                return
            else:
                tasks.append({
                    'path': os.path.join(directory, fn),
                    'meta': meta
                })
                sample_names.add(meta.get('sample_name'))
        if fail and not args.warn_existing:
            print ('Encountered errors.  Not running any identification.  To '
                   'skip samples that are already in the database use '
                   '--warn-existing.')
            return


    for task in tasks:
        study, sample = setup_sample(session, meta)
        process_sample(session=session, sample=sample, **task)
