import dnautils
import json
import logging
import os
import traceback

from Bio import SeqIO

from sqlalchemy.sql import exists

import sldb.common.config as config
from sldb.common.celery_app import local_worker
from sldb.common.log import logger, setup_default_logging
import sldb.common.modification_log as mod_log
from sldb.identification import add_as_noresult, add_as_sequence, get_result
from sldb.common.models import Sample, Sequence, Study, Subject
from sldb.identification import AlignmentException, SequenceAlignment
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines
from sldb.identification.tasks import (align_sequence, align_to_vties,
                                       collapse_sequences)
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups


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


def process_sample(path, session, meta, min_similarity, max_vties):
    study, sample = setup_sample(session, meta)
    unique_seqs = {}
    logger.info('Collapsing identical sequences')
    ftype = 'fasta' if path.endswith('.fasta') else 'fastq'
    for i, record in enumerate(SeqIO.parse(path, ftype)):
        seq = str(record.seq)

        try:
            unique_seqs[seq].ids.append(record.description)
        except KeyError:
            aln = SequenceAlignment(
                ids=[record.description],
                sequence=seq,
                quality=funcs.ord_to_quality(
                    record.letter_annotations.get('phred_quality')
                )
            )
            aln.pending_result = align_sequence.delay(aln)
            unique_seqs[seq] = aln

    mutations = []
    lengths = []
    logger.info('Waiting for workers to finish')
    for seq in unique_seqs.keys():
        try:
            all_ids = unique_seqs[seq].ids
            unique_seqs[seq] = unique_seqs[seq].await_result()
            unique_seqs[seq].ids = all_ids
            lengths.append(unique_seqs[seq].v_length)
            mutations.append(unique_seqs[seq].v_mutation_fraction)
        except AlignmentException as e:
            add_as_noresult(session, sample, unique_seqs[seq], str(e))
            del unique_seqs[seq]
    session.commit()

    logger.info('Queueing V-ties alignment')
    avg_mut = sum(mutations) / float(len(mutations))
    avg_len = sum(lengths) / float(len(lengths))

    sample.v_ties_mutations = avg_mut
    sample.v_ties_len = avg_len
    session.commit()

    for aln in unique_seqs.values():
        aln.pending_result = align_to_vties.delay(aln, avg_len, avg_mut)

    logger.info('Waiting for workers to finish')
    bucketed_seqs = {}
    for seq in unique_seqs.keys():
        try:
            unique_seqs[seq] = unique_seqs[seq].await_result()
            if (unique_seqs[seq].v_similarity < min_similarity or
                    len(unique_seqs[seq].v_genes) > max_vties):
                raise AlignmentException(
                    'V-match too low or too many V-ties'
                )
            bucket_key = (
                funcs.format_ties(unique_seqs[seq].v_genes, 'IGHV'),
                funcs.format_ties(unique_seqs[seq].j_genes, 'IGHJ'),
                unique_seqs[seq].cdr3_num_nts,
            )
            if bucket_key not in bucketed_seqs:
                bucketed_seqs[bucket_key] = {}

            try:
                bucketed_seqs[bucket_key][seq].ids += unique_seqs[seq].ids
            except KeyError:
                bucketed_seqs[bucket_key][seq] = unique_seqs[seq]
        except AlignmentException as e:
            add_as_noresult(session, sample, unique_seqs[seq], str(e))
    session.commit()

    collapsed_seqs = []
    for sequences in bucketed_seqs.values():
        collapsed_seqs.append(collapse_sequences.delay(sequences))

    for collapsed in collapsed_seqs:
        for aln in get_result(collapsed):
            add_as_sequence(session, sample, aln)
    session.commit()


def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    setup_default_logging()
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    samples = {}
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
            logger.critical('Metadata file not found')
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
                if args.warn_existing:
                    logger.warning('Sample %s already exists. Skipping.',
                                   meta.get('sample_name'))
                else:
                    logger.warning('Sample %s already exists. Cannot proceed.',
                                   meta.get('sample_name'))
                fail = True
            elif meta.get('sample_name') in samples:
                logger.critical('Sample %s exists more than once in metadata. '
                                'Cannot continue.', meta.get('sample_name'))
                return
            else:
                samples[meta.get('sample_name')] = {
                    'path': os.path.join(directory, fn),
                    'meta': meta
                }
        if fail and not args.warn_existing:
            logger.critical(
                'Encountered errors.  Not running any identification.  To '
                'skip samples that are already in the database use '
                '--warn-existing.'
            )
            return

    with local_worker(not args.no_worker, args.nproc):
        for task in samples.values():
            process_sample(
                session=session, min_similarity=args.min_similarity / 100.0,
                max_vties=args.max_vties, **task
            )
