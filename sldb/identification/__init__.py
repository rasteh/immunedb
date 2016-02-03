import re

import dnautils

from sldb.common.log import logger
from sldb.common.celery_app import get_result
from sldb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                Sequence)
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups


class AlignmentException(Exception):
    pass

ALIGNMENT_COPY_FIELDS = (
    'probable_indel_or_misalign', 'num_gaps', 'pad_length', 'v_match',
    'v_length', 'j_match', 'j_length', 'removed_prefix',
    'removed_prefix_qual', 'v_mutation_fraction', 'pre_cdr3_match',
    'pre_cdr3_match', 'post_cdr3_length', 'post_cdr3_match',
    'in_frame', 'functional', 'stop', 'copy_number', 'cdr3_nt',
    'cdr3_aa', 'cdr3_num_nts', 'sequence', 'quality', 'germline',
    'probable_indel_or_misalign'
)

class SequenceAlignment(object):
    INDEL_MISMATCH_THRESHOLD = .6
    INDEL_WINDOW = 30

    __slots__ = [
        'ids', 'germline', 'sequence', 'quality', 'pending_result', 'j_genes',
        'j_length', 'j_match', 'j_anchor_pos', 'v_genes', 'v_length',
        'v_match', 'germline_offset', 'v_mutation_fraction', 'removed_prefix',
        'removed_prefix', 'removed_prefix_qual', 'cdr3_start', 'cdr3_nt',

        'pre_cdr3_match', 'pre_cdr3_length', 'post_cdr3_match',
        'post_cdr3_length',
    ]

    def __init__(self, ids, sequence, quality=None):
        self.ids = ids
        self.sequence = sequence
        self.quality = quality
        self.removed_prefix = None
        self.removed_prefix_qual = None
        self.cdr3_start = CDR3_OFFSET

    @property
    def pad_length(self):
        return max(0, self.germline_offset)

    @property
    def copy_number(self):
        return len(self.ids)

    @property
    def cdr3_aa(self):
        return lookups.aas_from_nts(self.cdr3_nt)

    @property
    def cdr3_num_nts(self):
        return len(self.cdr3_nt)

    @property
    def v_similarity(self):
        return self.v_match / float(self.v_length)

    @property
    def num_gaps(self):
        return self.sequence[:CDR3_OFFSET].count('-')

    @property
    def in_frame(self):
        return len(self.cdr3_nt) % 3 == 0 and self.cdr3_start % 3 == 0

    @property
    def stop(self):
        return lookups.has_stop(self.sequence)

    @property
    def functional(self):
        return self.in_frame and not self.stop

    def await_result(self):
        return get_result(self.pending_result)

    @property
    def probable_indel_or_misalign(self):
        start = re.search('[ATCG]', self.sequence).start()
        germ = self.germline[start:self.cdr3_start]
        seq = self.sequence[start:self.cdr3_start]

        for i in range(0, len(germ) - self.INDEL_WINDOW + 1):
            dist = dnautils.hamming(germ[i:i+self.INDEL_WINDOW],
                                    seq[i:i+self.INDEL_WINDOW])
            if dist >= self.INDEL_MISMATCH_THRESHOLD * self.INDEL_WINDOW:
                return True

        return False


def add_as_noresult(session, sample, alignment, reason):
    try:
        session.bulk_save_objects([
            NoResult(
                seq_id=seq_id,
                sample_id=sample.id,
                sequence=alignment.sequence,
                quality=alignment.quality,
                reason=reason,
            ) for seq_id in alignment.ids
        ])
    except ValueError as e:
        logger.warning('Could not add NoResult: %s', e)


def add_as_sequence(session, sample, alignment):
    try:
        seq = Sequence(
            sample=sample,
            subject_id=sample.subject.id,
            seq_id=alignment.ids[0],

            partial=alignment.pad_length > 0,

            v_gene=funcs.format_ties(alignment.v_genes, 'IGHV'),
            j_gene=funcs.format_ties(alignment.j_genes, 'IGHJ'),
        )
        for field in ALIGNMENT_COPY_FIELDS:
            setattr(seq, field, getattr(alignment, field))

        session.flush()
        try:
            session.bulk_save_objects([
                DuplicateSequence(
                    sample_id=sample.id,
                    seq_id=seq_id,
                    duplicate_seq_ai=seq.ai
                    ) for seq_id in alignment.ids[1:]
            ])
        except ValueError as e:
            logger.warning('Could not add DuplicateSequence: %s', e)
    except ValueError as e:
        logger.warning('Could not add Sequence: %s', e)
        add_as_noresult(session, sample, alignment, str(e))
