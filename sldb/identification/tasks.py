import dnautils
import os
from functools import wraps

from Bio.Seq import Seq

from sldb.common.celery_app import app
from sldb.common.models import CDR3_OFFSET
from sldb.identification import AlignmentException
from sldb.identification.v_genes import (get_common_seq, find_v_position,
                                         VGene, VGermlines)
from sldb.identification.j_genes import JGermlines
from sldb.util import Map
from sldb.util.funcs import find_streak_position

MISMATCH_THRESHOLD = 3
v_germlines = None
j_germlines = None

def find_j(j_germlines, aln):
    '''Finds the location and type of J gene'''
    # Iterate over every possible J anchor.  For each germline, try its
    # full sequence, then exclude the final 3 characters at a time until
    # there are only MIN_J_ANCHOR_LEN nucleotides remaining.
    #
    # For example, the order for one germline:
    # TGGTCACCGTCTCCTCAG
    # TGGTCACCGTCTCCT
    # TGGTCACCGTCT

    for match, full_anchor, j_gene in j_germlines.get_all_anchors():
        aln.j_anchor_pos = aln.sequence.rfind(match)
        if aln.j_anchor_pos >= 0:
            break

        rc_sequence = str(Seq(aln.sequence).reverse_complement())
        aln.j_anchor_pos = rc_sequence.rfind(match)
        if aln.j_anchor_pos >= 0:
            aln.sequence = rc_sequence
            if aln.quality is not None:
                aln.quality = aln.quality[::-1]
            break
    else:
        raise AlignmentException('Could not find J anchor')

    # Get the full germline J gene
    j_full = j_germlines[j_gene]

    # Get the portion of the germline J in the CDR3
    germline_in_cdr3 = j_germlines.get_j_in_cdr3(j_gene)
    cdr3_end_pos = (
        aln.j_anchor_pos + j_germlines.anchor_len -
        j_germlines.upstream_of_cdr3
    )
    sequence_in_cdr3 = aln.sequence[cdr3_end_pos - len(germline_in_cdr3):
                                          cdr3_end_pos]
    if len(germline_in_cdr3) == 0 or len(sequence_in_cdr3) == 0:
        raise AlignmentException('Could not find sequence or germline in '
                                 'CDR3')

    # Get the extent of the J in the CDR3
    streak = find_streak_position(
        reversed(germline_in_cdr3),
        reversed(sequence_in_cdr3),
        MISMATCH_THRESHOLD)

    # Trim the J gene based on the extent in the CDR3
    if streak is not None:
        j_full = j_full[len(germline_in_cdr3) - streak:]

    # Find where the full J starts
    j_start = aln.j_anchor_pos + len(match) - len(j_full)

    # If the trimmed germline J extends past the end of the
    # sequence, there is a misalignment
    if len(j_full) != len(aln.sequence[j_start:j_start+len(j_full)]):
        raise AlignmentException('Germline extended past end of J')

    aln.j_genes = j_germlines.get_ties(j_gene, match)
    aln.j_length = len(j_full)
    return aln


def find_v(v_germlines, aln):
    aligned_v = VGene(aln.sequence)
    for anchor_pos in find_v_position(aln.sequence):
        aln.v_genes = None
        v_score = None
        for v, germ in sorted(v_germlines.iteritems()):
            try:
                dist, total_length = germ.compare(
                    aligned_v, aln.j_anchor_pos, MISMATCH_THRESHOLD
                )
            except:
                continue
            # Record this germline if it is has the lowest distance
            if dist is not None:
                if v_score is None or dist < v_score:
                    aln.v_genes = [v]
                    aln.v_length = total_length
                    germ_pos = germ.ungapped_anchor_pos
                    v_score = dist
                elif dist == v_score:
                    # Add the V-tie
                    aln.v_genes.append(v)

        if aln.v_genes is not None:
            # Determine the pad length
            aln.germline_offset = germ_pos - anchor_pos
            # Mutation ratio is the distance divided by the length of overlap
            aln.v_mutation_fraction = round(v_score / float(aln.v_length), 2)
            return aln

    raise AlignmentException('Could not find suitable V anchor')


def with_germlines(func):
    if os.environ.get('V_PATH'):
        v_germlines = VGermlines(os.environ.get('V_PATH'))
        j_germlines = JGermlines(
            os.environ.get('J_PATH'),
            int(os.environ.get('UPSTREAM')),
            int(os.environ.get('ANCHOR_LEN')),
            int(os.environ.get('MIN_ANCHOR_LEN'))
        )

        @wraps(func)
        def _wrapper(*args, **kwargs):
            return func(v_germlines, j_germlines, *args, **kwargs)
        return _wrapper

    return func


@app.task
@with_germlines
def align_sequence(v_germlines, j_germlines, aln):
    aln = find_j(j_germlines, aln)
    return find_v(v_germlines, aln)

@app.task
@with_germlines
def align_to_vties(v_germlines, j_germlines, aln, avg_len, avg_mut):
    aln.v_genes = v_germlines.get_ties(aln.v_genes, avg_len, avg_mut)
    # Set the germline to the V gene up to the CDR3
    aln.germline = get_common_seq(
        [v_germlines[v].sequence for v in aln.v_genes]
    )[:CDR3_OFFSET]
    # If we need to pad the sequence, do so, otherwise trim the sequence to
    # the germline length
    if aln.germline_offset >= 0:
        aln.sequence = 'N' * aln.germline_offset + aln.sequence
        if aln.quality is not None:
            aln.quality = ' ' * aln.germline_offset + aln.quality
    else:
        aln.removed_prefix = aln.sequence[:-aln.germline_offset]
        aln.sequence = aln.sequence[-aln.germline_offset:]
        if aln.quality is not None:
            aln.removed_prefix_qual = aln.quality[:-aln.germline_offset]
            aln.quality = aln.quality[-aln.germline_offset:]
    # Update the anchor positions after adding padding / trimming
    aln.j_anchor_pos += aln.germline_offset

    # Add germline gaps to sequence before CDR3 and update anchor positions
    for i, c in enumerate(aln.germline):
        if c == '-':
            aln.sequence = aln.sequence[:i] + '-' + aln.sequence[i:]
            if aln.quality is not None:
                aln.quality = aln.quality[:i] + ' ' + aln.quality[i:]
            aln.j_anchor_pos += 1

    j_germ = get_common_seq(
        map(reversed, [j_germlines[j] for j in aln.j_genes]))
    j_germ = ''.join(reversed(j_germ))

    # Calculate the length of the CDR3
    cdr3_len = (
        aln.j_anchor_pos + j_germlines.anchor_len -
        j_germlines.upstream_of_cdr3 - CDR3_OFFSET
    )

    if cdr3_len < 3:
        raise AlignmentException('CDR3 has no AAs')

    aln.j_anchor_pos += cdr3_len
    # Fill germline CDR3 with gaps
    aln.germline += '-' * cdr3_len
    aln.germline += j_germ[-j_germlines.upstream_of_cdr3:]
    # If the sequence is longer than the germline, trim it
    if len(aln.sequence) > len(aln.germline):
        aln.sequence = aln.sequence[:len(aln.germline)]
        if aln.quality is not None:
            aln.quality = aln.quality[:len(aln.germline)]
    elif len(aln.sequence) < len(aln.germline):
        aln.sequence += 'N' * (len(aln.germline) - len(aln.sequence))
        if aln.quality is not None:
            aln.quality += ' ' * (len(aln.germline) - len(aln.quality))

    # Get the pre-CDR3 germline
    pre_cdr3_germ = aln.germline[:CDR3_OFFSET]
    pre_cdr3_seq = aln.sequence[:CDR3_OFFSET]

    # If there is padding, get rid of it in the sequence and align the
    # germline
    if aln.germline_offset > 0:
        pre_cdr3_germ = pre_cdr3_germ[aln.germline_offset:]
        pre_cdr3_seq = pre_cdr3_seq[aln.germline_offset:]

    # Calculate the pre-CDR3 length and distance
    pre_cdr3_length = len(pre_cdr3_seq)
    pre_cdr3_match = pre_cdr3_length - dnautils.hamming(
        str(pre_cdr3_seq), str(pre_cdr3_germ))

    # Get the length of J after the CDR3
    post_cdr3_length = j_germlines.upstream_of_cdr3
    # Get the sequence and germline sequences after CDR3
    post_j = j_germ[-post_cdr3_length:]
    post_s = aln.sequence[-post_cdr3_length:]

    # Calculate their match count
    post_cdr3_match = post_cdr3_length - dnautils.hamming(post_j, post_s)

    aln.pre_cdr3_match = pre_cdr3_match
    aln.pre_cdr3_length = pre_cdr3_length
    aln.post_cdr3_match = post_cdr3_match
    aln.post_cdr3_length = post_cdr3_length

    aln.v_match = aln.v_length - dnautils.hamming(
        aln.germline[:CDR3_OFFSET],
        aln.sequence[:CDR3_OFFSET]
    )

    aln.j_match = aln.j_length - dnautils.hamming(
        aln.germline[-len(j_germ):],
        aln.sequence[-len(j_germ):]
    )
    aln.v_mutation_fraction = dnautils.hamming(aln.germline, aln.sequence)
    aln.cdr3_nt = aln.sequence[CDR3_OFFSET:CDR3_OFFSET + cdr3_len]

    return aln


@app.task
def collapse_sequences(alignments):
    alignments = sorted(alignments.values(), cmp=lambda a, b:
                       cmp(a.copy_number, b.copy_number), reverse=True)
    current_largest = 0
    while current_largest < len(alignments):
        larger = alignments[current_largest]
        for i in reversed(range(current_largest + 1, len(alignments))):
            smaller = alignments[i]
            if dnautils.equal(larger.sequence, smaller.sequence):
                larger.ids += smaller.ids
                del alignments[i]
        current_largest += 1

    return alignments
