import dnautils
import os
from functools import wraps

from Bio.Seq import Seq
from celery import Celery

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
app = Celery('identify', backend='redis://', broker='amqp://guest@localhost//')

def find_j(j_germlines, sequence, quality=None):
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
        j_anchor_pos = sequence.rfind(match)
        if j_anchor_pos >= 0:
            break

        rc_sequence = str(Seq(sequence).reverse_complement())
        j_anchor_pos = rc_sequence.rfind(match)
        if j_anchor_pos >= 0:
            sequence = rc_sequence
            if quality is not None:
                quality = quality[::-1]
            break
    else:
        raise AlignmentException('Could not find J anchor')

    # Get the full germline J gene
    j_full = j_germlines[j_gene]

    # Get the portion of the germline J in the CDR3
    germline_in_cdr3 = j_germlines.get_j_in_cdr3(j_gene)
    cdr3_end_pos = (
        j_anchor_pos + j_germlines.anchor_len -
        j_germlines.upstream_of_cdr3
    )
    sequence_in_cdr3 = sequence[cdr3_end_pos - len(germline_in_cdr3):
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
    j_start = j_anchor_pos + len(match) - len(j_full)

    # If the trimmed germline J extends past the end of the
    # sequence, there is a misalignment
    if len(j_full) != len(sequence[j_start:j_start+len(j_full)]):
        raise AlignmentException('Germline extended past end of J')

    return Map(
        anchor_pos=j_anchor_pos,
        genes=j_germlines.get_ties(j_gene, match),
        length=len(j_full)
    )


def find_v(sequence, j_alignment, v_germlines):
    for anchor_pos in find_v_position(sequence):
        aligned_v = VGene(sequence)
        best_vs = None
        v_score = None
        for v, germ in sorted(v_germlines.iteritems()):
            try:
                dist, total_length = germ.compare(
                    aligned_v, j_alignment.anchor_pos, MISMATCH_THRESHOLD
                )
            except:
                continue
            # Record this germline if it is has the lowest distance
            if dist is not None:
                if v_score is None or dist < v_score:
                    best_vs = [v]
                    v_length = total_length
                    germ_pos = germ.ungapped_anchor_pos
                    v_score = dist
                elif dist == v_score:
                    # Add the V-tie
                    best_vs.append(v)

        if best_vs is not None:
            # Determine the pad length
            pad_len = germ_pos - anchor_pos
            # Mutation ratio is the distance divided by the length of overlap
            mutation_fraction = round(v_score / float(v_length), 2)
            return Map(
                genes=best_vs,
                germ_pos=germ_pos,
                pad_len=pad_len,
                length=v_length,
                mutation_fraction=mutation_fraction
            )

    raise AlignmentException('Could not find suitable V anchor')


def align_to_germline(v_germlines, j_germlines, v_alignment, j_alignment,
                      sequence, quality, avg_len, avg_mut):
    if avg_len is not None and avg_mut is not None:
        best_vs = v_germlines.get_ties(v_alignment.genes, avg_len, avg_mut)
    # Set the germline to the V gene up to the CDR3
    germline = get_common_seq(
        [v_germlines[v].sequence for v in v_alignment.genes]
    )[:CDR3_OFFSET]
    # If we need to pad the sequence, do so, otherwise trim the sequence to
    # the germline length
    if v_alignment.pad_len >= 0:
        sequence = 'N' * v_alignment.pad_len + sequence
        if quality is not None:
            quality = (' ' * v_alignment.pad_len) + quality
    else:
        removed_prefix = sequence[:-v_alignment.pad_len]
        sequence = str(sequence[-v_alignment.pad_len:])
        if quality is not None:
            removed_prefix_qual = quality[:-v_alignment.pad_len]
            quality = quality[-v_alignment.pad_len:]
    # Update the anchor positions after adding padding / trimming
    j_alignment.anchor_pos += v_alignment.pad_len

    # Add germline gaps to sequence before CDR3 and update anchor positions
    for i, c in enumerate(germline):
        if c == '-':
            sequence = sequence[:i] + '-' + sequence[i:]
            if quality is not None:
                quality = quality[:i] + ' ' + quality[i:]
            j_alignment.anchor_pos += 1

    j_germ = get_common_seq(
        map(reversed, [j_germlines[j] for j in j_alignment.genes]))
    j_germ = ''.join(reversed(j_germ))
    # Calculate the length of the CDR3
    cdr3_len = (
        j_alignment.anchor_pos + j_germlines.anchor_len -
        j_germlines.upstream_of_cdr3 - CDR3_OFFSET
    )

    if cdr3_len < 3:
        raise AlignmentException('CDR3 has no AAs'.format(cdr3_len))

    j_alignment.anchor_pos += cdr3_len
    # Fill germline CDR3 with gaps
    germline += '-' * cdr3_len
    germline += j_germ[-j_germlines.upstream_of_cdr3:]
    # If the sequence is longer than the germline, trim it
    if len(sequence) > len(germline):
        sequence = sequence[:len(germline)]
        if quality is not None:
            quality = quality[:len(germline)]
    elif len(sequence) < len(germline):
        sequence += 'N' * (len(germline) - len(sequence))
        if quality is not None:
            quality += ' ' * (len(germline) - len(quality))

    # Get the pre-CDR3 germline
    pre_cdr3_germ = germline[:CDR3_OFFSET]
    pre_cdr3_seq = sequence[:CDR3_OFFSET]

    # If there is padding, get rid of it in the sequence and align the
    # germline
    if v_alignment.pad_len > 0:
        pre_cdr3_germ = pre_cdr3_germ[v_alignment.pad_len:]
        pre_cdr3_seq = pre_cdr3_seq[v_alignment.pad_len:]

    # Calculate the pre-CDR3 length and distance
    pre_cdr3_length = len(pre_cdr3_seq)
    pre_cdr3_match = pre_cdr3_length - dnautils.hamming(
        str(pre_cdr3_seq), str(pre_cdr3_germ))

    # Get the length of J after the CDR3
    post_cdr3_length = j_germlines.upstream_of_cdr3
    # Get the sequence and germline sequences after CDR3
    post_j = j_germ[-post_cdr3_length:]
    post_s = sequence[-post_cdr3_length:]

    # Calculate their match count
    post_cdr3_match = post_cdr3_length - dnautils.hamming(
        post_j, post_s)

    v_match = v_alignment.length - dnautils.hamming(
        germline[:CDR3_OFFSET],
        sequence[:CDR3_OFFSET]
    )

    j_match = j_alignment.length - dnautils.hamming(
        germline[-len(j_germ):],
        sequence[-len(j_germ):]
    )

    return Map(
        sequence=sequence,
        germline=germline,
        cdr3_len=cdr3_len,
        pre_cdr3_match=pre_cdr3_match,
        pre_cdr3_length=pre_cdr3_length,
        post_cdr3_match=pre_cdr3_match,
        post_cdr3_length=pre_cdr3_length,
    )


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
def align_sequence(v_germlines, j_germlines, sequence, quality=None):
    j_alignment = find_j(j_germlines, sequence, quality)
    v_alignment = find_v(sequence, j_alignment, v_germlines)
    return {
        'v': v_alignment.to_dict(),
        'j': j_alignment.to_dict()
    }

@app.task
@with_germlines
def align_to_vties(v_germlines, j_germlines, v_alignment, j_alignment,
                   sequence, quality, avg_len, avg_mut):
    final = align_to_germline(v_germlines, j_germlines, Map(v_alignment),
                              Map(j_alignment), sequence, quality, avg_len,
                              avg_mut)
    return final.to_dict()
