import itertools
import multiprocessing as mp
import re
import subprocess
import shlex
from sqlalchemy import func

import dnautils
import sldb.util.funcs as funcs
from sldb.identification import AlignmentException
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import VGermlines
from sldb.identification.vdj_sequence import VDJSequence
from sldb.common.models import (HashExtension, DuplicateSequence, NoResult,
                                Sample, Sequence)
import sldb.util.lookups as lookups
import sldb.util.concurrent as concurrent


def gaps_before(gaps, pos):
    return sum((e[1] for e in gaps if e[0] < pos))

def gap_positions(seq):
    gaps = []
    for diff in re.finditer('[-]+', seq):
        start, end = diff.span()
        gaps.append((start, end - start))
    return gaps


class LocalAlignmentWorker(concurrent.Worker):
    def __init__(self, v_germlines, j_germlines, align_path, max_deletions,
                 max_insertions):
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._align_path = align_path
        self._max_deletions = max_deletions
        self._max_insertions = max_insertions

        self._first_alleles = {
            name: v.sequence for name, v in self._v_germlines.iteritems()
            if int(name.split('*', 1)[1]) == 1
        }

    def do_task(self, args):
        # Find best aligned first allele
        v_align = self._align_seq_to_germs(args['seq_id'], args['seq'],
                                           self._first_alleles)

        if v_align is None or not self._alignment_passes(v_align['germ'],
                                                         v_align['seq']):
            print 'Bad V'
            return

        v_name = v_align['germ_name'].split('*', 1)[0]
        v_ties = {
            '|'.join(name): v for name, v in self._v_germlines.all_ties(
                args['avg_len'], args['avg_mut'], cutoff=False
            ).iteritems() if v_name in '|'.join(name)
        }

        v_align = self._align_seq_to_germs(
            args['seq_id'], v_align['seq'].replace('-', ''), v_ties
        )

        germ_insertions = gap_positions(v_align['germ'])
        imgt_gaps = [
            i + gaps_before(germ_insertions, i)
            for i, c in enumerate(v_ties[v_align['germ_name']]) if c == '-'
        ]
        for gap in imgt_gaps:
            v_align['germ'] = ''.join((
                v_align['germ'][:gap],
                '^',
                v_align['germ'][gap:]
            ))
            v_align['seq'] = ''.join((
                v_align['seq'][:gap],
                '-',
                v_align['seq'][gap:]
            ))

        cdr3_start = 0
        count = 0
        for i, c in enumerate(v_align['germ']):
            if c != '-':
                count += 1
                if count > 309:
                    cdr3_start = i
                    break
        v_align['germ'] = v_align['germ'].replace('^', '-')

        j_align = self._align_seq_to_germs(
            args['seq_id'], v_align['seq'][cdr3_start:], self._j_germlines,
        )

        if j_align is None:
            print 'Bad J'
            return

        final_germ = ''.join((
            v_align['germ'][:-len(j_align['germ'])],
            j_align['germ']
        ))
        final_seq = ''.join((
            v_align['seq'][:-len(j_align['seq'])],
            j_align['seq']
        ))
        final_germ = final_germ.rstrip('-')
        final_seq = final_seq[:len(final_germ)]
        cdr3_end = len(final_seq) - self._j_germlines.upstream_of_cdr3
        final_germ = ''.join((
            final_germ[:cdr3_start],
            '-' * (cdr3_end - cdr3_start),
            final_germ[cdr3_end:]
        ))

    def _align_seq_to_germs(self, seq_id, seq, germs):
        stdin = []
        for g_name, g_seq in germs.iteritems():
            stdin.append('>{}\n{}\n'.format(g_name, g_seq.replace('-', '')))
            stdin.append('>{}\n{}\n'.format(seq_id, seq))

        cmd = (
            '{} --match 2 --mismatch -2 --gapopen -5 --gapextend -2 '
            '--wildcard N 2 --printfasta --printscores --freestartgap '
            '--freeendgap --file -'
        ).format(self._align_path)
        proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        output, err = proc.communicate(''.join(stdin))
        regex = (
            r'(?P<germ_name>.+)\n'
            r'(?P<germ_padding>[-]*)(?P<germ>[ATCGN-]+)\n'
            r'.+\n'
            r'(?P<seq_padding>[-]*)(?P<seq>[ATCGN-]+)\n'
            r'score: (?P<score>\d+)\n'
        )

        best = None
        for match in re.finditer(regex, output, re.MULTILINE):
            match = match.groupdict()
            match['germ_offset'] = len(match['germ_padding'])
            match['seq_offset'] = len(match['seq_padding'])
            match['score'] = int(match['score'])
            if best is None or match['score'] > best['score']:
                best = match

        if best is None:
            return None
        if best['germ_offset'] > 0:
            best['seq'] = best['seq'][best['germ_offset']:]
        elif best['seq_offset'] > 0:
            best['seq'] = ('N' * best['seq_offset']) + best['seq']

        return best

    def _alignment_passes(self, germ, seq):
        if len(gap_positions(seq.strip('-'))) > self._max_deletions:
            return False
        if len(gap_positions(germ.strip('-'))) > self._max_insertions:
            return False
        return True


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3, 0, 0)

    mutation_cache = {}
    tasks = concurrent.TaskQueue()

    '''
    indels = session.query(NoResult).limit(1000)
    indels = session.query(NoResult).filter(
        NoResult.seq_id == 'M03592:1:000000000-ADANF:1:2119:10104:6789'
    )
    '''
    indels = session.query(Sequence).filter(
        Sequence.probable_indel_or_misalign == 1
    ).limit(100)
    noresults = session.query(NoResult).limit(100)
    total = indels.count()
    print 'Creating task queue for {} indels'.format(total)
    for i, seq in enumerate(itertools.chain(indels, noresults)):
        if seq.sample_id not in mutation_cache:
            mutation_cache[seq.sample_id] = session.query(
                Sample.v_ties_mutations,
                Sample.v_ties_len
            ).filter(Sequence.sample_id == seq.sample_id).first()
        avg_mut, avg_len = mutation_cache[seq.sample_id]
        session.expunge(seq)
        tasks.add_task({
            'num': i,
            'sample_id': seq.sample_id,
            'seq_id': seq.seq_id,
            'seq': seq.sequence.replace('-', '').strip('N'),
            'avg_mut': avg_mut,
            'avg_len': avg_len
        })

    workers = min(args.nproc, tasks.num_tasks)
    for i in range(0, workers):
        tasks.add_worker(LocalAlignmentWorker(
            v_germlines,
            j_germlines,
            args.align_path,
            args.max_deletions,
            args.max_insertions)
        )

    tasks.start()
