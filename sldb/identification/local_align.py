import multiprocessing as mp
import re
import subprocess
import shlex
from sqlalchemy import func

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
    def __init__(self, local_align_path, max_deletions, max_insertions):
        self._local_align_path = local_align_path
        self._max_deletions = max_deletions
        self._max_insertions = max_insertions

    def do_task(self, args):
        j_align = self._align_seq_to_germs(
            args['seq_id'], args['seq'], args['j_ties']
        )

        if j_align is None:
            print args['seq_id'], 'Bad J'
            return
        v_seq = args['seq'][:j_align['seq_pos']]
        v_align = self._align_seq_to_germs(args['seq_id'], v_seq,
                                           args['v_ties'])
        if v_align is None:
            print args['seq_id'], 'Bad V'
            return

        missing_start = v_align['seq_pos'] + v_align['seq_len']
        missing_len = (j_align['seq_pos'] - v_align['seq_pos']
            - len(v_align['seq']))
        missing = args['seq'][missing_start:missing_start + missing_len]

        new_seq = ''.join((
            v_align['seq'],
            missing,
            j_align['seq']
        ))

        v_stripped = v_align['germ'].replace('-', '')
        missing = v_align['germ_orig'][
            v_align['germ_pos'] + len(v_stripped):
            v_align['germ_pos'] + len(v_stripped) + missing_len
        ].lower()
        new_germ = ''.join((
            v_align['germ'],
            missing,
            j_align['germ']
        ))

        if v_align['germ_pos'] > 0:
            new_germ = v_align['germ_orig'][:v_align['germ_pos']] + new_germ
            diff = len(new_germ) - len(new_seq)
            if diff > 0:
                new_seq = ''.join((
                    v_align['seq_orig'][:v_align['seq_pos']][-diff:],
                    new_seq
                ))
        #TODO: Add back germline J if sequence is partial

	gapped_germ = args['v_ties'][tuple(v_align['germ_name'].split('|'))]
        germ_gaps = gap_positions(v_align['germ'])
        imgt_gaps = [
            i + gaps_before(germ_gaps, i)
            for i, c in enumerate(gapped_germ) if c == '-'
        ]
        for gap in imgt_gaps:
            new_germ = new_germ[:gap] + '-' + new_germ[gap:]
            new_seq = new_seq[:gap] + '-' + new_seq[gap:]

        cdr_start = 309 + v_align['germ'].count('-')
        cdr_end = (
            len(new_seq) - args['j_ties'].upstream_of_cdr3 -
            j_align['germ'].count('-')
        )

    def _align_seq_to_germs(self, seq_id, seq, germs):
        stdin = []
        for g_name, g_seq in germs.iteritems():
            if isinstance(g_name, tuple):
                g_name = '|'.join(g_name)
            stdin.append('>{}\n{}\n'.format(g_name, g_seq.replace('-', '')))
            stdin.append('>{}\n{}\n'.format(seq_id, seq))

        cmd = (
            '{} --match 4 --mismatch -2 --gapopen -5 --gapextend -2 '
            '--wildcard N 2 --maxhits 1 --file - '
            '--printfasta --printseq'
        ).format(self._local_align_path)
        proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        output, err = proc.communicate(''.join(stdin))
        regex = (
            r'==[^\n]+\n'
            r'(?P<germ_name>.+)\n'
            r'(?P<germ_orig>.+)\n'
            r'(?P<seq_name>.+)\n'
            r'(?P<seq_orig>.+)\n\n'
            r'hit .* (?P<score>.+)\n'
            r'\s*(?P<germ>[^\s]+)\s*'
            r'\[pos: (?P<germ_pos>\d+); len: (?P<germ_len>\d+)\]\n'
            r'\s*(?P<seq>[^\s]+)\s*'
            r'\[pos: (?P<seq_pos>\d+); len: (?P<seq_len>\d+)\]'
        )

        best = None
        for match in re.finditer(regex, output, re.MULTILINE):
            match = match.groupdict()
            match['score'] = int(match['score'])
            match['germ_pos'] = int(match['germ_pos'])
            match['germ_len'] = int(match['germ_len'])
            match['seq_pos'] = int(match['seq_pos'])
            match['seq_len'] = int(match['seq_len'])
            if best is None or match['score'] > best['score']:
                best = match
        return best


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3, 0, 0)

    mutation_cache = {}
    tasks = concurrent.TaskQueue()

    #indels = session.query(NoResult).limit(100)
    indels = session.query(NoResult).filter(
        NoResult.seq_id == 'M03592:1:000000000-ADANF:1:2119:10104:6789'
    )
    total = indels.count()
    fixed = 0
    print 'Creating task queue for {} indels'.format(total)
    for i, seq in enumerate(indels):
        if seq.sample_id not in mutation_cache:
            mutation_cache[seq.sample_id] = session.query(
                Sample.v_ties_mutations,
                Sample.v_ties_len
            ).filter(Sequence.sample_id == seq.sample_id).first()
        avg_mut, avg_len = mutation_cache[seq.sample_id]
        session.expunge(seq)
        tasks.add_task({
            'sample_id': seq.sample_id,
            'seq_id': seq.seq_id,
            'seq': seq.sequence,
            'v_ties': v_germlines.all_ties(avg_len, avg_mut, cutoff=False),
            'j_ties': j_germlines,
            'avg_mut': avg_mut,
            'avg_len': avg_len
        })


    workers = min(args.nproc, tasks.num_tasks)
    for i in range(0, workers):
        tasks.add_worker(LocalAlignmentWorker(
            args.local_align_path,
            args.max_deletions,
            args.max_insertions)
        )

    tasks.start()
