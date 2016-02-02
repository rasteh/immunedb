import csv
import re

from sldb.identification import (add_as_noresult, add_as_sequence, add_uniques,
                                 AlignmentException)
from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import VGermlines
from sldb.common.models import NoResult, Sample, Study, Subject
import sldb.util.funcs as funcs


class ImportException(Exception):
    pass


DEFAULT_MAPPINGS = {
    'full_sequence': 'V-D-J-REGION',
    'seq_id': 'Sequence ID',
    'v_gene': 'V-GENE and allele',
    'j_gene': 'J-GENE and allele',
    'copy_number': None
}


def _collapse_seqs(session, sample, reader, columns):
    seqs = {}
    for record in reader:
        if record[columns.full_sequence] is None:
            NoResult(seq_id=record[columns.seq_id], sample_id=sample.id)
            continue

        record[columns.full_sequence] = record[columns.full_sequence].replace(
            '.', '').upper()
        if record[columns.full_sequence] not in seqs:
            seqs[record[columns.full_sequence]] = {
                'record': record,
                'seq_ids': []
            }
        seqs[record[columns.full_sequence]]['seq_ids'].append(
            record[columns.seq_id]
        )
        if columns.copy_number in record:
            for dup in record[columns.copy_number]:
                seqs[record[columns.full_sequence]]['seq_ids'].append(
                    '{}_DUP_{}'.format(columns.seq_id, dup)
                )

    return seqs.values()


def read_file(session, handle, sample, v_germlines, j_germlines, columns):
    seqs = _collapse_seqs(session, sample, csv.DictReader(handle,
                          delimiter='\t'), columns)

    aligned_seqs = {}
    missed = 0
    total = 0
    for total, seq in enumerate(seqs):
        if total > 0 and total % 1000 == 0:
            print 'Finished {}'.format(total)
            session.commit()
        v_genes = set(
            re.findall('IGHV[^ ,]+', seq['record'][columns.v_gene])
        )
        j_genes = set(
            re.findall('IGHJ[^ ,]+', seq['record'][columns.j_gene])
        )
        v_genes = filter(lambda v: v in v_germlines, v_genes)
        j_genes = filter(lambda j: j in j_germlines, j_genes)

        vdj = VDJSequence(
            seq['seq_ids'], seq['record'][columns.full_sequence], v_germlines,
            j_germlines, force_vs=v_genes, force_js=j_genes
        )
        try:
            if len(v_genes) == 0 or len(j_genes) == 0:
                raise AlignmentException('No V or J gene in input')
            vdj.analyze()
            vdj.align_to_germline()
            if vdj.sequence in aligned_seqs:
                aligned_seqs[vdj.sequence].ids += vdj.ids
            else:
                aligned_seqs[vdj.sequence] = vdj
        except AlignmentException as e:
            add_as_noresult(session, vdj, sample)
            missed += 1
    print 'Aligned {} / {} sequences'.format(total - missed + 1, total)

    print 'Collapsing ambiguous character sequences'
    add_uniques(session, sample, aligned_seqs.values())
    session.commit()


def run_import(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    study, new = funcs.get_or_create(session, Study, name=args.study_name)

    if new:
        print 'Created new study "{}"'.format(study.name)
        session.commit()

    sample, new = funcs.get_or_create(session, Sample, name=args.sample_name,
                                      study=study)
    if new:
        sample.date = args.date
        print 'Created new sample "{}"'.format(sample.name)
        for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                    'ig_class', 'v_primer', 'j_primer'):
            setattr(sample, key, vars(args).get(key, None))
        subject, new = funcs.get_or_create(
            session, Subject, study=study,
            identifier=args.subject)
        sample.subject = subject
        session.commit()
    else:
        print 'Sample already exists'
        return

    with open(args.input_file) as fh:
        read_file(session, fh, sample, v_germlines, j_germlines, args)
