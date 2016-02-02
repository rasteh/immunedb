from sqlalchemy.sql import exists

import dnautils
from sldb.aggregation.tasks import collapse_subject_sequences
from sldb.common.celery_app import get_result
import sldb.common.config as config
from sldb.common.models import (Clone, Sample, Sequence, SequenceCollapse,
                                Subject)
import sldb.common.modification_log as mod_log


'''
class CollapseWorker(concurrent.Worker):
    """A worker for collapsing sequences without including positions where
    either sequences has an 'N'.
    :param Session session: The database session
    """
    def __init__(self, session):
        self._session = session
        self._tasks = 0

    def do_task(self, bucket_hash):
        seqs = self._session.query(
            Sequence.sample_id, Sequence.ai, Sequence.seq_id,
            Sequence.sequence, Sequence.copy_number
        ).filter(
            Sequence.bucket_hash == bucket_hash
        ).all()

        to_process = sorted([{
            'sample_id': s.sample_id,
            'ai': s.ai,
            'seq_id': s.seq_id,
            'sequence': s.sequence,
            'cn': s.copy_number
        } for s in seqs], key=lambda e: -e['cn'])

        while len(to_process) > 0:
            # Get the largest sequence in the list
            larger = to_process.pop(0)
            # Iterate over all smaller sequences to find matches
            instances = 1
            for i in reversed(range(0, len(to_process))):
                smaller = to_process[i]
                if dnautils.equal(larger['sequence'], smaller['sequence']):
                    # Add the smaller sequence's copy number to the larger
                    larger['cn'] += smaller['cn']
                    # If the smaller sequence matches the larger, collapse it
                    # to the larger
                    self._session.add(SequenceCollapse(**{
                        'sample_id': smaller['sample_id'],
                        'seq_ai': smaller['ai'],
                        'copy_number_in_subject': 0,
                        'collapse_to_subject_seq_ai': larger['ai'],
                        'collapse_to_subject_sample_id': larger['sample_id'],
                        'collapse_to_subject_seq_id': larger['seq_id'],
                        'instances_in_subject': 0
                    }))
                    instances += 1
                    # Delete the smaller sequence from the list to process
                    # since it's been collapsed
                    del to_process[i]

            # Update the larger sequence's copy number and "collapse" to itself
            self._session.add(SequenceCollapse(**{
                'sample_id': larger['sample_id'],
                'seq_ai': larger['ai'],
                'copy_number_in_subject': larger['cn'],
                'collapse_to_subject_sample_id': larger['sample_id'],
                'collapse_to_subject_seq_id': larger['seq_id'],
                'collapse_to_subject_seq_ai': larger['ai'],
                'instances_in_subject': instances,
            }))

        self._session.commit()
        self._tasks += 1
        if self._tasks > 0 and self._tasks % 100 == 0:
            self._print('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self._print('Committing collapsed sequences')
        self._session.commit()
        self._session.close()


'''
def record_collapse(collapsed_seqs):
    print collapsed_seqs

def run_collapse(session, args):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info=vars(args))
    subject_ids = []

    for subject in (args.subject_ids or map(
                lambda e: e.id, session.query(Subject.id).all()
                )):
        if session.query(Sample).filter(
                Sample.subject_id == subject,
                ~exists().where(
                    SequenceCollapse.sample_id == Sample.id
                )).first() is None:
            print 'Subject {} already collapsed.  Skipping.'.format(subject)
        else:
            print 'Resetting collapse info for subject {}'.format(subject)
            samples = session.query(Sample).filter(
                  Sample.subject_id == subject
            ).all()
            for sample in samples:
                session.query(SequenceCollapse).filter(
                    SequenceCollapse.sample_id == sample.id
                ).delete(synchronize_session=False)
            print 'Resetting clone info for subject {}'.format(subject)
            session.query(Clone).filter(Clone.subject_id == subject).delete()
            subject_ids.append(subject)
    session.commit()

    print 'Creating task queue to collapse {} subjects.'.format(
        len(subject_ids))

    for subject_id in subject_ids:
        buckets = session.query(
            Sequence.bucket_hash
        ).filter(
            Sequence.subject_id == subject_id
        ).group_by(
            Sequence.bucket_hash
        )
        tasks = []
        for bucket in buckets:
            sequences = session.query(
                Sequence.sample_id,
                Sequence.ai,
                Sequence.seq_id,
                Sequence.sequence,
                Sequence.copy_number
            ).filter(
                Sequence.bucket_hash == bucket.bucket_hash
            ).all()
            tasks.append(collapse_subject_sequences.delay(sequences))

        for task in tasks:
            record_collapse(get_result(task))
