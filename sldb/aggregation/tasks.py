import dnautils

from sldb.common.celery_app import app

@app.task
def collapse_subject_sequences(sequences):
    sequences = sorted(sequences, cmp=lambda a, b:
                       cmp(a.copy_number, b.copy_number))
    current_largest = 0
    finished = []
    while current_largest < len(sequences):
        larger = sequences[current_largest]
        instances = 1
        for i in reversed(range(len(sequences[current_largest + 1:]))):
            smaller = sequences[i]
            if dnautils.equal(larger.sequence, smaller.sequence):
                larger.ds += smaller.ids
                # Record the smaller collapsed sequence
                finished.append({
                    'sample_id': smaller.sample_id,
                    'seq_ai': smaller.ai,
                    'copy_number_in_subject': 0,
                    'collapse_to_subject_seq_ai': larger.ai,
                    'collapse_to_subject_sample_id': larger.sample_id,
                    'collapse_to_subject_seq_id': larger.seq_id,
                    'instances_in_subject': 0
                })
                instances += 1
                del sequences[i]

        finished.append({
            'sample_id': larger.sample_id,
            'seq_ai': larger.ai,
            'copy_number_in_subject': larger.copy_number,
            'collapse_to_subject_sample_id': larger.sample_id,
            'collapse_to_subject_seq_id': larger.seq_id,
            'collapse_to_subject_seq_ai': larger.ai,
            'instances_in_subject': instances,
        })

        current_largest += 1

    return sequences
