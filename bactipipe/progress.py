"""Small optional execution events, not estimates of elapsed-time progress."""
import json
import os
import time

PREFIX = 'BACTIPIPE_PROGRESS '
STAGES = ('inputs', 'read_qc', 'assembly', 'identification', 'assembly_qc', 'publication')


def emit(stage, state='running', *, sample=None, total=None, finished=None, failed=None, metrics=None):
    if os.environ.get('BACTIPIPE_PROGRESS_EVENTS') != '1':
        return
    if stage not in STAGES or state not in ('running', 'completed', 'failed', 'skipped'):
        return
    event = dict(schema=1, stage=stage, state=state, timestamp=time.time(), monotonic=time.monotonic())
    if metrics is not None:
        event['metrics'] = metrics
    if sample is not None:
        event['sample'] = str(sample)[:200]
    for key, value in (('total', total), ('finished', finished), ('failed', failed)):
        if value is not None:
            event[key] = int(value)
    try:
        os.write(1, (PREFIX + json.dumps(event, ensure_ascii=True) + '\n').encode())
    except (OSError, ValueError):
        pass
