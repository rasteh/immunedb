from sqlalchemy.sql.expression import cast
from sqlalchemy.types import Float, Integer


def get_regions(insertions):
    regions = [78, 36, 51, 30, 114]
    offset = 0
    if insertions is not None and len(insertions) > 0:
        for i, region_size in enumerate(regions):
            for (start, size) in insertions:
                if offset <= start < offset + region_size:
                    regions[i] += size
            offset += region_size
    return regions


def get_pos_region(regions, cdr3_len, pos):
    cdr3_start = sum(regions)
    j_start = cdr3_start + cdr3_len
    if pos >= j_start:
        return 'FR4'
    elif pos >= cdr3_start:
        return 'CDR3'

    total = 0
    for i, length in enumerate(regions):
        total += length
        if pos < total:
            rtype = 'FW' if i % 2 == 0 else 'CDR'
            rnum = (i // 2) + 1
            return '{}{}'.format(rtype, rnum)


def ord_to_quality(quality):
    if quality is None:
        return None
    return ''.join([' ' if q is None else chr(q + 33) for q in quality])


def periodic_commit(session, query, interval=100):
    for i, r in enumerate(query):
        if i > 0 and i % interval == 0:
            session.commit()
        yield r
    session.commit()


def get_or_create(session, model, **kwargs):
    """Gets or creates a record based on some kwargs search parameters"""
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance, False
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance, True


def find_streak_position(s1, s2, max_streak):
    '''Finds the first streak of max_streak characters where s1 does not equal
    s2

    For example if max_streak is 3:
        ATCGATCGATCGATCG
        ATCGATCGATCTTACG
                     ^--- Returned index
    '''
    streak = 0
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        streak = streak + 1 if c1 != c2 else 0
        if streak >= max_streak:
            return i
    return None


def format_ties(ties, name, strip_alleles=False):
    if ties is None:
        return None

    ties = [e.replace(name, '') for t in ties for e in t.split('|')]
    if strip_alleles:
        ties = [e.split('*', 1)[0] for e in ties]
    return '{}{}'.format(name, '|'.join(sorted(set(ties))))


def int_cast(agg):
    return cast(agg, Integer)
