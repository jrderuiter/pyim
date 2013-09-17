__author__ = 'j.d.ruiter'

def print_alignment(aln):
    query = aln['query_seq']
    target = aln['target_seq']

    query_index = aln['query_start']
    target_index = aln['target_start']

    query_str = query[:query_index]
    target_str = ' ' * (query_index - aln['target_start']) + target[:aln['target_start']]
    match_str = ' ' * query_index

    for stat in aln['alignment']:

        if stat == 'M' or stat == 'S':
            query_str = query_str + query[query_index]
            query_index = query_index + 1

            target_str = target_str + target[target_index]
            target_index = target_index + 1
        elif stat == 'D':
            query_str = query_str + query[query_index]
            query_index = query_index + 1

            target_str = target_str + '-'
        elif stat == 'I':
            query_str = query_str + '-'

            target_str = target_str + target[target_index]
            target_index = target_index + 1
        else:
            raise ValueError('Unknown character in alignment string')

        match_str = match_str + ('|' if stat == 'M' else ' ')

    query_str = query_str + query[query_index:]
    target_str = target_str + target[target_index:]

    print "Read:\t%s\tScore:\t%d" % (aln['query_name'], aln['score'])
    print "Target:\t%s\t\tIden:\t%3.2f%%" % (aln['target_name'], aln['identity'])

    print query_str
    print match_str
    print target_str
