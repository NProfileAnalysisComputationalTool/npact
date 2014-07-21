defaults = {
        'nucleotides': ['c', 'g'],

        # keys for extract.c
        'GeneDescriptorKey1': 'gene',
        'GeneDescriptorKey2': 'locus_tag',
        'GeneDescriptorSkip1': 0,
        'GeneDescriptorSkip2': 0,

        # keys for CG:
        'window_size': 201,
        'step': 51,
        'period': 3,

        # acgt_gamma:
        'skip_prediction': False,
        'significance': "0.001",

        ##allplots
        #http://docs.python.org/library/string.html#formatspec
        'first_page_title': None,
        'following_page_title': 'Page {0}',

        'bp_per_page': 50000,
        'start_base': 0,

        }


def initial(filename):
    config = defaults.copy()
    config['filename'] = filename
    config['basename'] = os.path.basename(filename)
    return config


def length(config):
    if length not in config:
        config['length'] = 0
    return config
