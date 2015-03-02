from __future__ import absolute_import
import logging
import math
from subprocess import PIPE

from pynpact.steps import extract, nprofile, acgt_gamma
from pynpact import binfile
from pynpact import capproc, parsing
from pynpact.util import Hasher, reducedict, which, replace_ext
from pynpact.steps import producer, enqueue


log = logging.getLogger('pynpact.steps.allplots')
statuslog = logging.getLogger('pynpact.statuslog')


BIN = binfile('Allplots')

KEYS = ['first_page_title', 'following_page_title',
        'length', 'startBase', 'endBase', 'period',
        'basesPerGraph', 'graphsPerPage', 'x-tics',
        'nucleotides', 'alternate_colors', 'basename']
FILE_KEYS = ['File_of_unbiased_CDSs',
             'File_of_conserved_CDSs',
             'File_of_new_CDSs',
             'File_of_published_rejected_CDSs',
             'File_of_stretches_where_CG_is_asymmetric',
             'File_of_published_accepted_CDSs',
             'File_of_potential_new_CDSs',
             'File_of_blocks_from_new_ORFs_as_cds',
             'File_of_blocks_from_annotated_genes_as_cds',
             'File_of_GeneMark_regions',
             'File_of_G+C_coding_potential_regions',
             'File_of_met_positions (e.g.:D 432)',
             'File_of_stop_positions (e.g.:D 432)',
             'File_of_tatabox_positions (e.g.:105.73 D 432 TATAAAAG)',
             'File_of_capbox_positions',
             'File_of_ccaatbox_positions',
             'File_of_gcbox_positions',
             'File_of_kozak_positions',
             'File_of_palindrom_positions_and_size',
             'File_list_of_nucleotides_in_200bp windows',
             'File_list_of_nucleotides_in_100bp windows']


def plan(config, executor):
    # unless extract is already disabled.
    if 'allplots' in config:
        return
    config['allplots'] = True
    if which('ps2pdf'):
        return convert_ps_to_pdf(config, executor)
    else:
        return combine_ps_files(config, executor)


def allplots(config, executor):
    after = []
    try:
        after.extend(extract.plan(config, executor))
    except:
        pass
    after.extend(nprofile.plan(config, executor))
    after.extend(acgt_gamma.plan(config, executor))

    parsing.length(config)
    parsing.first_page_title(config)
    parsing.following_page_title(config)
    parsing.endBase(config)

    h = Hasher()
    # Strip down to the config for this task only
    rconfig = reducedict(config, KEYS + FILE_KEYS)

    basesPerGraph = rconfig['basesPerGraph']
    graphsPerPage = rconfig['graphsPerPage']
    startBase = rconfig.pop('startBase')
    endBase = rconfig.pop('endBase')
    bp_per_page = rconfig['bp_per_page'] = basesPerGraph * graphsPerPage
    page_count = math.ceil(float(endBase - startBase) / bp_per_page)
    log.info("Generating %d pages of allplots", page_count)
    page_num = 1  # page number offset
    filenames = []
    waiton = []
    # per-page loop
    while startBase < endBase:
        pconfig = dict(rconfig.items())
        pconfig['page_num'] = page_num
        pconfig['startBase'] = startBase
        if startBase + bp_per_page < endBase:
            pconfig['endBase'] = startBase + bp_per_page
        else:
            pconfig['endBase'] = endBase
        h = Hasher().hashdict(pconfig).hashfiletime(BIN).hashfiletime(__file__)
        psname = parsing.derive_filename(config, h.hexdigest(), 'ps')
        filenames.append(psname)
        waiton.extend(enqueue(_ap, executor, pconfig, psname, after=after))
        page_num += 1
        startBase += bp_per_page

    # Finally set the output filenames into the master config dict
    config['psnames'] = filenames
    return waiton


@producer()
def _ap(pconfig, out):
    page_num = pconfig['page_num']
    cmd = [BIN, "-q", "--stdin"]
    if pconfig.get('alternate_colors'):
        cmd.append("-C")

    # add the rest of the required args
    cmd += [pconfig['startBase'],
            pconfig['bp_per_page'],
            pconfig['graphsPerPage'],
            pconfig['x-tics'],
            pconfig['period'],
            pconfig['endBase']]

    with capproc.guardPopen(
            cmd, stdin=PIPE, stdout=out, stderr=False,
            logger=log) as ap:
        # write the allplots.def file information through stdin.
        apdef = build_allplots_def(pconfig, page_num)
        log.debug("Writing Allplots.def:\n%s", apdef)
        ap.stdin.write(apdef)
        ap.stdin.close()
        ap.wait()
        log.debug("Finished allplots with rc %s", ap.returncode)


def build_allplots_def(pconfig, page_num):
    """Write the configuration file for Allplots that it's expecting on stdin

    Allplots used to check for 'Allplots.def' in the dcurrent
    directory but we modified it to be able to read this configuration
    on stdin, this writes in that format.

    """
    parsing.first_page_title(pconfig)
    parsing.following_page_title(pconfig)
    lines = []
    wl = lines.append

    # NB: the "Plot Title" is disregarded, but that line
    # should also contain the total number of bases
    wl("%s %d" % ("FOOBAR", pconfig['length']))

    # Nucleotide(s)_plotted (e.g.: `% CG`)
    wl(''.join(pconfig['nucleotides']))
    # First-Page title
    wl(pconfig['first_page_title'].format(page_num))
    # Title of following pages
    wl(pconfig['following_page_title'].format(page_num))
    # write the rest of the filenames that allplots might read from
    for k in FILE_KEYS:
        # allplots doesn't like blank lines so make sure there is at
        # least a dummy value on every line.
        wl(pconfig.get(k, "None"))
    return '\n'.join(lines)


def combine_ps_files(config, executor):
    after = allplots(config, executor)
    psnames = config['psnames']
    log.debug("Going to combine %d postscript files", len(psnames))
    combined_ps_name = parsing.derive_filename(
        config, Hasher().hashlist(psnames).hexdigest(), 'ps')
    config['combined_ps_name'] = combined_ps_name
    return enqueue(
        _combine_ps_files, executor, config, combined_ps_name, after=after)


@producer()
def _combine_ps_files(config, psout):
    psnames = config['psnames']
    # While combining, insert the special markers so that
    # it will appear correctly as many pages.
    statuslog.info("Combining pages into output.")
    log.debug("Combining %d ps files", len(psnames))

    psout.write("%!PS-Adobe-2.0\n\n")
    psout.write("%%Pages: {0}\n\n".format(len(psnames)))
    idx = 1
    for psf in psnames:
        psout.write("%%Page: {0}\n".format(idx))
        with open(psf, 'r') as infile:
            infile.readline()
            infile.readline()
            psout.write(infile.readline())
            for l in infile:
                psout.write(l)
        idx += 1


def convert_ps_to_pdf(config, executor):
    after = combine_ps_files(config, executor)
    combined_ps_name = config['combined_ps_name']
    pdf_filename = replace_ext(combined_ps_name, 'pdf')
    config['pdf_filename'] = pdf_filename
    return enqueue(_ps2pdf, executor, config, pdf_filename, after=after)


@producer()
def _ps2pdf(config, out):
    ps_filename = config['combined_ps_name']
    statuslog.info("Converting to PDF")
    cmd = ['ps2pdf', ps_filename, '-']
    capproc.capturedCall(
        cmd, stdout=out, logger=log, check=True)
    statuslog.info("Finished PDF conversion")
