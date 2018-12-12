#!/usr/bin/python
import sys
import os
import re
from subprocess import call
import ipdb
import re

def strip_page_group(file_name_in, file_name_out):
    # adopted from here:
    # https://tex.stackexchange.com/questions/76273/multiple-pdfs-with-page-group-included-in-a-single-page-warning
    file_in = open(file_name_in, 'r')
    file_out = open(file_name_out, 'w')

    page_group = None
    for line in file_in:
        # print("%s"% line)
        if page_group is None:
            if line.rstrip() == b"  /Group <<": #line.endswith(b"/Group <<"):
                #ipdb.set_trace()
                print('FOUND')
                page_group = [line]
            else:
                file_out.write(line)
                #stdout.write(line)
        else:
            page_group.append(line)
            if line.rstrip() == b"  >>":
                break
    else: # This is only executed if page group is NOT found.
        if page_group:
            stdout.write(b"".join(page_group))
            page_group = None
    for line in file_in: # copy the remainder to the out file
        file_out.write(line)

    file_in.close()
    file_out.close()

    if page_group:
        print(b"".join(page_group))
        return 0
    else:
        print(b"note: did not find page group\n")
        return 1



def main(arg_s):

    # parse the arguments to find the output file name
    out_file = None
    for arg in arg_s:
        if arg.startswith('--file'):
            print("in file %s"% arg.split('=')[1])
            in_file = arg.split('=')[1]
        elif arg.startswith('--export-pdf'):
            print("out file %s"% arg.split('=')[1])
            out_file = arg.split('=')[1]

    if out_file is None:
        ipdb.set_trace()
    CMD = ['inkscape'] + arg_s
    devnull = open(os.devnull, 'w')
    call(CMD, stdout=devnull, stderr=devnull)

    # Original pipeline:
    # qpdf --qdf input.pdf - | python strip_page_group.py | fix-qdf >output.pdf
    call(['qpdf', '--qdf', out_file, 'tmp1.qdf'])


    tmp_file1 = 'tmp1.qdf'
    tmp_file2 = 'tmp2.qdf'

    strip_page_group(tmp_file1, tmp_file2)
    of = open(out_file, 'w')
    call(['fix-qdf', tmp_file2], stdout=of)
    of.close()

if __name__ == '__main__':
    main(sys.argv[1:])
    # argv should look like
    # x = ['/home/arnold/matlab/afm_mpc_journal/latex/temp/svg2pdf.py',
    #'-z', '-D', '--export-latex', '--file=lqr_locus_constsig_0p9.svg',
    # '--export-pdf=TEST_out.pdf']
