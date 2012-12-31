import sys
# other system-wide
sys.path.append('venv/lib/python2.7/site-packages')

import termcolor
import datetime

from calibre.utils.config import JSONConfig


class Const:
    NCBI_EMAIL       = 'email_for_NCBI'
    CROSSREF_API_KEY = 'CROSSREF_API_KEY'


prefs = JSONConfig('plugins/CitationGetter')
prefs.defaults = {
    Const.NCBI_EMAIL: '',
    Const.CROSSREF_API_KEY: '',
    }



DLOG_OUT_FILE = "/tmp/cglog.txt"
DLOG_OUT_FILE = None

def dlog(*s):
    str_out = "\n".join([" %s %s %s\n" % (datetime.datetime.now(), termcolor.colored(">>>", "green", None, ["bold"]), line) for line in s]) + "\n"
    if DLOG_OUT_FILE:
        with open(DLOG_OUT_FILE, 'a') as ofile:
            ofile.write(str_out)
    else:
        sys.stderr.write(str_out)
        sys.stderr.flush()
