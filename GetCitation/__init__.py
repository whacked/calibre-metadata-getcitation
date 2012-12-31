import sys, os
import urllib2, re, time, datetime


import glob
sys.path.append(os.path.expanduser('~/dev/calibre/calibre-src/src'))
sys.path.append(glob.glob("/usr/local/lib/python2.7/dist-packages/biopy*.egg")[0])
from Bio import Entrez

# for PyQt
sys.path.append('/usr/lib/python2.7/dist-packages')
from PyQt4 import QtGui, QtCore
from PyQt4.Qt import *

# other system-wide
sys.path.append('venv/lib/python2.7/site-packages')
import termcolor


from calibre.ebooks.metadata.sources.base import Source
from calibre.ebooks.metadata.book.base import Metadata


# import xmltodict
from lxml import etree

from exampledata import QueryResponse


USE_TEST_DATA = False
DLOG_OUT_FILE = "/tmp/cglog.txt"
DLOG_OUT_FILE = None

def dlog(*s):
    str_out = "\n".join([" %s %s %s\n" % (str(datetime.datetime.now()), termcolor.colored(">>>", "green", None, ["bold"]), line) for line in s]) + "\n"
    if DLOG_OUT_FILE:
        with open(DLOG_OUT_FILE, 'a') as ofile:
            ofile.write(str_out)
    else:
        sys.stderr.write(str_out)
        sys.stderr.flush()

def calibre_dict():
    dc = {'authors': []}
    for key in ("title",
                "pubdate",
                "author_sort",
                "publisher",
                "issue",
                "abstract",
                "keywords",
                "volume",
                "pages",
                ):
        dc[key] = None
    return dc

def mapget(el, *argv):
    return map(el.find, argv)

class Const:
    NCBI_EMAIL       = 'email_for_NCBI'
    CROSSREF_API_KEY = 'CROSSREF_API_KEY'
class Validator:
    @staticmethod
    def email(s):
        return re.match(r'[^@]+@[^\.]+\.\w+', s) is not None

class CGetter:
    entry_identifier = None
    ls_required_pref_key = []
    def does_understand(self, s):
        raise Exception("You should not call me directly")
    def __init__(self):
        raise Exception("You should not directly instantiate the CGetter class")
    def resolve(self, identifier):
        raise Exception("You should not call me directly")

class PubMedGetter(CGetter):
    entry_identifier = "pmid"
    ls_required_pref_key = [(Const.NCBI_EMAIL, Validator.email),]
    
    def does_understand(self, s):
        return str.isdigit(s)
    def __init__(self):
        print "PMGetter"

    def resolve(self, identifier):
        if USE_TEST_DATA:
            xml = QueryResponse['pubmed.xml']
        else:
            Entrez.email = prefs[Const.NCBI_EMAIL]
            handle = Entrez.efetch(db='pubmed', id=str(identifier), retmode="xml", tool="BioPython, Bio.Entrez, calibre-ebook")
            xml = handle.read()
        # dc = xmltodict.parse(xml)
        dc = calibre_dict()
        dom = etree.XML(xml)

        for key, searchpath in (("title"    , '//ArticleTitle'),
                                ("publisher" , '//Journal/Title'),
                                ("issue"     , '//JournalIssue/Issue'),
                                ("abstract"  , '//AbstractText'),
                                ("keywords"  , '//Keyword'),
                                ("volume"    , '//JournalIssue/Volume'),
                                ("pages"     , '//MedlinePgn'),
                                ):
            el_match = dom.xpath(searchpath)
            if el_match:
                dc[key] = el_match[0].text

            dlog(key, searchpath)

        # authors is special
        for author in dom.xpath('//Author'):
            # lastname = el_author.find('LastName')
            # forename = el_author.find('ForeName')
            # initials = el_author.find('Initials')
            lastname, forename, initials = mapget(author, 'LastName', 'ForeName', 'Initials')
            dc['authors'].append(("%s %s" % (forename != None and forename.text,
                                             lastname != None and lastname.text)).strip())

        el_pubdate = dom.xpath('//PubDate')
        if len(el_pubdate):
            _y, _m, _d = mapget(el_pubdate[0], "Year", "Month", "Day")
            m = d = 1
            if _y is not None:
                y = int(_y.text)
                if _m is not None:
                    if _m.text.isdigit():
                        m = int(_m.text)
                    elif len(_m.text) is 3:
                        m = time.strptime(_m.text, "%b").tm_mon
                    if _d is not None:
                        d = int(_d.text)
                dc["pubdate"] = datetime.datetime(y, m, d)
                del y
            del _y,_m,_d, m,d
        return dc

class CrossRefGetter(CGetter):

    entry_identifier = "doi"
    ls_required_pref_key = [(Const.CROSSREF_API_KEY, Validator.email),]

    CROSSREF_TEMPLATE = "http://www.crossref.org/openurl/?id=doi:%s&noredirect=true&pid=%s&format=unixref"
    
    def does_understand(self, s):
        # ref: http://stackoverflow.com/questions/27910/finding-a-doi-in-a-document-or-page
        return re.match(r'\s*(10[.][0-9]{3,}(?:[.][0-9]+)*/(?:(?!["&\'<>])\S)+)\b', s)
    def __init__(self):
        dlog(self.__class__)

    def resolve(self, identifier):
        if USE_TEST_DATA:
            xml = QueryResponse['crossref.doi']
        else:
            CROSSREF_LOOKUP = self.CROSSREF_TEMPLATE % (identifier, prefs[Const.CROSSREF_API_KEY])
            req = urllib2.Request(url = CROSSREF_LOOKUP)
            res = urllib2.urlopen(req)
            xml = res.read()

        dc = calibre_dict()
        dom = etree.XML(xml)
        
        for key, searchpath in (("title"    , '//journal_article/titles/title'),
                                ("publisher" , '//journal//full_title'),
                                ("issue"     , '//journal/journal_issue/issue'),
                                ("volume"    , '//journal/journal_issue/journal_volume/volume'),
                                ):
            el_match = dom.xpath(searchpath)
            if el_match:
                dc[key] = el_match[0].text

        el_pages = dom.xpath('//journal/journal_article/pages')
        if len(el_pages):
            first_page, last_page = map(lambda el: el.text, mapget(el_pages[0], "first_page", "last_page"))
            dc['pages'] = "%s-%s" % (first_page, last_page)

        el_pubdate = dom.xpath('//journal//publication_date')
        if len(el_pubdate):
            _y, _m, _d = mapget(el_pubdate[0], "year", "month", "day")
            m = d = 1
            if _y is not None:
                y = int(_y.text)
                if _m is not None:
                    if _m.text.isdigit():
                        m = int(_m.text)
                    elif len(_m.text) is 3:
                        m = time.strptime(_m.text, "%b").tm_mon
                    if _d is not None:
                        d = int(_d.text)
                dc["pubdate"] = datetime.datetime(y, m, d)
                del y
            del _y,_m,_d, m,d
        dlog(dc)

        #books.author_sort
        for person_name in dom.xpath("//journal/journal_article/contributors/person_name"):
            given_name, surname = mapget(person_name, "given_name", "surname")
            dc['authors'].append(("%s %s" % (given_name != None and given_name.text,
                                             surname != None and surname.text)).strip())
        return dc        




from calibre.utils.config import JSONConfig
from PyQt4.Qt import QWidget, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit

prefs = JSONConfig('plugins/CitationGetter')
prefs.defaults = {
    Const.NCBI_EMAIL: '',
    Const.CROSSREF_API_KEY: '',
    }

class ConfigWidget(QWidget):

    config_help_message = 'Set your API keys and access information for online querying services'
    
    def __init__(self, plugin):
        super(ConfigWidget, self).__init__()
        
        self.vlayout = QVBoxLayout()
        self.setLayout(self.vlayout)

        self.dcinput = {}

        for key in prefs.defaults:
            hlayout = QHBoxLayout()
            txtlabel = QLabel(key.replace("_", " "))
            txtinput = QLineEdit()
            txtinput.setText(str(prefs[key]))
            txtlabel.setBuddy(txtinput)
            self.dcinput[key] = txtinput
            
            hlayout.addWidget(txtlabel)
            hlayout.addWidget(txtinput)
            self.vlayout.addLayout(hlayout)

    def save_settings(self):
        for key in prefs.defaults:
            prefs[key] = unicode(self.dcinput[key].text())



class GetCitation(Source):

    name                    = 'Citation Getter'
    description             = 'Query supported sources for citation information'

    action_type             = 'current'
    
    supported_platforms     = ['windows', 'osx', 'linux']
    author                  = 'Sir Skeleton'
    version                 = (0, 0, 1)
    minimum_calibre_version = (0, 8, 0)

    capabilities            = frozenset(['identify'])
    touched_fields          = frozenset(['title', 'authors', 'publisher', 'pubdate', 'series'])

    tp_Getter = (
        PubMedGetter(),
        CrossRefGetter(),
        )

    def config_widget(self):
        return ConfigWidget(self)

    def save_settings(self, config_widget):
        '''
        Save the settings specified by the user with config_widget.
        
        :param config_widget: The widget returned by :meth:`config_widget`.
        '''
        config_widget.save_settings()

    def identify(self, 
                 log, result_queue, abort, 
                 title=None, authors=None,
                 identifiers={}, timeout=30):

        dlog("IDENTIFIERS")
        dlog(identifiers)

        dc = {}
        for Getter in self.tp_Getter:
            dlog("USING GETTER: %s" % Getter)
            dlog("my prefs:", prefs)
            isvalid = True
            for pref_key, is_valid_check in Getter.ls_required_pref_key:
                dlog("validating: %s --> = %s" % (pref_key, prefs[pref_key]))
                if not is_valid_check(prefs[pref_key]):
                    dlog("INVALID: %s" % prefs[pref_key])
                    isvalid = False
                    break
            entry_id = identifiers.get(Getter.entry_identifier, None)
            dlog("---------------", entry_id)
            if entry_id and Getter.does_understand(entry_id):
                dc = Getter.resolve(entry_id)
                dlog("RESULT")
                dlog(dc)
                if dc['authors']: break

        dlog(dc)
        if not dc:
            dlog("not dc: returning")
            return None
        dlog("OK")
        
        mi = Metadata(dc["title"], dc["authors"])
        for attr in ("pubdate",
                     "publisher",
                     "issue",
                     "abstract",
                     "keywords",
                     "volume",
                     "pages",
                     ):
            if hasattr(mi, attr) and attr in dc:
                setattr(mi, attr, dc[attr])

        mi.print_all_attributes()
        self.clean_downloaded_metadata(mi)
        result_queue.put(mi)
        dlog("=" * 20)
        dlog(dir(self))
        dlog("=" * 20)

        return None



if __name__ == "__main__":
    from Queue import Queue
    from threading import Event
    import pprint

    USE_TEST_DATA = True
    if USE_TEST_DATA:
        print termcolor.colored("    >>> USING TEST DATA <<<", "yellow", None, ["bold"])

    identifiers = {'pmid': "15798944",
                   'doi': '10.1055/s-2005-867080'}
    if False:
        dc = {}
        for Getter in GetCitation.tp_Getter:
            entry_id = identifiers.get(Getter.entry_identifier, None)
            if entry_id and Getter.does_understand(entry_id):
                dlog('USING', Getter)
                dc = Getter.resolve(entry_id)
                break
        pprint.pprint(dc)
        sys.exit()
    
    app = QtGui.QApplication([])
    # ConfigWidget(None)
    cg = GetCitation(None)
    dc = cg.identify(None, Queue(), Event(), identifiers = identifiers)
    

