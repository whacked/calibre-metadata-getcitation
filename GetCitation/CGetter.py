import urllib2, re, time, datetime
from exampledata import QueryResponse
from common import *

# import xmltodict
from lxml import etree
import os, glob
sys.path.append(os.path.expanduser('~/dev/calibre/calibre-src/src'))
sys.path.append(glob.glob("/usr/local/lib/python2.7/dist-packages/biopy*.egg")[0])
from Bio import Entrez


class Validator:
    @staticmethod
    def email(s):
        return re.match(r'[^@]+@[^\.]+\.\w+', s) is not None

def mapget(el, *argv):
    return map(el.find, argv)
def mapget_text(el, *argv):
    return [got is not None and got.text or "" for got in mapget(el, *argv)]


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



class CGetter:
    getter_list = []

    USE_TEST_DATA = False

    entry_identifier = None
    ls_required_pref_key = []
    def does_understand(self, s):
        raise Exception("You should not call me directly")
    def __init__(self):
        dlog(self.__class__)
    def resolve(self, identifier):
        raise Exception("You should not call me directly")

    @classmethod
    def add(selfclass, getterclass):
        selfclass.getter_list.append(getterclass())

def process_pubdate(el_pubdate, *timekey):
    if len(el_pubdate):
        _year, _month, _day = mapget_text(el_pubdate[0], *timekey)
        month = day = 1
        if _year:
            if _month:
                if _month.isdigit():
                    month = int(_month)
                elif len(_month) is 3:
                    month = time.strptime(_month, "%b").tm_mon
                if _day:
                    day = int(_day)
            return datetime.datetime(int(_year), month, day)

@CGetter.add
class PubMedGetter(CGetter):
    entry_identifier = "pmid"
    ls_required_pref_key = [(Const.NCBI_EMAIL, Validator.email),]
    
    def does_understand(self, s):
        return str.isdigit(s)

    def resolve(self, identifier):
        if self.USE_TEST_DATA:
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

        # authors is special
        for author in dom.xpath('//Author'):
            # lastname = el_author.find('LastName')
            # forename = el_author.find('ForeName')
            # initials = el_author.find('Initials')
            lastname, forename, initials = mapget_text(author, 'LastName', 'ForeName', 'Initials')
            dc['authors'].append(("%s %s" % (forename, lastname)).strip())

        dc["pubdate"] = process_pubdate(dom.xpath('//PubDate'),
                                        "Year", "Month", "Day")
        return dc

# # @CGetter.add
# class PLoSGetter(CGetter):

#     entry_identifier = "doi"
#     ls_required_pref_key = [(Const.PLOS_API_KEY, lambda s: isinstance(s, basestring) and len(s) > 25)]
    
#     def does_understand(self, s):
#         if not re.match(r'\s*(10[.][0-9]{3,}(?:[.][0-9]+)*/(?:(?!["&\'<>])\S)+)\b', s):
#             return False
#         if s.split(".")[0].lower() in ("pone",):
#             return False
#         return True

#     def resolve(self, identifier):
#         pass

@CGetter.add
class CrossRefGetter(CGetter):

    entry_identifier = "doi"
    ls_required_pref_key = [(Const.CROSSREF_API_KEY, Validator.email),]

    CROSSREF_TEMPLATE = "http://www.crossref.org/openurl/?id=doi:%s&noredirect=true&pid=%s&format=unixref"
    
    def does_understand(self, s):
        # ref: http://stackoverflow.com/questions/27910/finding-a-doi-in-a-document-or-page
        return re.match(r'\s*(10[.][0-9]{3,}(?:[.][0-9]+)*/(?:(?!["&\'<>])\S)+)\b', s)

    def resolve(self, identifier):
        if self.USE_TEST_DATA:
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
            first_page, last_page = mapget_text(el_pages[0], "first_page", "last_page")
            dc['pages'] = "%s-%s" % (first_page, last_page)

        dc["pubdate"] = process_pubdate(dom.xpath('//journal//publication_date'),
                                        "year", "month", "day")
        
        #books.author_sort
        for person_name in dom.xpath("//journal/journal_article/contributors/person_name"):
            given_name, surname = mapget_text(person_name, "given_name", "surname")
            dc['authors'].append(("%s %s" % (given_name, surname)).strip())
        return dc        


