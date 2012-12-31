import sys, os

# for PyQt
sys.path.append('/usr/lib/python2.7/dist-packages')
from PyQt4 import QtGui, QtCore
from PyQt4.Qt import *


from calibre.ebooks.metadata.sources.base import Source
from calibre.ebooks.metadata.book.base import Metadata

from common import *

from CGetter import *


from PyQt4.Qt import QWidget, QHBoxLayout, QVBoxLayout, QLabel, QLineEdit

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
        for Getter in CGetter.getter_list:
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
    import termcolor

    CGetter.USE_TEST_DATA = True
    if CGetter.USE_TEST_DATA:
        print termcolor.colored("    >>> USING TEST DATA <<<", "yellow", None, ["bold"])

    identifiers = {'pmid': "15798944",
                   'doi': '10.1055/s-2005-867080'}
    if True:
        dc = {}
        for Getter in CGetter.getter_list:
            entry_id = identifiers.get(Getter.entry_identifier, None)
            if entry_id and Getter.does_understand(entry_id):
                dlog('USING', Getter)
                dc = Getter.resolve(entry_id)
            pprint.pprint(dc)
    else:

        app = QtGui.QApplication([])
        # ConfigWidget(None)
        cg = GetCitation(None)
        dc = cg.identify(None, Queue(), Event(), identifiers = identifiers)
