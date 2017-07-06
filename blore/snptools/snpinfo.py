import collections
import numpy as np
import re

SNPINFO_FIELDS = ['rsid', 'bp_location', 'ref_allele', 'alt_allele']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

    #@property
    #def numid(self):
    #    return int(re.sub("^.*[sv-]", "", self.rsid))

    #def __repr__(self):
    #    parent_string = super(SnpInfo, self).__repr__()
    #    return '%s [numid=%i]' %(parent_string, self.numid)
