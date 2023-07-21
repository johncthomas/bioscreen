
from bioscreen._imports import *

import attrs

logger = logging.getLogger("bioscreen.classes")
logging.basicConfig()


StrOptional = Optional[str]
CollectionOptional = Optional[Collection]

ComparisonsResults = pd.DataFrame

SampleGroup = str

class SigCols(enum.Enum):
    fdr = AMap(dict(key='FDR', short='FDR', long='FDR'))
    fdr10 = AMap(dict(key='FDR10', short='-log10(FDR)', long='-log10(FDR)'))

@attrs.define
class StatCol:
    key:str
    short:str=None
    long:str=None
    _formatted:str=None
    original:str = None

    @property
    def formatted(self):
        if self._formatted is None:
            return self.long
        return self._formatted


class StatColumns(AMap):
    def __init__(self, columns=Collection[StatCol]):
        cols = {col.key:col for col in columns}
        super().__init__(cols)

    @property
    def short(self) -> AMap:
        return AMap({k:col.short for k, col in self.items()})

    @property
    def long(self) -> AMap:
        return AMap({k: col.long for k, col in self.items()})

    @property
    def original(self) -> AMap:
        return AMap({k: col.original for k, col in self.items()})

    def original_to_good(self) -> dict:
        return {col.original:col.key for col in self.values()}


