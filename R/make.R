
from <- "../../metrumresearchgroup/devmodels/"
to <- file.path("inst", "models")
cp <- function(from_dir,from_file,to_dir,to_file=from_file) {
  from <- file.path(from,from_dir,from_file)
  to_dir <- file.path(to)
  to_file <- file.path(to_dir,to_file)
  file.copy(from,to_file,overwrite=TRUE)

}

cp("cipro", "cipro.cpp", "cipro")
cp("gcsf", "gcsf.cpp", "gcsf")
cp("opg", "opg.cpp", "opg")
cp("epo", "epo.cpp", "epo")
cp("azithro", "azithro.cpp", "azithro")
cp("sunitinib", "sunit.cpp", "sunit")
cp("moxi", "moxi.cpp", "moxi")
cp("secukinumab", "secukinumab.cpp", "secukinumab")
cp("pembrolizumab", "pembro.cpp", "pembro","pembro.cpp")
cp("meropenem", "meropenem.cpp", "meropenem")
cp("conway", "conway.cpp", "hiv2", "hiv2.cpp")
cp("ddi", "csa.cpp", "csa")
cp("ddi", "yoshikado.cpp", "pitavddi", "pitavddi.cpp")
cp("rifampicin", "rifampicin.cpp", "rifampicin")
cp("rifampicin", "rifampicin_midazolam.cpp", "rifampicin_midazolam")
cp("rifampicin", "midazolam.cpp", "midazolam")
cp("pbpk", "pbpk.txt", "pbpk", "pbpk.cpp")
