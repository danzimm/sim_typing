
models = []
customTemplates = []

def type_for_name(name):
  return models[name]['type']

def model_names():
  return models.keys()

def specials():
  return customTemplates

SNANAidxTable = {
  104: "CSP-2006ep",
  103: "CSP-2004gv",
  105: "CSP-2007Y",
  202: "SDSS-000020",
  203: "SDSS-002744",
  212: "s11-2006jo",
  234: "SDSS-019323",
  21: "SNLS-04D11a",
  218: "SDSS-017548",
  22: "SNLS-04D4jv",
  101: "CSP-2004fe",
  102: "CSP-2004gg",
  205: "SDSS-004012",
  207: "SDSS-013195",
  211: "SDSS-014475",
  217: "SDSS-015475",
  201: "s11-2004hx",
  204: "SDSS-003818",
  208: "SDSS-013376",
  210: "SDSS-014450",
  213: "SDSS-014599",
  214: "SDSS-015031",
  215: "SDSS-015320",
  216: "SDSS-015339",
  219: "SDSS-017564",
  220: "SDSS-017862",
  221: "SDSS-018109",
  222: "SDSS-018297",
  223: "SDSS-018408",
  224: "SDSS-018441",
  225: "SDSS-018457",
  226: "SDSS-018590",
  227: "SDSS-018596",
  228: "SDSS-018700",
  229: "SDSS-018713",
  230: "SDSS-018734",
  231: "SDSS-018793",
  232: "SDSS-018834",
  233: "SDSS-018892",
  235: "SDSS-020038",
  206: "SDSS-012842",
  209: "SDSS-013449"
}

def SNANAidx_to_model(idx):
  if idx in SNANAidxTable:
    return SNANAidxTable[idx]
  return ""

def modelmap_initialize_modeldefs(tag=None):
  global models, customTemplates
  modeldefsfile = 'modeldefs' + ('' if tag is None else '_' + tag)
  modeldefs = __import__(modeldefsfile, globals(), locals(), ["models", "customTemplates"], -1)
  models = modeldefs.models
  customTemplates = modeldefs.customTemplates

