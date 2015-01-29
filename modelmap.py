
def type_for_name(name):
  return {
      'salt2-extended': 'SN Ia',
      's11-2004hx': 'SN IIL/P',
      's11-2005lc': 'SN IIP',
      's11-2005hl': 'SN Ib',
      's11-2005hm': 'SN Ib',
      's11-2005gi': 'SN IIP',
      's11-2006fo': 'SN Ic',
      's11-2006jo': 'SN Ib',
      's11-2006jl': 'SN IIP',
      'CSP-2006ep': 'SN Ib',
      'SDSS-017548': 'SN Ic',
      'SDSS-000018': 'SN IIP'
  }[name]

def model_names():
  return [
      's11-2004hx',
      's11-2006fo',
      's11-2005hl',
      'salt2-extended',
      's11-2005lc',
      's11-2005hm',
      's11-2005gi',
      's11-2006jo',
      's11-2006jl',
      'CSP-2006ep',
      'SDSS-017548',
      'SDSS-000018',
      'SDSS-004012',
      'CSP-2004fe',
      'SDSS-018457'
  ]

def specials():
  return model_names()[9:]

SNANAidxTable = {
  104: "CSP-2006ep",
  103: "CSP-2004gv",
  105: "CSP-2007Y",
  202: "SDSS-000020",
  203: "SDSS-002744",
  212: "SDSS-014492",
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
  201: "SDSS-000018",
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

