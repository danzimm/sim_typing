
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
      'CSP-2006ep': 'SN Ib'
  }[name]

def model_names():
  return [
      'salt2-extended',
      's11-2004hx',
      's11-2005lc',
      's11-2005hl',
      's11-2005hm',
      's11-2005gi',
      's11-2006fo',
      's11-2006jo',
      's11-2006jl',
      'CSP-2006ep'
  ]

def specials():
  return [
      'CSP-2006ep'
  ]
