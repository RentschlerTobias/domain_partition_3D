import json,sys
f=sys.argv[1]
for line in open(f):
    try: d=json.loads(line)
    except: continue
    if d.get('type')!='user': continue
    m=d.get('message',{})
    c=m.get('content')
    if isinstance(c,str): t=c
    elif isinstance(c,list):
        t=' '.join(x.get('text','') for x in c if isinstance(x,dict) and x.get('type')=='text')
    else: continue
    if not t.strip() or t.startswith('<'): continue
    print('----',d.get('timestamp','')[:16]); print(t[:1500])
