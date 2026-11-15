import json,sys
f=sys.argv[1]; t0=sys.argv[2]; t1=sys.argv[3]
for line in open(f):
    try: d=json.loads(line)
    except: continue
    if d.get('type')!='assistant': continue
    ts=d.get('timestamp','')
    if not (t0<=ts<=t1): continue
    c=d.get('message',{}).get('content')
    if not isinstance(c,list): continue
    t=' '.join(x.get('text','') for x in c if isinstance(x,dict) and x.get('type')=='text')
    if len(t.strip())<200: continue
    print('====',ts[:16]); print(t[:4000])
