"""Lower indexed serializers into structural reader layouts.

Unsupported control flow is retained as an explicit error instruction, never
silently dropped. This module runs at scheme generation time, not in readers.
"""
import json
import re


def generate(root, texts, indexed, contracts, ops, balanced, split_args, declarations, infer_type):
    layouts, parents, records = {}, {}, {}
    def make_record(body, decls, stream):
        return extract_record(body, decls, stream, balanced, split_args, infer_type)
    def normalize(owner):
        owner = owner.split(' void ')[-1]
        return re.sub(r'\s+', '', owner)
    for key, value in indexed.items():
        records[normalize(key)] = value
    # Capture inline serializers with the declarations of their enclosing class.
    for path, text in texts.items():
        classes = []
        for m in re.finditer(r'\bclass\s+(?:CASADI_EXPORT\s+)?(\w+)(?:\s*:[^{;]+)?\s*\{', text):
            end = balanced(text, m.end()-1)
            classes.append((m.start(), end, m[1], text[m.end():end]))
            parent = re.search(r':\s*public\s+([\w:<>, ]+)', m[0])
            if parent:
                parents[m[1]] = normalize(parent[1].rstrip('{ '))
        for m in re.finditer(r'\bvoid\s+(serialize_node|serialize_type|serialize_body|serialize)\s*\(SerializingStream\s*&\s*(\w+)\)\s*const\s*(?:override\s*)?\{', text):
            enclosing = [c for c in classes if c[0] < m.start() < c[1]]
            if not enclosing: continue
            cls = max(enclosing, key=lambda c:c[0])
            body = text[m.end():balanced(text,m.end()-1)].strip()
            records[cls[2]+'::'+m[1]] = make_record(body, declarations(cls[3]), m[2])
        # Helpers which serialize data on behalf of a class.
        for m in re.finditer(r'static\s+void\s+(pack_\w+)\s*\(SerializingStream\s*&\s*(\w+)[^{]+\{', text):
            body=text[m.end():balanced(text,m.end()-1)].strip()
            records[m[1]]=make_record(body,declarations(body),m[2])

    def resolve_type(field):
        t=field['cpp_type']
        if t is None:
            types=contracts.get(field.get('name'),{}).get('cpp_types',[])
            if len(types)==1: t=next(iter(types))
        return t

    def compile_record(record):
        body=record['cpp_body']
        events=[]
        for f in record['fields']:
            if f['operation']=='version':
                step={'op':'version','name':f['name'],'value':int(f['expression'])}
            else:
                step={'op':'field','name':f['name'],'type':resolve_type(f),'bind':f['expression']}
                if f.get('name_expression'): step['name_expression']=f['name_expression']
            events.append((f['offset'],step))
        for c in record.get('calls',[]):
            events.append((c['offset'],{'op':'call','layout':normalize(c['owner']+'::'+c['method'])}))
        for m in re.finditer(r'\b(pack_\w+)\s*\(s\s*,',body):
            end=balanced(body,m.end()-1 if body[m.end()-1]=='(' else body.index('(',m.start()),'(',')')
            args=split_args(body[body.index('(',m.start())+1:end])
            events.append((m.start(),{'op':'call','layout':m[1],'params':{'d':json.loads(args[1])}}))
        for m in re.finditer(r'\b[\w:]+(?:\.\w+|->\w+)?\s*\([^;]*?\bs\b[^;]*?\)\s*;', body):
            if not any(body.rfind(';', 0, m.start())+1 <= position < m.end() for position, _ in events):
                events.append((m.start(), {'op':'unsupported','reason':'unlowered serialization call: '+m[0]}))
        ranges=[]
        for m in re.finditer(r'\b(if|for|while|switch)\s*\(',body):
            end=balanced(body,m.end()-1,'(',')'); condition=body[m.end():end].strip()
            start=end+1
            while start<len(body) and body[start].isspace():start+=1
            if start<len(body) and body[start]=='{': finish=balanced(body,start)
            else: finish=body.find(';',start)
            if finish<0:continue
            ranges.append((m.start(),start,finish,m[1],condition))
        def region(start,end):
            children=[r for r in ranges if start<=r[0]<end and not any(start<=p[0]<r[0]<p[2]<=end for p in ranges)]
            ordered=[(p,s) for p,s in events if start<=p<end and not any(r[0]<=p<=r[2] for r in children)]
            for p,a,b,kind,condition in children:
                inner=region(a+1,b)
                if not inner: continue  # Computation without serialized output.
                if kind=='if': step={'op':'if','condition':condition,'body':inner}
                elif kind=='for' and ':' in condition:
                    step={'op':'repeat','count':condition.split(':',1)[1].strip()+'.size()','body':inner}
                else:step={'op':'unsupported','reason':kind+' ('+condition+')'}
                ordered.append((p,step))
            return [s for _,s in sorted(ordered,key=lambda x:x[0])]
        result=region(0,len(body))
        # An unhandled else branch must not be read unconditionally.
        if re.search(r'\belse\b',body):
            return [{'op':'unsupported','reason':'unlowered else branch'}]
        return result

    for key, record in records.items():
        layouts[key]=compile_record(record)
    functions={}
    text=texts[root/'casadi/core/function_internal.cpp']
    for name,cls in re.findall(r'\{"([^"]+)",\s*(\w+)::deserialize\}',text):functions[name]=cls
    mx={str(ops[op]) if op in ops else op:normalize(cls) for op,cls in re.findall(r'\{(OP_\w+|-1),\s*([\w<>:, ]+)::deserialize\}',texts[root/'casadi/core/mx_node.cpp'])}
    # The native dispatcher uses these byte-layout families for arithmetic nodes.
    math=texts[root/'casadi/core/calculus.hpp']
    def family(kind):
        macro=re.search(r'#define CASADI_MATH_'+kind+r'_BUILTIN[^\n]*(?:\n[^\n]*\\)*\n[^\n]*',math)[0]
        method=re.search(r'bool casadi_math<T>::is_'+kind.lower()+r'\([^)]*\)\s*\{',math)
        body=math[method.end():balanced(math,method.end()-1)]
        return set(re.findall(r'case (OP_\w+)',macro+body))
    unary,binary=family('UNARY'),family('BINARY')
    sx={str(ops['OP_PARAMETER']):'SymbolicSX',str(ops['OP_CALL']):'CallSX','-1':'OutputSX'}
    for name in unary:
        mx[str(ops[name])]='UnaryMX';sx[str(ops[name])]='UnarySX'
    for name in binary:
        mx[str(ops[name])]='BinaryMX';sx[str(ops[name])]='BinarySX'
    generic_text=texts[root/'casadi/core/generic_type.hpp']
    names=re.findall(r'\bOT_\w+',generic_text[generic_text.index('enum TypeID'):generic_text.index('};',generic_text.index('enum TypeID'))])
    generic={}
    for m in re.finditer(r'typedef\s+GenericTypeInternal<(OT_\w+),\s*([\s\S]*?)>\s+\w+Type;',texts[root/'casadi/core/generic_type.cpp']):
        generic[str(names.index(m[1]))]=normalize(m[2])
    plugins = {}
    for text in texts.values():
        creator = re.search(r'plugin->creator\s*=\s*(\w+)::creator', text)
        name = re.search(r'plugin->name\s*=\s*"([^"]+)"', text)
        if creator and name:
            cls = creator[1]
            base = cls
            seen = set()
            while base not in seen:
                seen.add(base)
                for function, owner in functions.items():
                    if owner == base:
                        plugins.setdefault(function, {})[name[1]] = cls
                base = parents.get(base, '').split('<')[0]
                if not base: break
    programs = type_programs(layouts, parents, functions, mx, sx, generic, ops, plugins)
    stream = texts[root/'casadi/core/serializing_stream.cpp']
    for match in re.finditer(r'void SerializingStream::pack\(const (\w+)& \w+\)\s*\{', stream):
        body=stream[match.end():balanced(stream,match.end()-1)]
        decoration=re.search(r"decorate\('([^']+)'\)",body)
        if match[1] in programs:
            programs[match[1]]['shared']='shared_pack' in body
            if decoration:programs[match[1]]['decoration']=decoration[1]
    return {'file_types': {'0':'Sparsity','2':'DM','4':'Linsol','5':'Function','6':'GenericType','7':'casadi_int','8':'double','9':'std::string','10':'std::vector<Sparsity>','12':'std::vector<DM>','15':'std::vector<Function>','16':'std::vector<GenericType>','17':'std::vector<casadi_int>','18':'std::vector<double>','19':'std::vector<std::string>'}, 'types': programs, 'version':1,'layouts':dict(sorted(layouts.items())), 'parents':dict(sorted(parents.items())),
            'dispatch':{'Function':functions,'MX':mx,'SXElem':sx,'GenericType':generic}}



def extract_record(body, decls, stream, balanced, split_args, infer_type):
    fields=[];calls=[]
    for m in re.finditer(r'\b'+stream+r'\.(pack|version)\s*\(',body):
        end=balanced(body,m.end()-1,'(',')');args=split_args(body[m.end():end])
        if len(args)!=2:continue
        name=json.loads(args[0]) if re.fullmatch(r'"[^"\\]*"',args[0]) else None
        fields.append({'operation':m[1],'name':name,'name_expression':args[0],
                       'expression':args[1],'cpp_type':infer_type(args[1],decls),'offset':m.start()})
    for m in re.finditer(r'([\w:<>, ]+)::(serialize_body|serialize_type|delayed_serialize_members)\s*\('+stream+r'\)',body):
        calls.append({'owner':m[1].strip(),'method':m[2],'offset':m.start()})
    return {'cpp_body':body,'fields':fields,'calls':calls}


def type_programs(layouts, parents, functions, mx, sx, generic, ops, plugins):
    def field(name, t, bind=None):
        return {'op':'field','name':name,'type':t,'bind':bind or name}
    def call(name, **params):
        return {'op':'call','layout':name, **({'params':params} if params else {})}
    def select(bind, cases):
        return {'op':'select','bind':bind,'cases':cases}
    def body(cls):
        # Concrete template substitutions remain layout metadata, not runtime logic.
        if cls in ('MXFunction','SXFunction'):
            mat='MX' if cls=='MXFunction' else 'SX'
            return [call(cls+'::serialize_body',MatType=mat)]
        return [call(cls+'::serialize_body')]
    function_cases={name:body(cls) for name,cls in functions.items()}
    for name, cls in functions.items():
        prefix = layouts.get(cls+'::serialize_type', [])
        prefix = [s for s in prefix if s != {'op':'call','layout':'FunctionInternal::serialize_type'}]
        if prefix and all(s['op']=='field' for s in prefix):
            selector = next((s for s in prefix if s.get('bind')=='class_name()'), None)
            if selector:
                candidates={cls}
                for child in parents:
                    parent=child;seen=set()
                    while parent and parent not in seen:
                        seen.add(parent)
                        if parent==cls:candidates.add(child)
                        parent=parents.get(parent, '').split('<')[0]
                function_cases[name]=prefix+[select('class_name()', {c:body(c) for c in sorted(candidates)})]
            else:
                function_cases[name]=prefix+body(cls)

    for family, registrations in plugins.items():
        function_cases[family]=[field('PluginInterface::plugin_name','std::string','plugin'),
                                select('plugin',{name:body(cls) for name,cls in registrations.items()})]
    types={
        'Function':{'shared':True,'body':[field('Function::null','bool','null'),
            {'op':'if','condition':'!null','body':[
                field('FunctionInternal::base_function','std::string','base'),
                select('base',function_cases)]}]},
        'Sparsity':{'shared':True,'body':[field('SparsityInternal::compressed','std::vector<casadi_int>')]},
        'Slice':{'body':[call('Slice::serialize')]},
        'DM':{'body':[field('Matrix::sparsity','Sparsity'),field('Matrix::nonzeros','std::vector<double>')]},
        'SX':{'body':[field('Matrix::sparsity','Sparsity'),field('Matrix::nonzeros','std::vector<SXElem>')]},
        'GenericType':{'body':[field('GenericType::type','int','type'),select('type',{
            k:[field('GenericType::d',v)] for k,v in generic.items()})]},
    }
    # Expand dispatch families to their serialized layouts, without evaluating nodes.
    mx_cases={k:body(v) for k,v in mx.items()}
    for k,v in mx.items():
        if v=='BinaryMX':
            mx_cases[k]=[field('BinaryMX::scalar_flags','char')]+body('BinaryMX<ScX,ScY>')
    for name,stem in [('OP_GETNONZEROS','GetNonzeros'),('OP_SETNONZEROS','SetNonzeros'),('OP_ADDNONZEROS','SetNonzeros')]:
        suffix='<Add>' if stem=='SetNonzeros' else ''
        mx_cases[str(ops[name])]=[field(stem+'::type','char','subtype'),select('subtype',{
            str(ord(c)):body(stem+variant+suffix) for c,variant in [('a','Vector'),('b','Slice'),('c','Slice2')]
        })]
    constants={str(ord(c)):[call('MXNode::serialize_body')] for c in ['0','1','-']}
    for char,typ in [('D','double'),('I','casadi_int')]:
        constants[str(ord(char))]=[field('Constant::value',typ),call('MXNode::serialize_body')]
    for char,cls in [('a','ConstantDM'),('f','ConstantFile'),('p','ConstantPool')]:
        constants[str(ord(char))]=body(cls)
    mx_cases[str(ops['OP_CONST'])]=[field('ConstantMX::type','char','subtype'),select('subtype',constants)]
    sx_cases={k:[call(v+'::serialize_node')] for k,v in sx.items()}
    constant_cases={}
    for cls in ['RealtypeSX','IntegerSX','ZeroSX','OneSX','MinusOneSX','InfSX','MinusInfSX','NanSX']:
        key=cls+'::serialize_node'
        if key not in layouts: continue
        prefix=layouts[key][0]
        literal=prefix.get('bind','')
        if len(literal)==3 and literal[0]=="'":
            constant_cases[str(ord(literal[1]))]=layouts[key][1:]
    sx_cases[str(ops['OP_CONST'])]=[field('ConstantSX::type','char','subtype'),select('subtype',constant_cases)]
    types['MX']={'shared':True,'body':[field('MXNode::op','int','op'),select('op',mx_cases)]}
    types['SXElem']={'shared':True,'body':[field('SXNode::op','casadi_int','op'),select('op',sx_cases)]}
    resources={cls:body(cls) for cls in ['ZipMemResource','ZipResource','DirResource']}
    types['Resource']={'shared':True,'body':[{'op':'version','name':'ResourceInternal','value':1},
        field('ResourceInternal::type','std::string','resource_type'),select('resource_type',resources)]}
    for name, tag in {'Function':'F','Sparsity':'S','Slice':'S','SXElem':'E','MX':'X','GenericType':'G','Resource':'R'}.items():
        types[name]['decoration'] = tag
    return types
