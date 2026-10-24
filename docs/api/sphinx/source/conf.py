# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'CasADi'
copyright = '2018, Joel Andersson, Joris Gillis, Greg Horn'
author = 'Joel Andersson, Joris Gillis, Greg Horn'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = '3.4.4'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosectionlabel', 
    #'sphinxcontrib.osexample',
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

add_module_names = False
toc_object_entries = False  # sidebar lists sections only
autodoc_typehints = "none"
autodoc_default_options = {'exclude-members': 'thisown'}

# SWIG proxies: keep documented members of the casadi bases, drop object/junk
def skip_member(app, what, name, obj, skip, options):
    if name in ('thisown', 'this') or name.startswith('_'):
        return True
    return skip

# Docstrings arrive from SWIG (signature block) + doxy2swigX (prose, provenance links).
import re
PROVENANCE = ('Doc source', 'Implementation', 'Extra doc')


def _split_top(text):
    """Split a See-also list on commas outside parentheses: 'a, Sparsity(n, m)'."""
    parts, depth, cur = [], 0, ''
    for ch in text:
        depth += {'(': 1, ')': -1}.get(ch, 0)
        if ch == ',' and depth == 0:
            parts.append(cur); cur = ''
        else:
            cur += ch
    return parts + [cur]

def _xref(target, context):
    """Link a See-also name to the class member or module function it denotes, if any."""
    import casadi
    target = re.sub(r'^casadi::', '', target).replace('::', '.')  # C++ scoped names, e.g. Function::transform
    bare = re.sub(r'\(.*\)$', '', target).strip()
    parts = context.split('.')
    cls = parts[2] if len(parts) > 3 else None  # casadi.casadi.<Class>.<member>
    if '.' not in bare and cls and hasattr(getattr(casadi, cls, None), bare):
        return ':py:obj:`%s <casadi.casadi.%s.%s>`' % (target, cls, bare)
    obj = casadi
    for p in bare.split('.'):
        obj = getattr(obj, p, None)
    if obj is not None:
        return ':py:obj:`%s <casadi.casadi.%s>`' % (target, bare)
    return '``%s``' % target

def _skip_first_literal(lines):
    """Lines without the first '::' literal block (SWIG's signature list)."""
    for i, l in enumerate(lines):
        if l.strip() == '::':
            ind = len(l) - len(l.lstrip()); j = i + 1
            while j < len(lines) and (not lines[j].strip() or len(lines[j]) - len(lines[j].lstrip()) > ind):
                j += 1
            return lines[:i] + lines[j:]
    return lines

def _fix_literals(line):
    """Repair inline literals produced from Doxygen \\p / backticks so RST accepts them."""
    if line.lstrip().startswith(('>>>', '>> ', '...', '|', '+')):  # examples and grid tables
        return line
    line = line.replace('````', '')               # \\p a-by-\\p b: adjacent spans merge
    parts = line.split('``')
    if len(parts) % 2 == 0:                        # unbalanced: leave untouched
        return line
    for k in range(0, len(parts), 2):              # plain text: RST metacharacters are literal
        t = re.sub(r'(?<!\\)\|', r'\\|', parts[k])                    # |x| is no substitution
        t = re.sub(r'(?<![\\*])\*(?!\*)', r'\\*', t)                  # x*y is no emphasis; **bold** stays
        parts[k] = re.sub(r'(?<=[A-Za-z0-9])_(?=[^\w`]|$)', r'\\_', t)  # name_ is no link
    for k in range(1, len(parts), 2):              # odd parts are literal contents
        lit = parts[k].strip()
        m = re.match(r'^(.*?)([.,;:]+)$', lit)
        if m and m.group(1):                       # trailing punctuation goes outside
            lit = m.group(1)
            parts[k + 1] = m.group(2) + parts[k + 1]
        parts[k] = lit
        if parts[k - 1] and re.search(r'\w$', parts[k - 1]):
            parts[k - 1] += ' '
        if parts[k + 1] and re.match(r'^[\w(]', parts[k + 1]):
            parts[k + 1] = ' ' + parts[k + 1]
    return '``'.join(parts)

def _fix_prose(lines):
    """Apply _fix_literals to prose only: skip literal blocks ('::') and directive bodies."""
    out, skip_ind = [], None
    for line in lines:
        ind = len(line) - len(line.lstrip())
        if skip_ind is not None:
            if not line.strip() or ind > skip_ind:
                out.append(line); continue
            skip_ind = None
        out.append(_fix_literals(line) if line.strip() else line)
        if line.rstrip().endswith('::') or line.lstrip().startswith('.. '):
            skip_ind = ind
    return out

def _merge_overload_groups(lines):
    """SWIG autodoc appends, per overload group, '.......' + '::' signature block + doc +
    '.............'. Drop empty groups; print a lone group's doc as is; label several."""
    start = [i for i, l in enumerate(lines) if re.match(r'^\s*\.{7}\s*$', l)]
    if not start:
        return lines
    head, groups, i = lines[:start[0]], [], start[0]
    while i < len(lines):
        if not re.match(r'^\s*\.{7}\s*$', lines[i]):
            i += 1
            continue
        j = i + 1
        while j < len(lines) and not re.match(r'^\s*\.{7,}\s*$', lines[j]):
            j += 1
        block, i = lines[i + 1:j], j + (1 if j < len(lines) and re.match(r'^\s*\.{13}\s*$', lines[j]) else 0)
        sigs, k = [], 0
        while k < len(block) and not block[k].strip(): k += 1
        if k < len(block) and block[k].strip() == '::':
            k += 1
            while k < len(block) and not block[k].strip(): k += 1
            while k < len(block) and block[k].strip():
                sigs.append(block[k].strip()); k += 1
        body = block[k:]
        if not any(l.strip() for l in body):
            continue
        # Doxygen copies a shared doc onto every overload; fold identical copies (they differ
        # only in their provenance lines) into one group listing all their signatures
        key = '\n'.join(l.strip() for l in body if l.strip() and not re.match(
            r'^(Doc source|Implementation|Extra doc):?\s*(https?://\S+)?$|^https?://\S+$', l.strip()))
        for g in groups:
            if g[2] == key:
                g[0].extend(x for x in sigs if x not in g[0])
                break
        else:
            groups.append((list(sigs), body, key))
    # SWIG also prints the (first) doc in the head: a group repeating the head adds nothing
    head_key = '\n'.join(l.strip() for l in _skip_first_literal(head) if l.strip() and not re.match(
        r'^(Doc source|Implementation|Extra doc):?\s*(https?://\S+)?$|^https?://\S+$', l.strip()))
    groups = [g for g in groups if g[2] != head_key]
    out = list(head)
    for sigs, body, _ in groups:
        if len(groups) > 1 and sigs:
            out += ['', '.. rst-class:: overload', '',
                    '**Overload%s:** %s' % ('s' if len(sigs) > 1 else '',
                                           ', '.join('``%s``' % s for s in sigs)), '']
        out += body
    return out


import json, os
_SIDECAR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_doctest_sidecar.json')
try:
    _SIDECAR = json.load(open(_SIDECAR_PATH))
except (OSError, ValueError):
    _SIDECAR = {}
_LEXER = {'matlab': 'matlabsession', 'julia': 'jlcon', 'javascript': 'javascript'}

def _tag_prose(out, prose, pad, cls):
    """Wrap the paragraphs in prose that end out (up to blank lines) in a container tagged cls."""
    end = len(out)
    while end and not out[end - 1].strip():
        end -= 1
    start = end
    while start:
        k = start
        while k and out[k - 1].strip():
            k -= 1
        if not k < start or ' '.join(l.strip() for l in out[k:start]) not in prose:
            break
        start = k
        while start and not out[start - 1].strip():
            start -= 1
    while start < end and not out[start].strip():
        start += 1
    if start < end:
        body = [pad + '   ' + l[len(pad):] if l.strip() else '' for l in out[start:end]]
        out[start:] = [pad + '.. container:: ' + cls, ''] + body + ['']

def _tag_examples(lines):
    """Turn '::' blocks of >>> examples into tagged pycon blocks; add their translations."""
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line.strip() != '::':
            out.append(line); i += 1; continue
        ind = len(line) - len(line.lstrip())
        j = i + 1
        while j < len(lines) and (not lines[j].strip() or len(lines[j]) - len(lines[j].lstrip()) > ind):
            j += 1
        body = lines[i + 1:j]
        content = [l for l in body if l.strip()]
        if not content or not content[0].lstrip().startswith('>>>'):
            out += lines[i:j]; i = j; continue
        key = '\n'.join(l.strip() for l in content)
        extra = _SIDECAR.get(key, [])
        prose = {' '.join(l.strip() for l in t.split('\n'))
                 for lang, parts in extra if lang == 'python' for _, t in parts}
        extra = [e for e in extra if e[0] != 'python']
        langs = {lang for lang, _ in extra}
        pad = ' ' * ind
        cls = 'doctest doctest-python' + ''.join(' has-' + l for l in sorted(langs))
        if prose:
            _tag_prose(out, prose, pad, cls)
        out += [pad + '.. code-block:: pycon', pad + '   :class: ' + cls, '']
        cind = min(len(l) - len(l.lstrip()) for l in content)
        out += [pad + '   ' + l[cind:] if l.strip() else '' for l in body]
        for lang, parts in extra:
            sub = pad
            if len(parts) > 1 or parts[0][0] != 'doctest':  # prose too: one tagged container
                out += ['', pad + '.. container:: doctest doctest-' + lang]
                sub = pad + '   '
            for kind, text in parts:
                if kind == 'text':
                    out += [''] + [sub + l for l in text.split('\n')]
                    continue
                out += ['', sub + '.. code-block:: ' + (_LEXER.get(lang, 'text') if kind == 'doctest' else 'text'),
                        sub + '   :class: doctest doctest-' + lang, '']
                out += [sub + '   ' + l if l.strip() else '' for l in text.split('\n')]
        out.append('')
        i = j
    return out

def process_docstring(app, what, name, obj, options, lines):
    lines[:] = _merge_overload_groups(lines)
    out, links = [], []
    first, sig_done = True, False
    for i, line in enumerate(lines):
        s = line.strip()
        if not s:
            out.append(line)
            continue
        if first:  # the brief is indented 2 deeper than the rest: would render as a block quote
            first = False
            if line.startswith(' '):
                out.append(s)
                continue
        if s == '[INTERNAL]':
            continue
        m = re.match(r'^(Julia|JavaScript|MATLAB) example:$', s)
        if m:  # doxy2swigX label for a non-Python \doctest block; render as a tagged code block
            lang = {'Julia': 'julia', 'JavaScript': 'javascript', 'MATLAB': 'matlab'}[m.group(1)]
            j = i + 1
            while j < len(lines) and not lines[j].strip(): j += 1
            if j < len(lines) and lines[j].strip() == '::':
                indent = lines[j][:len(lines[j]) - len(lines[j].lstrip())]
                lines[j] = indent + '.. code-block:: ' + lang
                lines.insert(j + 1, indent + '   :class: doctest doctest-' + lang)
            continue
        if s == '::':
            nxt = [l for l in lines[i+1:] if l.strip()]
            if not nxt or not nxt[0].startswith(' ' * (len(line) - len(line.lstrip()) + 1)):
                continue  # SWIG emits an empty literal block for classes without a signature
            for k in range(i + 1, len(lines)):  # SWIG output-argument placeholders
                if lines[k].strip() and not lines[k].startswith(line[:len(line) - len(line.lstrip())] + ' '):
                    break
                lines[k] = re.sub(r'\s*\bOUTPUT\d*\b', '', lines[k]).replace(' ,', ',')
            if not sig_done:  # the first literal block is SWIG's overload list
                sig_done = True
                out.append('.. rst-class:: signatures')
                out.append('')
        m = re.match(r'^See also:\s*(.*)$', s)
        if m:  # doxygen \sa list -> cross references, resolved on the live module
            indent = line[:len(line) - len(line.lstrip())]
            out.append(indent + '**See also:** ' + ', '.join(
                _xref(x.strip(), name) for x in _split_top(m.group(1)) if x.strip()))
            continue
        m = re.match(r'^(Doc source|Implementation|Extra doc): *(https?://\S+)?$', s)
        if m:
            links.append([m.group(1), m.group(2)])
            continue
        if links and links[-1][1] is None and re.match(r'^https?://\S+$', s):
            links[-1][1] = s
            continue
        out.append(line)
    if links:  # merged overload docstrings repeat the links: keep the first of each URL
        seen, uniq = set(), []
        for k, u in links:
            if u and k not in seen:
                seen.add(k); uniq.append((k, u))
        out += ['', '.. rst-class:: provenance', '',
                ' | '.join('`%s <%s>`__' % (k, u) for k, u in uniq), '']
    lines[:] = _fix_prose(_tag_examples(out))

def process_signature(app, what, name, obj, options, signature, return_annotation):
    if signature == '(*args)':  # SWIG wrapper; the real overloads are listed in the docstring
        return ('(...)', None)

def setup(app):
    app.connect('autodoc-skip-member', skip_member)
    app.connect('autodoc-process-docstring', process_docstring)
    app.connect('autodoc-process-signature', process_signature)


numfig = True
numfig_secnum_depth = (2)
numfig_format = {'figure': 'Figure %s', 'table': 'Table %s', 'code-block': 'Code Block %s', 'section': 'Section %s'}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#

source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'
smartquotes = False  # option names and strings keep their straight quotes


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'CasADidoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'CasADi.tex', 'CasADi Documentation',
     'Joel Andersson, Joris Gillis, Greg Horn', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'casadi', 'CasADi Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'CasADi', 'CasADi Documentation',
     author, 'CasADi', 'One line description of project.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
