#!/usr/bin/env python
"""Doxygen XML to SWIG docstring converter.

Usage:

  doxy2swig.py [options] input.xml output.i

Converts Doxygen generated XML files into a file containing docstrings
that can be used by SWIG-1.3.x.  Note that you need to get SWIG
version > 1.3.23 or use Robin Dunn's docstring patch to be able to use
the resulting output.

input.xml is your doxygen generated XML file and output.i is where the
output will be written (the file will be clobbered).

"""
#
#
# This code is implemented using Mark Pilgrim's code as a guideline:
#   http://www.faqs.org/docs/diveintopython/kgp_divein.html
#
# Author: Prabhu Ramachandran
# License: BSD style
#
# Thanks:
#   Johan Hake:  the include_function_definition feature
#   Bill Spotz:  bug reports and testing.
#   Sebastian Henschel:   Misc. enhancements.
#   Michael Thon:  Modified to my liking (June 2015)
#

from xml.dom import minidom
import re
import textwrap
import sys
import os.path
import optparse


def my_open_read(source):
    if hasattr(source, "read"):
        return source
    else:
        return open(source)


def my_open_write(dest):
    if hasattr(dest, "write"):
        return dest
    else:
        return open(dest, 'w')

class Doxy2SWIG:
    """Converts Doxygen generated XML files into a file containing
    docstrings that can be used by SWIG-1.3.x that have support for
    feature("docstring").  Once the data is parsed it is stored in
    self.pieces.

    """

    def __init__(self, src, include_function_definition=True, quiet=False):
        """Initialize the instance given a source object.  `src` can
        be a file or filename.  If you do not want to include function
        definitions from doxygen then set
        `include_function_definition` to `False`.  This is handy since
        this allows you to use the swig generated function definition
        using %feature("autodoc", [0,1]).

        """
        f = my_open_read(src)
        self.my_dir = os.path.dirname(f.name)
        self.xmldoc = minidom.parse(f).documentElement
        f.close()

        self.pieces = []
        self.pieces.append('\n// File: %s\n' %
                           os.path.basename(f.name))

        self.space_re = re.compile(r'\s+')
        self.lead_spc = re.compile(r'^(%feature\S+\s+\S+\s*?)"\s+(\S)')
        self.multi = 0
        self.ignores = ['inheritancegraph', 'param', 'listofallmembers',
                        'innerclass', 'name', 'declname', 'incdepgraph',
                        'invincdepgraph', 'programlisting', 'type',
                        'references', 'referencedby', 'location',
                        'collaborationgraph', 'reimplements',
                        'reimplementedby', 'derivedcompoundref',
                        'basecompoundref', 'argsstring']
        #self.generics = []
        self.include_function_definition = include_function_definition
        self.quiet = quiet
        self.textwidth = 80
        self.indent = 0
        self.listitem = ''

    def generate(self):
        """Parses the file set in the initialization.  The resulting
        data is stored in `self.pieces`.

        """
        self.parse(self.xmldoc)

    def parse(self, node):
        """Parse a given node.  This function in turn calls the
        `parse_<nodeType>` functions which handle the respective
        nodes.

        """
        pm = getattr(self, "parse_%s" % node.__class__.__name__)
        pm(node)

    def parse_Document(self, node):
        self.parse(node.documentElement)

    def get_Text(self, node):
        txt = node.firstChild.data
        txt = txt.replace('\\', r'\\')
        txt = txt.replace('"', r'\"')
        m = self.space_re.match(txt)
        if m and len(m.group()) == len(txt):
            return ''
        else:
            return txt

    def parse_Text(self, node):
        txt = node.data
        txt = txt.replace('\\', r'\\')
        txt = txt.replace('"', r'\"')
        # ignore pure whitespace
        m = self.space_re.match(txt)
        if m and len(m.group()) == len(txt):
            pass
        else:
            self.add_text(txt)

    def parse_Element(self, node):
        """Parse an `ELEMENT_NODE`.  This calls specific
        `do_<tagName>` handers for different elements.  If no handler
        is available the `generic_parse` method is called.  All
        tagNames specified in `self.ignores` are simply ignored.

        """
        name = node.tagName
        ignores = self.ignores
        if name in ignores:
            return
        attr = "do_%s" % name
        if hasattr(self, attr):
            handlerMethod = getattr(self, attr)
            handlerMethod(node)
        else:
            self.generic_parse(node)
            #if name not in self.generics: self.generics.append(name)

    def parse_Comment(self, node):
        """Parse a `COMMENT_NODE`.  This does nothing for now."""
        return

    def add_text(self, value):
        """Adds text corresponding to `value` into `self.pieces`."""
        if isinstance(value, (list, tuple)):
            self.pieces.extend(value)
        else:
            self.pieces.append(value)

    def get_specific_subnodes(self, node, name, recursive=0):
        """Given a node and a name, return a list of subnodes matching the
        name. Search recursively for `recursive` levels.
        """
        children = [x for x in node.childNodes if x.nodeType == x.ELEMENT_NODE]
        ret = [x for x in children if x.tagName == name]
        if recursive > 0:
            for x in children:
                ret.extend(self.get_specific_subnodes(x, name, recursive-1))
        return ret

    def get_specific_nodes(self, node, names):
        """Given a node and a sequence of strings in `names`, return a
        dictionary containing the names as keys and child
        `ELEMENT_NODEs`, that have a `tagName` equal to the name.

        """
        nodes = [(x.tagName, x) for x in node.childNodes
                 if x.nodeType == x.ELEMENT_NODE and
                 x.tagName in names]
        return dict(nodes)

    def generic_parse(self, node, pad=0):
        """A Generic parser for arbitrary tags in a node.

        Parameters:

         - node:  A node in the DOM.
         - pad: `int` (default: 0)

           If 0 the node data is not padded with newlines.  If 1 it
           appends a newline after parsing the childNodes.  If 2 it
           pads before and after the nodes are processed.  Defaults to
           0.

        """
        npiece = 0
        if pad:
            npiece = len(self.pieces)
            if pad == 2:
                self.add_text('\n')
        for n in node.childNodes:
            self.parse(n)
        if pad:
            if len(self.pieces) > npiece:
                self.add_text('\n')

    def get_type(self, node):
        pieces, self.pieces = self.pieces, []
        type = self.get_specific_subnodes(node, 'type')
        if type:
            self.generic_parse(type[0])
        type_str = ''.join(self.pieces)
        self.pieces = pieces
        return type_str

    def surround_parse(self, node, pre_char, post_char):
        self.add_text(pre_char)
        self.generic_parse(node)
        self.add_text(post_char)
    
    do_ref = generic_parse
    do_formula = generic_parse
    
    def do_linebreak(self, node):
        self.add_text('  ')
    
    def do_emphasis(self, node):
        self.surround_parse(node, '_', '_')

    def do_bold(self, node):
        self.surround_parse(node, '**', '**')
    
    def do_heading(self, node):
        pieces, self.pieces = self.pieces, []
        level = int(node.attributes['level'].value)
        self.generic_parse(node)
        if level == 1:
            # self.pieces.insert(0,'\n')
            self.add_text(['\n', len(''.join(self.pieces).strip()) * '='])
        elif level == 2:
            self.add_text(['\n', len(''.join(self.pieces).strip()) * '-'])
        elif level >= 3:
            self.pieces.insert(0, level * '#' + ' ')
        self.add_text('')
        pieces.extend(self.pieces)
        self.pieces = pieces

    def do_computeroutput(self, node):
        self.surround_parse(node, '`', '`')

    def make_constructor_list(self, constructor_nodes, classname):
        if constructor_nodes == []:
            return
        self.add_text(['\n', 'Constructors','\n','============'])
        for n in constructor_nodes:
            defn_str = classname
            argsstring = self.get_specific_subnodes(n, 'argsstring')
            if argsstring:
                defn_str = defn_str + self.get_Text(argsstring[0])
            self.add_text('\n')
            self.add_line_with_subsequent_indent(defn_str)
            pieces, self.pieces = self.pieces, ['']
            self.indent += 4
            for sn in n.childNodes:
                if sn not in self.get_specific_nodes(n, ('definition', 'name')).values():
                    self.parse(sn)
            lines = ''.join(self.pieces).split('\n')
            for i in range(len(lines)):
                if lines[i] != '':
                    lines[i] = 4 * ' ' + lines[i]
            pieces.append('\n'.join(lines))
            self.pieces = pieces
            self.indent -= 4

    def make_attribute_list(self, node):
        have_attributes = False
        for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
            if n.attributes['kind'].value == 'variable' and n.attributes['prot'].value == 'public':
                have_attributes = True
        if not have_attributes:
            return
        self.add_text(['\n', 'Attributes','\n','=========='])
        for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
            if n.attributes['kind'].value == 'variable' and n.attributes['prot'].value == 'public':
                name = self.get_Text(self.get_specific_subnodes(n, 'name')[0])
                self.add_text(['\n', name, ' : '])
                self.add_text(self.get_type(n))
                self.add_text(['\n', ''])
                pieces, self.pieces = self.pieces, ['']
                self.indent += 4
                self.parse(self.get_specific_subnodes(n, 'briefdescription')[0])
                self.parse(self.get_specific_subnodes(n, 'detaileddescription')[0])
                lines = ''.join(self.pieces).split('\n')
                for i in range(len(lines)):
                    if lines[i] != '':
                        lines[i] = 4 * ' ' + lines[i]
                pieces.append('\n'.join(lines))
                self.pieces = pieces
                self.indent -= 4

    def do_compounddef(self, node):
        kind = node.attributes['kind'].value
        if kind in ('class', 'struct'):
            prot = node.attributes['prot'].value
            if prot != 'public':
                return
            names = ('briefdescription',
                     'detaileddescription', 'includes')
            subnode = self.get_specific_subnodes(node, 'compoundname')
            self.add_text('\n\n')
            classdefn = self.get_Text(subnode[0])
            classname = classdefn.split('::')[-1]
            self.add_text('%%feature("docstring") %s "\n' % classdefn)
            constructor_nodes = []
            for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
                defn = self.get_specific_subnodes(n, 'definition')
                if defn and self.get_Text(defn[0]) == classdefn + '::' + classname:
                    constructor_nodes.append(n)
            for n in constructor_nodes:
                defn_str = classname
                argsstring = self.get_specific_subnodes(n, 'argsstring')
                if argsstring:
                    defn_str = defn_str + self.get_Text(argsstring[0])
                self.add_line_with_subsequent_indent(defn_str)

            first = self.get_specific_nodes(node, names)
            for n in ('briefdescription','detaileddescription'):
                if n in first and n != 'compoundname':
                    self.parse(first[n])

            self.make_constructor_list(constructor_nodes, classname)
            self.make_attribute_list(node)
            subnode = self.get_specific_subnodes(node, 'includes')
            if subnode:
                self.parse(subnode[0])
            self.add_text(['";', '\n'])
            for n in node.childNodes:
                names = ('compoundname', 'briefdescription','detaileddescription', 'includes')
                if n not in self.get_specific_nodes(node,names).values():
                    self.parse(n)
        elif kind in ('file', 'namespace'):
            nodes = node.getElementsByTagName('sectiondef')
            for n in nodes:
                self.parse(n)

    def do_includes(self, node):
        self.add_text('\nC++ includes: ')
        self.generic_parse(node, pad=1)

    def do_itemizedlist(self, node):
        if self.pieces != []:
            self.add_text('\n')
        listitem = self.listitem
        if self.listitem == '':
            self.listitem = '*'
        else:
            self.listitem = '-'
        self.generic_parse(node)
        self.listitem = listitem

    def do_orderedlist(self, node):
        if self.pieces != []:
            self.add_text('\n')
        listitem = self.listitem
        self.listitem = 0
        self.generic_parse(node)
        self.listitem = listitem

    def do_listitem(self, node):
        try:
            self.listitem = int(self.listitem) + 1
            item = str(self.listitem) + '. '
        except:
            item = str(self.listitem) + ' '
        self.shifted_parse(node, item)
    
    def shifted_parse(self, node, item, indent = ''):
        """Adds the following to self.pieces:
        <item><indent><parsed and wrapped line 1>
        <shift><parsed and wrapped line 2>
        ...
        <shift><parsed and wrapped line n>
        where shift is just ' ' * len(item)"""
        pieces, self.pieces = self.pieces, [indent, '']
        self.indent += len(item)
        self.generic_parse(node)
        lines = ''.join(self.pieces).split('\n')
        for i in range(len(lines)):
            if lines[i] != '':
                lines[i] = len(item) * ' ' + lines[i]
        pieces.append(item + '\n'.join(lines)[len(item):])
        self.pieces = pieces
        self.indent -= len(item)

    def add_line_with_subsequent_indent(self, value, indent=4):
        if isinstance(value, (list, tuple)):
            line = ''.join(value)
        else:
            line = value
        line = line.strip()
        width = self.textwidth-self.indent-indent
        wrapped_lines = textwrap.wrap(line[indent:], width=width)
        for i in range(len(wrapped_lines)):
            if wrapped_lines[i] != '':
                wrapped_lines[i] = indent * ' ' + wrapped_lines[i]
        self.pieces.append(line[:indent] + '\n'.join(wrapped_lines)[indent:] + '  \n')

    def textwrap_para(self, dont_end_with_newline = True):
        """Wrap self.pieces to a width of self.textwidth - self.indent, but preserve
        newlines and "markdown newlines"."""
        width = self.textwidth - self.indent
        pieces, self.pieces = self.pieces, []
        for line in ''.join(pieces).splitlines():
            keep_markdown_newline = line[-2:] == '  '
            wrapped_line = textwrap.wrap(line, width=width, break_long_words=False)
            if wrapped_line == []:
                wrapped_line = ['']
            if keep_markdown_newline:
                wrapped_line[-1] = wrapped_line[-1] + '  '
            for wl in wrapped_line:
                self.pieces.extend([wl, '\n'])
        if dont_end_with_newline:
            self.pieces = self.pieces[:-1]
        else:
            self.pieces = self.pieces[:-1] + ['  \n']

    def do_para(self, node):
        if self.pieces[-1:] == ['']:
            pieces = []
        else:
            if self.pieces != []:
                self.add_text('\n')
            pieces, self.pieces = self.pieces, []
        self.generic_parse(node)
        self.textwrap_para(self.pieces[-1:] == [''])
        pieces.extend(self.pieces)
        self.pieces = pieces

    def do_parameterlist(self, node):
        if self.pieces != []:
            if self.pieces[-1][-1:] == '\n':
                self.add_text('\n')
            else:
                self.add_text('  \n\n')
        text = 'unknown'
        for key, val in node.attributes.items():
            if key == 'kind':
                if val == 'param':
                    text = 'Parameters'
                elif val == 'exception':
                    text = 'Exceptions'
                elif val == 'retval':
                    text = 'Returns'
                else:
                    text = val
                break
        self.add_text([text, '\n', len(text) * '-', '\n'])
        self.generic_parse(node)

    def do_parameteritem(self, node):
        pieces, self.pieces = self.pieces, []
        self.generic_parse(node)
        pieces.extend(self.pieces)
        self.pieces = pieces

    def do_parameternamelist(self, node):
        self.generic_parse(node)
        self.add_text(' : ')
    
    def do_parametername(self, node):
        if self.pieces != []:
            self.add_text(', ')
        try:
            data = self.get_Text(node)
        except AttributeError:  # perhaps a <ref> tag in it
            data = self.get_Text(node.firstChild)
        self.add_text(data)

    def do_parameterdescription(self, node):
        names = ''.join(self.pieces)
        self.pieces = []
        self.shifted_parse(node, '% 4s' % names[:4], names[4:])

    def do_parameterdefinition(self, node):
        self.generic_parse(node, pad=1)

    def do_detaileddescription(self, node):
        self.generic_parse(node)

    def do_briefdescription(self, node):
        self.generic_parse(node)

    def do_memberdef(self, node):
        prot = node.attributes['prot'].value
        id = node.attributes['id'].value
        kind = node.attributes['kind'].value
        tmp = node.parentNode.parentNode.parentNode
        compdef = tmp.getElementsByTagName('compounddef')[0]
        cdef_kind = compdef.attributes['kind'].value

        if prot != 'public':
            return
        first = self.get_specific_nodes(node, ('definition', 'name'))
        argsstring = self.get_specific_subnodes(node, 'argsstring')
        if argsstring:
            argsstring = self.get_Text(argsstring[0])
        else:
            argsstring = ''
        name = self.get_Text(first['name'])
        type = self.get_type(node)
        if name[:8] == 'operator':  # Don't handle operators yet.
            return
        if not 'definition' in first or kind in ['variable', 'typedef']:
            return

        function_definition = name + argsstring
        if type != '' and type != 'void':
            function_definition = function_definition + ' -> ' + type
        self.add_text('\n')
        self.add_text('%feature("docstring") ')

        anc = node.parentNode.parentNode
        if cdef_kind in ('file', 'namespace'):
            ns_node = anc.getElementsByTagName('innernamespace')
            if not ns_node and cdef_kind == 'namespace':
                ns_node = anc.getElementsByTagName('compoundname')
            if ns_node:
                ns = self.get_Text(ns_node[0])
                self.add_text(' %s::%s "\n' % (ns, name))
            else:
                self.add_text(' %s "\n' % (name))
        elif cdef_kind in ('class', 'struct'):
            # Get the full function name.
            anc_node = anc.getElementsByTagName('compoundname')
            cname = self.get_Text(anc_node[0])
            self.add_text(' %s::%s "\n' % (cname, name))

        if self.include_function_definition:
            self.add_line_with_subsequent_indent(function_definition)
            self.add_text('\n')
        self.add_text('')
        for n in node.childNodes:
            if n not in first.values():
                self.parse(n)
        self.add_text(['";', '\n'])

    def do_definition(self, node):
        data = self.get_Text(node)
        self.add_text('%s "\n%s' % (data, data))

    def do_sectiondef(self, node):
        kind = node.attributes['kind'].value
        if kind in ('public-func', 'func', 'user-defined', ''):
            self.generic_parse(node)

    def do_header(self, node):
        """For a user defined section def a header field is present
        which should not be printed as such, so we comment it in the
        output."""
        data = self.get_Text(node)
        self.add_text('\n/*\n %s \n*/\n' % data)
        # If our immediate sibling is a 'description' node then we
        # should comment that out also and remove it from the parent
        # node's children.
        parent = node.parentNode
        idx = parent.childNodes.index(node)
        if len(parent.childNodes) >= idx + 2:
            nd = parent.childNodes[idx + 2]
            if nd.nodeName == 'description':
                nd = parent.removeChild(nd)
                self.add_text('\n/*')
                self.generic_parse(nd)
                self.add_text('\n*/\n')

    def do_simplesect(self, node):
        kind = node.attributes['kind'].value
        if kind in ('date', 'rcs', 'version'):
            return
        if self.pieces != []:
            self.add_text('\n')
        if kind == 'warning':
            self.shifted_parse(node, 'WARNING: ')
        elif kind == 'see':
            self.shifted_parse(node, 'See: ')
        elif kind == 'return':
            self.add_text(['Returns:', '\n', '--------', '\n'])
            self.shifted_parse(node, '    ')
        else:
            self.generic_parse(node)

    def do_argsstring(self, node):
        self.generic_parse(node, pad=1)

    def do_member(self, node):
        kind = node.attributes['kind'].value
        refid = node.attributes['refid'].value
        if kind == 'function' and refid[:9] == 'namespace':
            self.generic_parse(node)

    def do_doxygenindex(self, node):
        self.multi = 1
        comps = node.getElementsByTagName('compound')
        for c in comps:
            refid = c.attributes['refid'].value
            fname = refid + '.xml'
            if not os.path.exists(fname):
                fname = os.path.join(self.my_dir,  fname)
            if not self.quiet:
                print("parsing file: %s" % fname)
            p = Doxy2SWIG(fname, self.include_function_definition, self.quiet)
            p.generate()
            self.pieces.extend(p.pieces)

    def write(self, fname):
        o = my_open_write(fname)
        o.write(''.join(self.pieces))
        o.write('\n')
        o.close()

def convert(input, output, include_function_definition=True, quiet=False):
    p = Doxy2SWIG(input, include_function_definition, quiet)
    p.generate()
    p.write(output)


def main():
    usage = __doc__
    parser = optparse.OptionParser(usage)
    parser.add_option("-n", '--no-function-definition',
                      action='store_true',
                      default=False,
                      dest='func_def',
                      help='do not include doxygen function definitions')
    parser.add_option("-q", '--quiet',
                      action='store_true',
                      default=False,
                      dest='quiet',
                      help='be quiet and minimize output')

    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("error: no input and output specified")

    convert(args[0], args[1], not options.func_def, options.quiet)


if __name__ == '__main__':
    main()
