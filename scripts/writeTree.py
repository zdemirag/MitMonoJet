from argparse import ArgumentParser
import re
import os

def branchType(code):
    if code == 'C':
        return 'Text_t const*'
    elif code == 'B':
        return 'Char_t'
    elif code == 'b':
        return 'UChar_t'
    elif code == 'S':
        return 'Short_t'
    elif code == 's':
        return 'UShort_t'
    elif code == 'I':
        return 'Int_t'
    elif code == 'i':
        return 'UInt_t'
    elif code == 'L':
        return 'Long64_t'
    elif code == 'l':
        return 'ULong64_t'
    elif code == 'F':
        return 'Float_t'
    elif code == 'D':
        return 'Double_t'
    elif code == 'O':
        return 'Bool_t'


PACKAGE = 'MitMonoJet/DataFormats'
FULLPATH = os.environ['CMSSW_BASE'] + '/src/MitMonoJet/DataFormats'

argParser = ArgumentParser(description = 'Generate C++ code for a flat tree')
argParser.add_argument('config', metavar = 'CONFIG')

args = argParser.parse_args()

objPat = re.compile('^\\[([^:\\]]+)(?:|\\:(MAX=[0-9]+))\\]')
brnPat = re.compile('^([a-zA-Z_][a-zA-Z0-9_]*)/([CBbSsIiLlFDO])')
fncPat = re.compile('^[a-zA-Z_][^ ]* +[a-zA-Z_][a-zA-Z0-9_]* *\\([^\\)]*\\) *\{ *return +[^;]+; *\}')
incPat = re.compile('#include [^ ]+')

includes = []
objs = [] # keep object declarations in order of appearance in the input
defs = {}
sizes = {}
currentObj = ''

namespace = os.path.basename(args.config)
namespace = namespace[:namespace.find('.')]

with open(args.config) as branchList:
    for line in branchList:
        line = line.strip()
        if objPat.match(line):
            matches = objPat.match(line)
            currentObj = matches.group(1)
            objs.append(currentObj)
            defs[currentObj] = ([], []) # branches and functions
            if matches.group(2):
                sizes[currentObj] = int(matches.group(2).replace('MAX=', ''))
            else:
                sizes[currentObj] = 32

        elif brnPat.match(line):
            if not currentObj:
                raise RuntimeError('Branch given before object name')

            matches = brnPat.match(line)
            defs[currentObj][0].append((matches.group(1), matches.group(2)))

        elif fncPat.match(line):
            if not currentObj:
                raise RuntimeError('Branch given before object name')

            defs[currentObj][1].append(line)

        elif incPat.match(line):
            includes.append(line)
        else:
            if line and not line.startswith('%'):
                print 'Skipping unrecognized pattern:', line


if not os.path.isdir(FULLPATH + '/interface'):
    os.makedirs(FULLPATH + '/interface')
if not os.path.isdir(FULLPATH + '/src'):
    os.makedirs(FULLPATH + '/src')
if not os.path.exists(FULLPATH + '/BuildFile.xml'):
    with open(FULLPATH + '/BuildFile.xml', 'w') as buildFile:
        buildFile.write('<use name="root"/>\n')
        buildFile.write('<export>\n')
        buildFile.write('  <lib name="1"/>\n')
        buildFile.write('</export>\n')

try:
    eventDef = defs.pop('Event')
except KeyError:
    eventDef = ([('runNumber', 'i'), ('lumiNumber', 'i'), ('eventNumber', 'i')], [])

arrNames = []
for obj in objs:
    if obj.lower().endswith('s') or obj.lower().endswith('x'):
        mod = 'es'
    else:
        mod = 's'
    arrNames.append(obj.lower() + mod)

# Object headers
with open(FULLPATH + '/interface/Objects.h', 'w') as header:
    header.write('#ifndef ' + namespace + '_Objects_h\n')
    header.write('#define ' + namespace + '_Objects_h\n')
    for inc in includes:
        header.write(inc + '\n')

    header.write('#include "Rtypes.h"\n\n')

    header.write('namespace ' + namespace + ' {\n\n')

    for obj in objs:
        header.write('  class ' + obj + 'Collection;\n')

    header.write('\n')

    for obj in objs:
        header.write('  class ' + obj + ' {\n')
        header.write('  public:\n')
        header.write('    ' + obj + '(' + obj + 'Collection&, UInt_t idx);\n')

        header.write('\n')
        
        for func in defs[obj][1]:
            header.write('    ' + func + '\n')

        if len(defs[obj][1]):
            header.write('\n')

        for brName, brType in defs[obj][0]:
            header.write('    ' + branchType(brType) + '& ' + brName + ';\n')

        header.write('  };\n\n')

    header.write('}\n\n')

    header.write('#endif\n')

# Collection haders
with open(FULLPATH + '/interface/Collections.h', 'w') as header:
    header.write('#ifndef ' + namespace + '_Collections_h\n')
    header.write('#define ' + namespace + '_Collections_h\n')
    header.write('#include "' + PACKAGE + '/interface/Objects.h"\n\n')

    header.write('class TTree;\n\n')

    header.write('namespace ' + namespace + ' {\n\n')

    for obj in objs:
        header.write('  class ' + obj + 'Collection {\n')
        header.write('  public:\n')
        header.write('    static UInt_t const NMAX = ' + str(sizes[obj]) + ';\n')
        header.write('    typedef ' + namespace + '::' + obj + ' value_type;\n')
        header.write('    typedef value_type& reference;\n')
        header.write('    typedef value_type const& const_reference;\n')
        header.write('    typedef ' + namespace + '::' + obj + '* iterator;\n')
        header.write('    typedef ' + namespace + '::' + obj + ' const* const_iterator;\n\n')

        header.write('    ' + obj + 'Collection();\n')
        header.write('    ~' + obj + 'Collection();\n\n')

        header.write('    reference at(UInt_t idx);\n')
        header.write('    const_reference at(UInt_t idx) const;\n')
        header.write('    void clear() { resize(0); }\n')
        header.write('    void resize(UInt_t size);\n')
        header.write('    iterator begin() { return array_; }\n')
        header.write('    const_iterator begin() const { return array_; }\n')
        header.write('    iterator end() { return array_ + size; }\n')
        header.write('    const_iterator end() const { return array_ + size; }\n\n')

        header.write('    void setAddress(TTree&);\n')
        header.write('    void book(TTree&);\n\n')

        header.write('    UInt_t size = 0;\n')
        for brName, brType in defs[obj][0]:
            header.write('    ' + branchType(brType) + ' ' + brName + '[NMAX];\n')

        header.write('\n')

        header.write('  private:\n')
        header.write('    value_type* array_;\n\n')

        header.write('  };\n\n')

    header.write('}\n\n')
    header.write('#endif\n')

# Event header
with open(FULLPATH + '/interface/Event.h', 'w') as header:
    header.write('#ifndef ' + namespace + '_Event_h\n')
    header.write('#define ' + namespace + '_Event_h\n')
    header.write('#include "' + PACKAGE + '/interface/Collections.h"\n\n')
    header.write('namespace ' + namespace + ' {\n\n')
    header.write('  class Event {\n')
    header.write('  public:\n')

    for func in eventDef[1]:
        header.write('    ' + func + '\n')
    
    if len(eventDef[1]):
        header.write('\n')

    for brName, brType in eventDef[0]:
        header.write('    ' + branchType(brType) + ' ' + brName + ';\n')

    header.write('\n')
    
    for iO, obj in enumerate(objs):
        header.write('    ' + obj + 'Collection ' + arrNames[iO] + ';\n')

    header.write('\n')
    header.write('    void setAddress(TTree&);\n')
    header.write('    void book(TTree&);\n')
    header.write('  };\n\n')
    header.write('}\n\n')
    header.write('#endif\n')

# Object source
with open(FULLPATH + '/src/Objects.cc', 'w') as src:
    src.write('#include "' + PACKAGE + '/interface/Objects.h"\n\n')
    src.write('#include "' + PACKAGE + '/interface/Collections.h"\n\n')

    for obj in objs:
        src.write(namespace + '::' + obj + '::' + obj + '(' + obj + 'Collection& col, UInt_t idx) :\n')
        for iB, (brName, brType) in enumerate(defs[obj][0]):
            src.write('  ' + brName + '(col.' + brName + '[idx])')
            if iB == len(defs[obj][0]) - 1:
                src.write('\n')
            else:
                src.write(',\n')

        src.write('{\n}\n\n')

# Collection source
with open(FULLPATH + '/src/Collections.cc', 'w') as src:
    src.write('#include "' + PACKAGE + '/interface/Collections.h"\n')
    src.write('#include "TTree.h"\n')
    src.write('#include <stdexcept>\n')
    src.write('#include <memory>\n\n')

    for obj in objs:
        src.write(namespace + '::' + obj + 'Collection::' + obj + 'Collection() :\n')
        src.write('  array_(std::allocator<' + obj + '>().allocate(NMAX))\n')
        src.write('{\n')
        src.write('  for (unsigned iP(0); iP != NMAX; ++iP)\n')
        src.write('    new (array_ + iP) ' + obj + '(*this, iP);\n')
        src.write('}\n\n')

        src.write(namespace + '::' + obj + 'Collection::~' + obj + 'Collection()\n')
        src.write('{\n')
        src.write('  std::allocator<' + obj + '>().deallocate(array_, NMAX);\n')
        src.write('}\n\n')

        src.write(namespace + '::' + obj + 'Collection::reference\n')
        src.write(namespace + '::' + obj + 'Collection::at(UInt_t _idx)\n')
        src.write('{\n')
        src.write('  if (_idx < size)\n')
        src.write('    return array_[_idx];\n')
        src.write('  throw std::out_of_range("' + obj + 'Collection::at");\n')
        src.write('}\n\n')

        src.write(namespace + '::' + obj + 'Collection::const_reference\n')
        src.write(namespace + '::' + obj + 'Collection::at(UInt_t _idx) const\n')
        src.write('{\n')
        src.write('  if (_idx < size)\n')
        src.write('    return array_[_idx];\n')
        src.write('  throw std::out_of_range("' + obj + 'Collection::at");\n')
        src.write('}\n\n')

        src.write('void\n')
        src.write(namespace + '::' + obj + 'Collection::resize(UInt_t _size)\n')
        src.write('{\n')
        src.write('  if (size < NMAX) {\n')
        src.write('    size = _size;\n')
        src.write('    return;\n')
        src.write('  }\n')
        src.write('  throw std::length_error("' + obj + 'Collection::resize");\n')
        src.write('}\n\n')

        prefix = obj.lower()

        src.write('void\n')
        src.write(namespace + '::' + obj + 'Collection::setAddress(TTree& _tree)\n')
        src.write('{\n')
        src.write('  _tree.SetBranchAddress("' + prefix + '.size", &size);\n')
        for brName, brType in defs[obj][0]:
            src.write('  _tree.SetBranchAddress("' + prefix + '.' + brName + '", ' + brName + ');\n')
        src.write('}\n\n')

        src.write('void\n')
        src.write(namespace + '::' + obj + 'Collection::book(TTree& _tree)\n')
        src.write('{\n')
        src.write('  _tree.Branch("' + prefix + '.size", &size, "size/i");\n')
        for brName, brType in defs[obj][0]:
            src.write('  _tree.Branch("' + prefix + '.' + brName + '", ' + brName + ', "' + brName + '[' + prefix + '.size]/' + brType + '");\n')
        src.write('}\n\n')

# Event source
with open(FULLPATH + '/src/Event.cc', 'w') as src:
    src.write('#include "' + PACKAGE + '/interface/Event.h"\n')
    src.write('#include "TTree.h"\n\n')
    src.write('void\n')
    src.write(namespace + '::Event::setAddress(TTree& _tree)\n')
    src.write('{\n')
    for brName, brType in eventDef[0]:
        src.write('  _tree.SetBranchAddress("' + brName + '", &' + brName + ');\n')

    src.write('\n')
    
    for arrName in arrNames:
        src.write('  ' + arrName + '.setAddress(_tree);\n')

    src.write('}\n\n')

    src.write('void\n')
    src.write(namespace + '::Event::book(TTree& _tree)\n')
    src.write('{\n')
    for brName, brType in eventDef[0]:
        src.write('  _tree.Branch("' + brName + '", &' + brName + ', "' + brName + '/' + brType + '");\n')

    src.write('\n')
    
    for arrName in arrNames:
        src.write('  ' + arrName + '.book(_tree);\n')

    src.write('}\n\n')
