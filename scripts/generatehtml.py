#!/usr/bin/env python
import os,sys
import getopt
import commands


#print usage
def usage() :
    print ' '
    print 'generateJSONplotterFromList.py [options]'
    print '  -o : output file (plotter.json by default'
    print ' '

#parse the options
try:
     # retrive command line options
     shortopts  = "o:h?"
     opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
     # print help information and exit:
     print "ERROR: unknown options in argument %s" % sys.argv[1:]
     usage()
     sys.exit(1)


output=''
for o,a in opts:
    if o in("-?", "-h"):
        usage()
        sys.exit(1)
    elif o in('-o'): output=a

if(len(output)==0) :
    usage()
    sys.exit(1)


status, pwd = commands.getstatusoutput('which generatehtml.py')
pwd = pwd.split('RootCoreBin')[0]+'DMHiggsAnalysis/scripts/'

status, user=commands.getstatusoutput('whoami')

fromfile = pwd+'index_template.html'
tofile = open(output+'/index.html',"w")
status, date = commands.getstatusoutput('date')

status, firstname=commands.getstatusoutput('phonebook `whoami` -t firstname')
firstname=firstname.split(';')[0]
status, surname=commands.getstatusoutput('phonebook `whoami` -t surname')
surname=surname.split(';')[0]
fullname = firstname+' '+surname

with open(fromfile) as fp:

    for line in fp:
        if 'Created on DateToInsert by' in line:
                line = line.replace('DateToInsert',date)
		line = line.replace('NameToInsert',fullname)
	tofile.writelines(line)

tofile.close()

