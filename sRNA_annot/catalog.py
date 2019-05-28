#-*-coding:utf-8-*-
import os
import sys
from optparse import OptionParser


def main():
    """
    %prog [options]
    """
    parser = OptionParser()
    parser.add_option('-d', '--database', help="database for the project")
    parser.add_option('-p','--prefix',help="the prefix")
    parser.add_option('-m','--mismatch',help="mismatch")
    parser.add_option('-n', '--novel_threshold', help="novel count threshold")
    opts, args = parser.parse_args()
    if (opts.database == None or opts.prefix == None or opts.mismatch == None):
        print('\033[0;31;40m%s\033[0m' % "Warning: database, prefix, mismatch and novel_threshold must be given\n")
        sys.exit(parser.print_help())

def catalog(database,prefix,mimatch,novel):
	OA = open(prefix+".catalog.xls",'w')
	OA.writelines("sRNA id\tCount\tType\tDescription\n")


