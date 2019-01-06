#!/usr/bin/env python

from sys import path
path.append("../../src/modules/python/lncrnas/")
import dataunite as du

args = du.parsing_args()
table = du.factorise(args)
#counter = 0
#for k,v in table.iteritems():
#    print k + " : " + "\t".join(v)
#    counter += 1
#    #if counter > 150: break

if du.write_out(table, "output.txt"):
    print "[SUCCESS!] - Data was stored at output.txt file"
