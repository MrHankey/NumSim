import sys
import csv

#read params
filename = sys.argv[1]
length = sys.argv[2]
height = sys.argv[3]
res_x = sys.argv[4]
res_y = sys.argv[5]
type = sys.argv[6]

alpha = 0

if ( type == 2 ):
    alpha = sys.argv[7]
    
if ( type > 2 or type < 0):
    print "Error: invalid type specified"
    exit()
    

    
