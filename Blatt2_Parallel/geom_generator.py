import sys
import csv
import math
from math import sqrt

length = 0
height = 0
res_x = 0 
res_y = 0
rows = []    

def bresenham(x0, y0, x1, y1):
  dx =  abs(x1-x0)
  sx = 1 if ( x0<x1 ) else -1
  dy = -abs(y1-y0)
  sy = 1 if y0<y1 else -1
  err = dx+dy, e2; # error value e_xy

  while(1):
    setPixel(x0,y0);
    
    if (x0==x1 and y0==y1):
        break
    
    e2 = 2*err;
    if (e2 > dy): 
        err += dy
        x0 += sx # e_xy+e_x > 0 
    if (e2 < dx):
        err += dx
        y0 += sy # e_xy+e_y < 0
        
def draw_simple_karman():
    x = int(math.floor(res_y/2))
    y = int(math.floor(res_y/2))
    br = 4
    #rows[x][y] = 1
    
    for offset in range (0, int(math.floor((x/2)/sqrt(2))) + 1):
        for off_x in range (-1, 3):
            rows[y+offset][x-offset + off_x] = 1
            rows[y-offset][x+offset + off_x] = 1
            
    offset = int(math.floor((x/2)/sqrt(2)))
    for off_x in range (-1, 2):
        rows[y+offset+1][x-offset + off_x] = 1
        rows[y-offset-1][x+offset + 1 + off_x] = 1
    
    return rows


#read params
if __name__ == "__main__":
    filename = sys.argv[1]
    inputfilename = sys.argv[2]

    type = int(sys.argv[3])


    fobj = open(inputfilename)
    i = 1
    for line in fobj:
	if (i==4):
            res_x = int(line.rstrip()) +2
	if (i==5):
            res_y = int(line.rstrip()) +2
 	if (i==6):
            length = int(line.rstrip())
	if (i==7):
            height = int(line.rstrip())
        i = i+1
    fobj.close()
    
    alpha = 0
    
    if ( type == 2 ):
        alpha = sys.argv[4]
        
    if ( type > 2 or type < 0):
        print "Error: invalid type specified"
        exit()
    
    #generate normal tunnel
    for y in range(0, res_y):
        row = []
        for col in range(0, res_x):
            val = 0;
            if ( y == 0 or y == res_y - 1 ):
                val = 1
            elif col == 0:
                if type == 2:
                    val = 5
                else:
                    val = 3
            elif col == res_x - 1:
                val = 4
            
            row.append(val)        
        rows.append(row)
    
    #generate obstacle
    if ( type == 1):
        for y in range(0, res_y):
            if y >= res_y/2:
                for col in range(0, res_x):
                    if col < res_y/2:
                        rows[y][col] = 1
                        
    
    #generate karman
    if ( type == 2):
        rows = draw_simple_karman()
                
        
    
    f = open(filename, 'wb')
    writer = csv.writer(f)
    writer.writerows(rows)
    f.close()
    
