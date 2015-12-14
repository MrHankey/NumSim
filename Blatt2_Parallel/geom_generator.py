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


#generate parameter file
def generateParamFile(fileInG, fileInP, fileOut):
    # Initialize
    dataG = []
    dataP = []
    fInG = open(fileInG, 'r')
    for line in fInG:
        dataG.append(line.rstrip())
    fInG.close()
    
    fInP = open(fileInP, 'r')
    for line in fInP:
        dataP.append(line.rstrip())
    fInP.close()
    
    
    writeString  = ""
    writeString += "Geometrie:      " + "iMax=" + dataG[3] + ", jMax="  + dataG[4] + ", xLength="  + dataG[5] + ", yLength=" + dataG[6] + "\n" ##
    writeString += "Zeitsteuerung:  " + "tEnd=" + dataP[4] + ", tau="   + dataP[6] + ", deltaT="   + dataP[3]                           + "\n" ##
    writeString += "Solver:         " + "eps="  + dataP[5] + ", omega=" + dataP[1] + ", alpha="    + dataP[2] + ", iterMax=" + dataP[7] + "\n" ##
    writeString += "Kraefte und RE: " + "GX="   + str(0)   + ", GY="    + str(0)   + ", RE="       + '10000' + ', REstufe=100' +        + "\n" ##
    writeString += "Anfangswerte:   " + "UI="   + dataG[0] + ", VI="    + dataG[1] + ", PI="       + dataG[2]                           + "\n" ##
    
    
    fOut = open(fileOut,'w') 
    fOut.write(writeString)
    fOut.close()


#read params
if __name__ == "__main__":
    filename = sys.argv[1]
    in_geom = sys.argv[2]
    in_param = sys.argv[3]

    type = int(sys.argv[4])

	#generate parameter file
    #generateParamFile(inputfilename,"parameter.txt","param.swag")
    
    f = open(in_param, 'r+')
    data = []
    for line in f:
        data.append(line.rstrip())
        
     
    #text = f.read()
    f.seek(0)
    if type == 2:
        data[0] = '10000'
    elif type == 3 :
        data[0] = '1000'
    else:
        data[0] = '100'
        
    writeString = ''
    for line in data:
        writeString += line + '\n'
        
    f.write(writeString)
    
    f.truncate()
    f.close()
	

    fobj = open(in_geom)
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
            
    if ( type > 3 or type < 0):
        print "Error: invalid type specified"
        exit()
    
    #generate driven cavity
    if ( type == 3):
        for y in range(0, res_y):
            row = []
            for col in range(0, res_x):
                val = 0;
                if ( y == 0 ):
                    val = 5
                elif col == 0 or col == res_x - 1 or y == res_y - 1:
                    val = 1
                
                row.append(val)
            rows.append(row)
    else:
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
    
