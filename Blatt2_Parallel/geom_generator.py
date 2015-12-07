import sys
import csv    

#read params
if __name__ == "__main__":
    filename = sys.argv[1]
    length = float(sys.argv[2])
    height = float(sys.argv[3])
    res_x = int(sys.argv[4])
    res_y = int(sys.argv[5])
    type = int(sys.argv[6])
    
    alpha = 0
    
    if ( type == 2 ):
        alpha = sys.argv[7]
        
    if ( type > 2 or type < 0):
        print "Error: invalid type specified"
        exit()
    
    rows = []
    #generate normal tunnel
    for y in range(0, res_y):
        row = []
        for col in range(0, res_x):
            val = 0;
            if ( y == 0 or y == res_y - 1 ):
                val = 1
            elif col == 0:
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
        for y in range(0, res_y):
            if y <= res_y/2:
                for x in range(0, res_x):
                    if col < res_y/2:
                        rows[y][col] = 1
        
    
    f = open(filename, 'wb')
    writer = csv.writer(f)
    writer.writerows(rows)
    f.close()
    
