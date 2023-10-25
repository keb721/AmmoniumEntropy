##############################
# Code to create ordered ice #
##############################

import numpy as np

pos = np.genfromtxt("ice_XI_211.xyz", skip_header=4, skip_footer=1)

newname = "ordered_o211.xyz"

x = pos[-3, 1] ; y = pos[-2, 2] ; z = pos[-1, 3]
ix = 1/x ; iy = 1/y ; iz = 1/z

uv = [x, y, z]
iv = [ix, iy, iz]

pos  = pos[:-3, 1:]
Opos = pos[::3, :]

newpos = np.zeros((len(pos), 3))
newpos[::3, :] = pos[::3, :]

conns = np.zeros((len(Opos), 4))
conns[:, :] = -1

# Now we want to make sure that hydrogens are either primarily in the +x, +y or -z directions

def min_image(vec, uv, iv):
    for i in range(3):
        vec[i] = vec[i] - uv[i]*round(vec[i]*iv[i])
    return vec


for i in range(len(pos)):
    if i%3 == 0 :
        continue
    mindist = 40
    minO = 40
    for j in range(len(Opos)):
        if j == int(np.floor(i/3)):
            O = int(np.floor(i/3))
            continue
        diff = pos[i, :] - pos[3*j, :]
        diff = min_image(diff, uv, iv)
        dist = np.sqrt(np.dot(diff, diff))
        if dist < mindist:  
            mindist = dist
            minO = j
    if i%3==1:
        conns[O,0] = minO
    else:
        conns[O,1] = minO
    
    if conns[minO, 2] == -1:
        conns[minO, 2] = O
    else:
        conns[minO, 3] = O

for i in range(len(Opos)):
    #print("Oxygen", i)
    for j in conns[i, :]:
        if j == - 1:
            continue
        tmp = Opos[int(j), :] - Opos[i, :]
        tmp = min_image(tmp, uv, iv)
        if tmp[2] < -2:
            #print("Adding", tmp)
            newpos[3*i+1] = newpos[3*i] +1/3*tmp
            #print(np.where(conns[int(j), :] == i))
            conns[int(j), np.where(conns[int(j), :] == i)[0][0]] = -1
            # Now add the one that's +x, and the others will fall into place
            for k in conns[i, :]:
                tmp = Opos[int(k), :] - Opos[i, :]
                tmp = min_image(tmp, uv, iv)
                if tmp[0] > 2:
                    #print("Adding", tmp)
                    newpos[3*i+2] =  newpos[3*i] +1/3*tmp
                    #print(k, conns[int(k), :], np.where(conns[int(k), :] == i))
                    conns[int(k), np.where(conns[int(k), :] == i)[0][0]] = -1
                    conns[i, :] = -1
                    break
                    

#print(newpos)
#print(conns) 


for i in range(len(Opos)):
    #print("Oxygen", i)
    curr = 1
    for j in conns[i, :]:
        if j == - 1:
            continue
        tmp = Opos[int(j), :] - Opos[i, :]
        tmp = min_image(tmp, uv, iv)
        #print("Adding", tmp)
        newpos[3*i+curr] = newpos[3*i] +0.9752*tmp/np.sqrt(np.dot(tmp, tmp))
        curr += 1


#print(newpos)


f = open(newname, 'w')
f.write("# \n# \n"+str(len(newpos))+"\n%PBC\n")
for i in range(len(newpos)):
    if i%3==0:
        f.write("O \t" +str(newpos[i, 0])+"\t"+str(newpos[i, 1])+"\t"+str(newpos[i, 2])+"\n")
    else:
        f.write("H \t" +str(newpos[i, 0])+"\t"+str(newpos[i, 1])+"\t"+str(newpos[i, 2])+"\n")

f.write("\n x \t "+str(x)+"0 \t 0 \t")
f.write("\n y \t 0 \t "+str(y)+"0")
f.write("\n z \t 0 \t 0 \t"+str(z))
f.write("\n offset \t 0 \t 0 \t"+str(-1.05))   # Allow offset for NH4

f.close()
