#!/usr/bin/env python3

import numpy as np
import os 
import math

#### Read Inputs ####

## Input paths into a list
localpath = os.getcwd()
inputpath = os.path.join(localpath,"Inputs")
os.chdir(inputpath)


## Determine number of components to be rotated
parents = open('parent.txt','r')  
parentnames = parents.read().split()
num_of_components = len(parentnames)
parents.close()

## Assign values of all inputs 
deflections =open('deflections.txt')
deflec_val = deflections.read().split()
deflec = np.array(deflec_val)
deflec = deflec.astype(np.float)

hp1_coordinates = open('hinge_points1.txt','r')
hp1_val = hp1_coordinates.read().split()
hp1 = np.array(hp1_val)
hp1 = hp1.reshape(3,num_of_components) #Shape to matrix
hp1 = hp1.astype(np.float)

hp2_coordinates = open('hinge_points2.txt','r')
hp2_val = hp2_coordinates.read().split()
hp2 = np.array(hp2_val)
hp2 = hp2.reshape(3,num_of_components) # Shanpe to matrix 
hp2 = hp2.astype(np.float)

## Check to make sure there are the correct number of values in each input 
if len(deflec_val) != len(parentnames):
    print("Numeber of components and deflections do not match")  
    exit() 
      
if len(hp1_val) != len(parentnames)*3:
    print("Numeber of components and first hinge points do not match")  
    exit() 
     
if len(hp2_val) != len(parentnames)*3:
    print("Numeber of components and second hinge points do not match")  
    exit() 
         
os.chdir(localpath)
#### End Inputs ####

world = np.array([0,0,0]) #define world orgin
h = np.subtract(hp2,hp1)



for plcv in range(len(parentnames)):
    
    
    
    hi = h[plcv,:]  # Define hinge point 1 for stl 
    p2 = hp2[plcv,:]  # Define hinge point 2 for stl 
    filename = parentnames[plcv] #STL being opened
    
    ### Read in values/attributes of mesh file ###
    
    meshdir = os.path.join(localpath, "Components(stl files)")
    meshpath = os.path.join(meshdir,filename)
    mesh_stl = open(meshpath,'r')    #Open stl file
    data=mesh_stl.readlines()        #Read all lines of stl file
    lst=[]                      #Create empty list so self-updating is possible
    nst=[]
    nfacet=0
    for i in range(len(data)):
        xx='vertex' in data[i] #See if the word 'vertex' exists within each line
        yy= 'normal' in data[i] #See if the word 'normal' exists within each line
        if xx==1:  
            lst.append(data[i]) #If so, add to the list
            vertices=np.zeros([len(lst),3])    #Create matrix for vertices
        elif yy==1:
             nst.append(data[i])
             normalStore = np.zeros([len(nst),3])
            

    for i in range(len(lst)):
        st=lst[i]
        sta=st.strip()          #Remove all trailing spaces in line
        sta=sta.replace('vertex','').strip()    #Remove the word 'vertex'
        sta=np.array(sta.split())         #Split remaining string into 1x3 vector containing vertices
        vertices[i,:]=sta                 #Add vertices to matrix

    for i in range(len(nst)):
        nt=nst[i]
        nta=nt.strip()          #Remove all trailing spaces in line
        nta=nta.replace('facet normal','').strip()    #Remove the word 'facet normal'
        nta=np.array(nta.split())         #Split remaining string into 1x3 vector containing vertices
        normalStore[i,:]=nta                 #Add vertices to matrix
        
   
    for i in range(len(data)):
        j='endfacet' in data[i]
        if j==1:
            nfacet=nfacet+1         #Find the number of facets

    v1 = vertices[0::3] #Store first vertex for normal calculations
    
    ### Deflection ###
     
    magh = math.sqrt(hi[0]**2 + hi[1]**2 + hi[2]**2)
    
    ## Direction cosines to world body directions
    alph = (180/np.pi)*math.acos(hi[0]/magh)
    bet = (180/np.pi)*math.acos(hi[1]/magh)
    gam = (180/np.pi)*math.acos(hi[2]/magh)
    
    ## Account for any singularities in DCM's 
    if alph == 90:
        alph = 90-.0000001
    else: 
        pass 

    if gam == 180:
        gam = 180 - 0.000001
    elif gam == 0:
        gam = 0 + 0.000001
    else:
        pass
    
    delta = [p2 - world] #Translation vector from world orgin to a point on hinge line
    delta =np.array(delta)
    delta = np.transpose(delta)
    #First rotation angle 
    theta1 = (180/np.pi)*math.acos(((math.cos(alph*np.pi/180)**2) + (math.sin(gam*np.pi/180)**2) - (math.cos(bet*np.pi/180)**2))/(2*math.cos(alph*np.pi/180)*math.sin(gam*np.pi/180)))   
    #Second rotation angle
    theta2 = 90 - gam
    
    #===== Determine octant to which the hinge line vector points =====#
    
    if hi[0] >= 0 and hi[1] >= 0 and hi[2]>=0:
        pass #Do Nothing
    elif hi[0] < 0 and hi[1] >= 0 and hi[2] >= 0:
        pass #Do Nothing 
        #OCTANT 3
    elif hi[0] < 0 and hi[1] < 0 and hi[2] >= 0:
        theta1 = -theta1;
        #OCTANT 4
    elif hi[0] >= 0 and hi[1] < 0 and hi[2] >= 0:
        theta1 = -theta1;
        #OCTANT 5
    elif hi[0] >= 0 and hi[1] >= 0 and hi[2] < 0:
        theta2 = -theta2;
        #OCTANT 6
    elif hi[0] < 0 and hi[1] >= 0 and hi[2] < 0:
        theta2 = -theta2;
        #OCTANT 7
    elif hi[0] < 0 and hi[1] < 0 and hi[2] < 0:
        theta1 = -theta1;
        theta2 = -theta2;
        #OCTANT 8
    elif hi[0] >= 0 and hi[1] < 0 and hi[2] < 0:
        theta1 = -theta1;
        theta2 = -theta2;
    else:
        pass 
    
    #====== Calculate a point at the end of surface normals =====#
    Px = np.zeros(len(v1))
    Py = np.zeros(len(v1))
    Pz = np.zeros(len(v1))
    for normlcv in range(len(v1)):
        rv1i =v1[normlcv,:] - world
        
        
        #Equation of line in direction of normal vector, staring at
        #vertex 1; Equivalently the position vector from the origin
        #to point P on a line passing through the normal vector and
        #vertex 1
        rpi = rv1i + normalStore[normlcv,:]
        
        Px[normlcv] = world[0] + rpi[0]
        Py[normlcv] = world[1] + rpi[1]
        Pz[normlcv] = world[2] + rpi[2]
    
    
    #======================= Perform a rotation ===========================#
    
    thetaH = deflec[plcv]
    
    #Rotation Matrix 1
    L1 = [[math.cos(theta1*np.pi/180),math.sin(theta1*np.pi/180),0],
           [-1*math.sin(theta1*np.pi/180),math.cos(theta1*np.pi/180),0],
           [0,0,1]]
    L1 = np.array(L1)
    #Rotation Matrix 2
    L2 = [[math.cos(theta2*np.pi/180),0,math.sin(theta2*np.pi/180)],
           [0,1,0],
           [-1*math.sin(theta2*np.pi/180),0,math.cos(theta2*np.pi/180)]]
    L2 = np.array(L2)
    #Rotation Matrix 3       
    L3 = [[1,0,0],
          [0,math.cos(thetaH*np.pi/180),math.sin(thetaH*np.pi/180)],
          [0,-1*math.sin(thetaH*np.pi/180),math.cos(thetaH*np.pi/180)]]
    L3 = np.array(L3)
    #Execute Rotation of vertex points

    newVertices = np.zeros((len(vertices),len(vertices[0])))

    for rowlcv in range(len(vertices)):
        for collcv in range(len(vertices[0])):
            if collcv==0:
                x=vertices[rowlcv,collcv]
            elif collcv == 1:
                y=vertices[rowlcv,collcv]
            elif collcv == 2:
                z=vertices[rowlcv,collcv]
    
        Xb = [[x],[y],[z]]
        Xb = np.array(Xb)
        ## Rotation 
        newVertices[rowlcv,:] = np.transpose(np.linalg.inv(L1)@np.linalg.inv(L2)@L3@L2@L1@(Xb-delta) + delta)
    
    normal = np.zeros((len(v1),len(v1[0]))) 
    for rotnlcv in range(len(v1)):
        XbP = [[Px[rotnlcv]],[Py[rotnlcv]],[Pz[rotnlcv]]] 
        XbP = np.array(XbP) 
        Xbv1 = [[v1[rotnlcv,0]],[v1[rotnlcv,1]],[v1[rotnlcv,2]]]
        Xbv1 = np.array(Xbv1)
        ## Rotation 
        newP = np.linalg.inv(L1)@np.linalg.inv(L2)@L3@L2@L1@(XbP-delta) + delta
        newv1= np.linalg.inv(L1)@np.linalg.inv(L2)@L3@L2@L1@(Xbv1-delta)+ delta
    
        newnorm = newP-newv1
        
        unitnormal = newnorm/np.linalg.norm(newnorm)
        normal[rotnlcv,:] = np.transpose(unitnormal)
        
    #======================= Create new STL file ===========================#
    ## Create file path
    if filename.endswith('.STL'):
        filename = filename[:-4]
    elif filename.endswith('.stl'):
        filename = filename[:-4]
    newfilename = filename+str(thetaH)+'.stl'
    outdir = os.path.join(localpath,'Rotated_Files(Output)')
    filedir = os.path.join(outdir, newfilename)
    
    ## Create and write to file
    outfile = open(filedir,'w')
    
    title = ['solid',filename,'\n']
    title = " ".join(title)
    outfile.write(title)
    normal = normal.tolist()
    newVertices = newVertices.tolist() 

    vertcount = 0
        
    
    for i in range(nfacet):
        normalstring = ' '.join(map(str,normal[i]))
        outfile.write('     facet normal %s\n'% normalstring)
        outfile.write('          outer loop\n')
        for clcv in range(3):
            newVertstring= ' '.join(map(str,newVertices[vertcount]))
            outfile.write('               vertex %s\n' % newVertstring)
            if clcv == 2 and i == len(normalStore):
                pass
            else:
                vertcount = vertcount +1 #to write vertices in groups of 3
        outfile.write('          endloop \n')
        outfile.write('     endfacet \n')        
        
                
    outfile.write('endsolid')
    
    outfile.close() 
                   
                   
                
    
 