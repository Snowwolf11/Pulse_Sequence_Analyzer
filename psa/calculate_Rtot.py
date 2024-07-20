import numpy as np

def _calculate_Rtot(PS, T, Umax, offset, inpoFact):
    Rtot = np.eye(3)  # Rotation Matrix is created
    
    for i,v in enumerate(PS): 
        Ux = (v[0]*inpoFact/100)*Umax*np.cos(v[1]*(np.pi/180))    #x-partiotion of force
        Uy = (v[0]*inpoFact/100)*Umax*np.sin(v[1]*(np.pi/180))    #y-partiotion of force
        Uz=offset                                                  #z-partiotion of force
        n=np.array([Ux, Uy, Uz])                                            #Einheitsvector = rotation-axis is created (= direction of force)
        alpha = -2*np.pi*np.linalg.norm(n)*T/(inpoFact)  #angle of rotation  (how far is rotated)  3.6?
        cosa=np.cos(alpha);
        sina=np.sin(alpha);
        mcosa=1-cosa;
        if np.linalg.norm(n) != 0:
            n=n/np.linalg.norm(n)
            Rn=np.array([[n[0]**2*mcosa+cosa, n[0]*n[1]*mcosa-n[2]*sina, n[0]*n[2]*mcosa+n[1]*sina],      #current Rotation Matrix
                         [n[0]*n[1]*mcosa+n[2]*sina,  n[1]**2*mcosa+cosa, n[1]*n[2]*mcosa-n[0]*sina],
                         [n[2]*n[0]*mcosa-n[1]*sina,  n[2]*n[1]*mcosa+n[0]*sina, n[2]**2*mcosa+cosa]])
              
        else: 
            Rn=np.eye(3)
              #Total rotation Matrix for the i'th vector is created
        Rtot = np.dot(Rtot, Rn)  # Total rotation Matrix for the i'th vector is created
    
    return Rtot
