
import numpy as np
import scipy.special as spe
import scipy.sparse as sps



def KnotInsertionCoefficients(U,u,p):
    """ Returns the knot insertion coefficients """
    k = findKnotSpan(u, U, p)
    l = len(U)-p-1 
    alpha = np.zeros(l+1)
    alpha[:k-p+1] = 1 
    alpha[k-p+1:k+1] = (u-U[k-p+1:k+1])/(U[k+1:k+p+1]-U[k-p+1:k+1])
    alpha[k+1:]  =  0 
    return alpha, k  
def KnotInsertionOperatorOneKnot(U,u,p):
    """ Consttucts the uni-variate knot refinement operator for one knot insertion """
    alpha, k = KnotInsertionCoefficients(U,u,p)
    nbf = len(alpha)-1
    indexI     = np.kron(np.arange(nbf),np.ones(2))
    indexJ     = np.c_[np.arange(nbf),np.arange(1,nbf+1)].ravel() 
    nnz_values = np.c_[alpha[:-1],1-alpha[1:]].ravel() 
    M = sps.csc_matrix(( nnz_values, (indexI,indexJ)), shape = (nbf,nbf+1))
    return M,k
def KnotInsertionOperatorMultipleKnots(U,um,p):
    """ Constructs the uni-variate knot refinement operator for multiple knot insertion""" 
    if len(um)==0:
        raise ValueError('Must insert at least one knot')
    Uc = np.copy(U)
    M,k = KnotInsertionOperatorOneKnot(Uc,um[0],p)
    Uc = np.insert(Uc,k+1,um[0])
    if len(um>1):
        for i in range(1,len(um)):
            Mo,k = KnotInsertionOperatorOneKnot(Uc,um[i],p)
            M = M.dot(Mo)
            Uc = np.insert(Uc,k+1,um[i])
    return M 
def KnotInsertionOperator2d(p,U,um,V,vm):
    """
    returns the 2d refinement opertor  by the tensor product operation
    Coarse knot vectors in Xi and Eta directions: cXi,cEta 
    Fine knot vectors in Xi and Eta directions: rXi,rEta
    """
    Cu  = KnotInsertionOperatorMultipleKnots(U,um,p)
    Cv  = KnotInsertionOperatorMultipleKnots(V,vm,p)

    C = sps.kron(Cv,Cu) 
    return C 

def bezier_extraction_nurbs_1d(U,m,p):
    C    = []
    a    = p+1
    b    = a+1
    nb   = 1
    C.append(np.identity(p+1))
    while b < m:
        C.append(np.identity(p+1)) # Initialize the next extraction operator.
        i = b
        # Count multiplicity of the knot at location b.
        while b < m and U[b] == U[b-1]:
            b = b+1
        mult = b-i+1
        if mult < p:
            #Use (10) to compute the alphas.
            numer = U[b-1]-U[a-1]
            alphas = [0 for x in range(p)]
            for j in range(p,mult,-1):
                alphas[j-mult-1] = numer / (U[a+j-1]-U[a-1])
            r = p-mult
            #Update the matrix coefficients for r new knots
            for j in range(1,r+1):
                save = r-j+1
                s = mult+j
                for k in range(p+1,s,-1):
                    alpha = alphas[k-s-1]
                    #The following line corresponds to (9).
                    C[-2][:,k-1] = alpha*C[-2][:,k-1] + (1.0-alpha)*C[-2][:,k-2]
                if b < m:
                    #Update overlapping coefficients of the next operator.
                    C[-1][save-1:j+save,save-1] = C[-2][p-j:p+1,p]
            nb = nb + 1 # Finished with the current operator.
            if b < m:
                #Update indices for the next operator.
                a = b
                b = b+1
    C.pop()
    return C, nb          
    
    
# def Oslo1(p, coarsekn, finekn, rf ):  
#     """ Knot insertion coefficients 
#     using Oslo algorithm 
#     Taken from Article Multi level Bezier extraction
#     for hierarchical local refinement of IGA """
#     # rf in range(m-p-1)
#     cf = findKnotSpan(finekn[rf],coarsekn,p) 
#     b = 1 
#     for k in range(1,p+1):
#         t1 = coarsekn[cf-k+1:cf+1]
#         t2 = coarsekn[cf+1:cf+k+1]
#         x = finekn[rf+k]
#         w = (x-t1)/(t2-t1)
#         b = np.r_[(1-w)*b,0] + np.r_[0,w*b]
#     return b 


# def KnotInsertionOperator(p,coarsekn,finekn):
#     m = len(finekn)
#     rf = 0 
#     C = Oslo1(p,coarsekn,finekn,rf)
#     for rf in range(1,m-p-1):
#         b = Oslo1(p,coarsekn,finekn,rf)
#         C = np.c_[C,b]
#     return C 
# def KnotInsertionOperator2d(p,cXi,cEta,fXi,fEta):
#     """
#     returns the 2d refinement opertor 
#     Coarse knot vectors in Xi and Eta directions: cXi,cEta 
#     Fine knot vectors in Xi and Eta directions: rXi,rEta
#     """
#     Cxi  = KnotInsertionOperator(p,cXi,fXi)
#     Ceta = KnotInsertionOperator(p,cEta,fEta)

#     C = np.kron(Ceta,Cxi) 
#     return C 

#def MultipleKnotInsertionsOperator(p,coarsekn,knots):
#    coarseTemp = coarsekn.copy()  
#    k = knots[0]
#    cf = findKnotSpan(k,coarseTemp,p) 
#    finekn = np.insert(coarseTemp, cf+1, k)   
#    C = OneKnotInsertionOperator(p,coarseTemp,finekn)
#    coarseTemp = finekn 
#    print(C)
#    for i in range(len(knots)-1):
#            k = knots[i+1]
#            cf = findKnotSpan(k,coarseTemp,p) 
#            finekn = np.insert(coarseTemp, cf+1, k)   
#            print(k,cf,finekn, coarseTemp)
#            print(OneKnotInsertionOperator(p,coarseTemp,finekn))
#            C = C @ OneKnotInsertionOperator(p,coarseTemp,finekn)
#            coarseTemp = finekn 
#    return C 


def BernsteinRef(u,p):
    """ Evaluates Bernsein basis functions defined on [-1,1]"""
    i = np.arange(p+1)
    ur = np.repeat(u,p+1).reshape((len(u),p+1))
    N = ((1-ur)**(p-i)*(1+ur)**i)*spe.comb(p,i)/(2**p) 
    N1 = ((1-ur)**(p-1-(i-1))*(1+ur)**(i-1))*spe.comb(p-1,i-1)/(2**(p-1)) 
    N2 = ((1-ur)**(p-1-(i))*(1+ur)**(i))*spe.comb(p-1,i)/(2**(p-1)) 
    dN = 0.5*p*(N1-N2)
    return N,dN  
 
    

def nubsconnect(p,n):
    """
    It returns the connectivity table INE for NURBS 
    nel : element index in the total connectivity 
    nen : number of functions per support 
    p : degree of functions 
    n : number of functions = number of control points 
    
    fonction qui donne la table de connectivite INE pour un NURBS
    Numerotation en arriere pour trouver directement NURBS coordinate
    nel nombre d'element total
    nen nombre de fonctions par support
    p degre des fnc
    n nombre de fonctions (nbre de points de controle) """
    nel = n-p
    nen = p+1
    IEN = np.zeros((nen,nel),dtype='int64')
    for i in range(nen):
        for j in range(nel):
            IEN[i,j]=p+j-i
    return IEN 

def connect_2D(n_xi,n_eta,p,q):
    """
    Connectivity between 2D elements and control points
    Each parametric element is supported by (p+1)*(q+1) basis function
    """
    nbf_elem = (p+1)*(q+1)
    n_elems = n_xi*n_eta
    IEN = np.zeros((n_elems,nbf_elem))
    k=0
    for j in range(n_eta):
        for i in range(n_xi):
            t=0
            # Local functions numbers support (p+1)*(q+1)
            for s in range(q+1):
                for r in range(p+1):
                    IEN[k,t] = (n_xi+p)*(j+s)+i+r
                    t+=1   
            k+=1
    return IEN
         
#%% 

def GaussLegendre(n):
#     [nodes,weigths]=GaussLegendre(n)
#
# Generates the abscissa and weights for a Gauss-Legendre quadrature.
# Reference:  Numerical Recipes in Fortran 77, Cornell press.

    xg = np.zeros(n)                                           # Preallocations.
    wg = np.zeros(n)   
    m = (n+1)/2
    #import pdb; pdb.set_trace()
    for ii in range(int(m)): # for ii=1:m
        z = np.cos(np.pi*(ii+1-0.25)/(n+0.5))                        # Initial estimate.
        z1 = z+1
        while np.abs(z-z1)>np.finfo(np.float64).eps:
            p1 = 1
            p2 = 0
            for jj in range(n): #for jj = 1:n
                p3 = p2
                p2 = p1
                p1 = ((2*jj+1)*z*p2-(jj)*p3)/(jj+1)      # The Legendre polynomial.
            
            pp = n*(z*p1-p2)/(z**2-1)                        # The L.P. derivative.
            z1 = z
            z = z1-p1/pp
        
        xg[ii] = -z                                   # Build up the abscissas.
        xg[-1-ii] = z
        wg[ii] = 2/((1-z**2)*(pp**2))                     # Build up the weights.
        wg[-1-ii] = wg[ii]
    
    
    return xg,wg

def GaussTetrahedron(n):
    """
    Gauss integration on a tetrahedron using barycentric coordinates 
    returns the barycentric coordinates of the integration points 
    and the weight coefficient:
    Attention: we must multiply the returned weight by the volume of the integrated tetrahedron 
    Taken from Alexander Ern and Jean-Luc Guermond book on finite elements 
    return c: the barycentric coordinates 
           wc: weight coefficient that must be multiplied by the tetrahedon volume 
    """
    if n==1:
        # Deg = 1: 1 point 
        c = np.array([[1/4,1/4,1/4,1/4]])
        wc = np.array([1])
    elif n==2:
        # Deg =2: 4 points
        a = (5-np.sqrt(5))/20
        c = np.array([[1-3*a,a,a,a],
                      [a,1-3*a,a,a],
                      [a,a,1-3*a,a],
                      [a,a,a,1-3*a]])
        wc = np.array([1/4,1/4,1/4,1/4])
    elif n==3:
        # Deg 3: 5 points 
        c = np.array([[1/4,1/4,1/4,1/4],
                      [1/2,1/6,1/6,1/6],
                      [1/6,1/2,1/6,1/6],
                      [1/6,1/6,1/2,1/6],
                      [1/6,1/6,1/6,1/2]])
        wc = np.array([-4/5, 9/20, 9/20, 9/20,9/20])
    elif n == 4 or n== 5 :
        # Deg 5: 15 points 
        a1 = (7-np.sqrt(15))/34
        a2 = (7+np.sqrt(15))/34
        a = (10-2*np.sqrt(5))/40
        w1 = (2665+14*np.sqrt(15))/37800
        w2 = (2665-14*np.sqrt(15))/37800
        w  = 10/189
        c = np.array([[1/4,1/4,1/4,1/4],
                      [1-2*a1,a1,a1,a1],
                      [a1,1-2*a1,a1,a1],
                      [a1,a1,1-2*a1,a1],
                      [a1,a1,a1,1-2*a1],
                      [1-2*a2,a2,a2,a2],
                      [a2,1-2*a2,a2,a2],
                      [a2,a2,1-2*a2,a2],
                      [a2,a2,a2,1-2*a2],
                      [1/2-a,1/2-a,a,a],
                      [1/2-a,a,1/2-a,a],
                      [1/2-a,a,a,1/2-a],
                      [a,1/2-a,1/2-a,a],
                      [a,a,1/2-a,1/2-a],
                      [a,1/2-a,1/2-a,a]])
        wc = np.array([16/135, w1,w1,w1,w1, w2,w2,w2,w2, w,w,w,w,w,w])
    else : 
        raise ValueError('WARNING: Maximum number of Gauss points for tetrahedon reached' )
    
    return c,wc 
def GaussTriangle(n):
    """ Gauss integration on the 
    reference triangle """   
    if n== 1 :
        # Deg = 1 : 1 point 
        x = 1/3*np.array([[1,1]])   
        w = np.array([1/2])
    elif n == 2 :
        # Deg = 2 : 3 points 
        x = 1/6*np.array([[1,1],
                    [4,1],
                    [1,4]])
        w =  1/6*np.array([ 1,1,1 ])
    elif n == 3 : 
        # Deg =3 : 4 points 
        x = np.array([  np.array([1/3 ,1/3]), 
                        1/5*np.array([1, 1]),
                        1/5*np.array([3, 1]),
                        1/5*np.array([1, 3]) ] )  
        w=1/96*np.array([-27 , 25 , 25 , 25]); 
    elif n == 4 : 
        # Deg = 4 : 6 points 
        a=0.445948490915965;
        b=0.091576213509771;
        x = np.array([ [a, a],
                       [1-2*a, a],
                       [a, 1-2*a], 
                       [b, b],
                       [1-2*b, b],
                       [b, 1-2*b] ])
        w = np.r_[  0.111690794839005*np.array([1 , 1 , 1]),  0.054975871827661*np.array([1 , 1 , 1])   ] 
    elif n== 5 : 
        # Deg = 5 :  7 points 
        a=(6+np.sqrt(15))/21;
        b=4/7-a;
        x= np.array([  [1/3, 1/3 ],
                       [a , a], 
                       [1-2*a , a],
                       [a , 1-2*a],
                       [b , b],
                       [1-2*b , b],
                       [b , 1-2*b] ] ) 
        w = np.r_[ np.array([9/80]), (155+np.sqrt(15))/2400*np.array([1 , 1 , 1]), (31/240-(155+np.sqrt(15))/2400)*np.array([1 , 1 , 1])  ]

    elif n == 6 :
        # Deg = 6 : 12 points 
        w = np.array(  [0.05839314 ,
                        0.05839314 ,
                        0.05839314 ,
                        0.02542245 ,
                        0.02542245 ,
                        0.02542245 ,
                        0.04142554 ,
                        0.04142554 ,
                        0.04142554 ,
                        0.04142554 , 
                        0.04142554 ,
                        0.04142554] ) 
          
        x = np.array( [ [0.5014265  ,   0.2492867], 
                 [0.2492867  , 0.5014265] , 
                 [0.2492867  , 0.2492867] , 
                 [0.8738220  , 0.06308901] ,
                 [0.06308901 , 0.8738220] , 
                 [0.06308901 , 0.06308901] , 
                 [0.6365025  , 0.05314505] , 
                 [0.6365025  , 0.3103525] , 
                 [0.05314505 ,  0.6365025] , 
                 [0.05314505  , 0.3103525] , 
                 [0.3103525 ,  0.6365025] , 
                 [0.3103525  , 0.05314505] ] ) 
    elif n == 7 : 
        # Deg 7 : 13 points 
        w= np.array( [-0.07478502 , 
                       0.08780763 , 
                       0.08780763 , 
                       0.08780763 ,
                       0.02667362 ,
                       0.02667362 ,
                       0.02667362 ,
                       0.03855688 ,
                       0.03855688 ,
                       0.03855688 ,
                       0.03855688 ,
                       0.03855688 , 
                       0.03855688 ] )
 
        x= np.array([ [0.3333333,   0.3333333 ] , 
                     [0.4793081  ,  0.2603460],
                     [0.2603460  ,  0.4793081],
                     [0.2603460  ,  0.2603460],
                     [0.8697398  ,  0.06513010],
                     [0.06513010 ,  0.8697398],
                     [0.06513010 ,  0.06513010],
                     [0.6384442  ,  0.04869032],
                     [0.6384442  ,  0.3128655],
                     [0.04869032 ,  0.6384442],
                     [0.04869032 ,  0.3128655],
                     [0.3128655  ,  0.6384442],
                     [0.3128655  ,  0.04869032] ] ) 
    elif n== 8 : 
        # Deg 8 : 16 points 
        w= np.array( [0.07215780 ,
                      0.04754582,
                      0.04754582,
                      0.04754582,
                      0.01622925,
                      0.01622925,
                      0.01622925,
                      0.05160869,
                      0.05160869,
                      0.05160869,
                      0.01361516,
                      0.01361516,
                      0.01361516,
                      0.01361516,
                      0.01361516,
                      0.01361516] ) 
        x= np.array([ [0.3333333,     0.3333333], 
                     [0.08141482 ,   0.4592926 ]  ,
                     [0.4592926  ,   0.08141482 ] ,
                     [0.4592926  ,   0.4592926 ] ,
                     [0.8989055  ,   0.05054723 ] ,
                     [0.05054723 ,   0.8989055 ] ,
                     [0.05054723 ,   0.05054723 ] ,
                     [0.6588614  ,   0.1705693 ] ,
                     [0.1705693  ,   0.6588614 ] ,
                     [0.1705693  ,   0.1705693 ] ,
                     [0.008394777 ,  0.7284924 ] ,
                     [0.008394777 ,  0.2631128 ] ,
                     [0.7284924   ,  0.008394777 ] ,
                     [0.7284924   ,  0.2631128 ] ,
                     [0.2631128   ,  0.008394777 ] , 
                     [0.2631128   ,  0.7284924] ] ) 
    elif n == 9 : 
        # Deg 9 : 19 points 
        w= np.array( [0.04856790 ,
                      0.01566735 ,
                      0.01566735 ,
                      0.01566735,
                      0.03891377,
                      0.03891377,
                      0.03891377,
                      0.03982387,
                      0.03982387,
                      0.03982387,
                      0.01278884,
                      0.01278884,
                      0.01278884,
                      0.02164177,
                      0.02164177,
                      0.02164177,
                      0.02164177,
                      0.02164177,
                      0.02164177] ) 
        x=  np.array([ [0.3333333 ,   0.3333333 ], 
                       [0.02063496,   0.4896825  ],
                       [0.4896825 ,   0.02063496 ],
                       [0.4896825 ,   0.4896825  ],
                       [0.1258208 ,   0.4370896  ],
                       [0.4370896 ,   0.1258208  ],
                       [0.4370896 ,   0.4370896  ],
                       [0.6235929 ,   0.1882035  ],
                       [0.1882035 ,   0.6235929  ],
                       [0.1882035 ,   0.1882035  ],
                       [0.9105410 ,   0.04472951 ],
                       [0.04472951 ,  0.9105410  ],
                       [0.04472951 ,  0.04472951 ],
                       [0.03683841 ,  0.7411986  ],
                       [0.03683841 ,  0.2219630  ],
                       [0.7411986  ,  0.03683841 ],
                       [0.7411986  ,  0.2219630  ],
                       [0.2219630  ,  0.03683841 ],
                       [0.2219630  ,  0.7411986  ] ] )  
    elif n == 11 : 
        # Deg 11 : 27 points 
        w= np.array( [0.006829866 ,
            0.006829866 ,
            0.006829866 ,
            0.01809227 ,
            0.01809227 ,
            0.01809227 ,
            0.0004635032 ,
            0.0004635032 ,
            0.0004635032 ,
            0.02966149 ,
            0.02966149 ,
            0.02966149 ,
            0.03857477 ,
            0.03857477 ,
            0.03857477 ,
            0.02616856 ,
            0.02616856 ,
            0.02616856 ,
            0.02616856 ,
            0.02616856 ,
            0.02616856 ,
            0.01035383 ,
            0.01035383 ,
            0.01035383 ,
            0.01035383 ,
            0.01035383 ,
            0.01035383] ) 
        x= np.array([ [0.9352701    ,      0.03236495 ], 
                     [0.03236495 ,  0.9352701  ]
                     [0.03236495 ,  0.03236495 ]
                     [0.7612982  , 0.1193509  ]
                     [0.1193509  , 0.7612982  ]
                     [0.1193509  , 0.1193509  ]
                     [0.06922210 ,  0.5346110 ]
                     [0.5346110  , 0.06922210 ]
                     [0.5346110  , 0.5346110  ]
                     [0.5933802  , 0.2033099  ]
                     [0.2033099  , 0.5933802  ]
                     [0.2033099  , 0.2033099  ]
                     [0.2020614  , 0.3989693  ] 
                     [0.3989693  , 0.2020614  ]
                     [0.3989693  , 0.3989693  ]
                     [0.05017814 ,  0.5932012 ]
                     [0.05017814 ,  0.3566206 ]
                     [0.5932012  , 0.05017814 ]
                     [0.5932012  , 0.3566206  ]
                     [0.3566206  , 0.05017814 ]
                     [0.3566206  , 0.5932012  ]
                     [0.02102202 ,  0.8074890 ]
                     [0.02102202 ,  0.1714890 ]
                     [0.8074890  , 0.02102202 ]
                     [0.8074890  , 0.1714890  ]
                     [0.1714890  , 0.02102202 ] 
                     [0.1714890  , 0.8074890] ] )  
    else : 
        raise ValueError('WARNING: Maximum number of Gauss points for triangle reached (11)' ) 
    return x,w 


# def derbasisfuncVectorInput(p,U,u,nb_u_values,span,nders):
#     ders_matrix = np.zeros(((nders+1)*(nb_u_values),p+1))
#     for i in range(nb_u_values):
#         ders = derbasisfuns(span,p,U,nders,u[i])
#         ders_matrix[2*i:2*(i+1),:] = ders  
#     return ders_matrix 
# def derbasisfunsPython(p,U,u,nb_u_values,span,nders):
#     ders = derbasisfuncVectorInput(p,U,u,nb_u_values,span,nders)
#     N = ders[::2,:]
#     dN = ders[1::2,:]
#     return N,dN  
    

def Lagrange1D(n, xi):
    """1D Lagrange basis functions 
    # defined on the reference element [-1,1]"""
    if n == 2 :
        N = np.array([0.5*(1-xi),0.5*(xi+1)]) 
        dN = np.array([-0.5,0.5]) 
    else :
        raise ValueError('Only n=2 implemented')
    return [N,dN]

def Lagrange2D(n,type_e, xi, eta):
    """2D Lagrange basis functions 
    defined on reference elements ([-1,1] square or [0,1] triangle)
    n: number of basis functions 
    type_e: element type quad or triangle
    """
    if type_e !='quad' and type_e!='triangle':
        raise ValueError('Element type not considered')
    if type_e =='triangle':
        if n == 3:
            N = np.array([1-xi-eta,eta,xi])
            dNdxi = np.array([-1,0,1])
            dNdeta = np.array([-1,1,0])
        else : 
            raise ValueError('Only T3 implemented')
    elif type_e =='quad':
        if n== 4:
            N = np.array([0.25*(1-xi)*(1-eta),
                      0.25*(1-xi)*(1+eta),
                      0.25*(1+xi)*(1+eta),
                      0.25*(1+xi)*(1-eta)])
            dNdxi = np.array([-0.25*(1-eta),
                          -0.25*(1+eta),
                           0.25*(1+eta),
                           0.25*(1-eta)])
            dNdeta = np.array([-0.25*(1-xi),
                            0.25*(1-xi),
                            0.25*(1+xi),
                           -0.25*(1+xi)]) 
        else :
            raise ValueError('Only Q4 implemented')  
    return [N,dNdxi,dNdeta]




#%%
def findKnotSpan(u,U,p):
    """
    Finds the knots space of a given knot parameter u in
    the knot vector U corresponding to the degree p 
    """
    m = np.size(U)
    if u==U[m-p-1]:
        k=m-p-2
    else :
        k=np.max(np.where(u>=U)) 
    return k 

 

def findspan(n,p,u,U):               
    return findKnotSpan(u,U,p)                       
                     
def findspanUniformKnotVector(U,deg,l,u):
    if u==U[len(U)-deg-1]:
        return len(U)-deg-2
    return int(np.floor( (u-U[0])/l))+deg 
 

#%% Semble Ok. Bien penser à mettre les knots à ajouter sous forme d'un vecteur np.array (même s'il n'y en a qu'un)
def bspkntins(d,c,k,u):
    ''' Function Name: 
    #  
    #   bspkntins - Insert knots into a univariate B-Spline. 
    #  
    # Calling Sequence: 
    #  
    #   [ic,ik] = bspkntins(d,c,k,u) 
    #  
    # Parameters: 
    #  
    #   d	: Degree of the B-Spline. 
    #  
    #   c	: Control points, matrix of size (dim,nc). 
    #  
    #   k	: Knot sequence, row vector of size nk. 
    #  
    #   u	: Row vector of knots to be inserted, size nu 
    #  
    #   ic	: Control points of the new B-Spline, of size (dim,nc+nu) 
    #  
    #   ik	: Knot vector of the new B-Spline, of size (nk+nu) 
    #  
    # Description: 
    #  
    #   Insert knots into a univariate B-Spline. This function provides an 
    #   interface to a toolbox 'C' routine. '''
    mc,nc = c.shape
    nu = len(u) 
    nk = len(k) 
                                                         #  
                                                         # int bspkntins(int d, double *c, int mc, int nc, double *k, int nk, 
                                                         #               double *u, int nu, double *ic, double *ik) 
                                                         # { 
                                                         #   int ierr = 0; 
                                                         #   int a, b, r, l, i, j, m, n, s, q, ind; 
                                                         #   double alfa; 
                                                         # 
                                                         #   double **ctrl  = vec2mat(c, mc, nc); 
    ic = np.zeros((mc,nc+nu))                            #   double **ictrl = vec2mat(ic, mc, nc+nu); 
    ik = np.zeros(nk+nu) 
                                                         # 
    n = c.shape[1] - 1                                   #   n = nc - 1; 
    r = len(u) - 1                                       #   r = nu - 1; 
                                                         # 
    m = n + d + 1                                        #   m = n + d + 1; 
    a = findspan(n, d, u[0], k)                          #   a = findspan(n, d, u[0], k); 
    b = findspan(n, d, u[r], k)                          #   b = findspan(n, d, u[r], k); 
    b+=1                                                 #   ++b; 
                                                         # 
    for q in range(mc):                                  #   for (q = 0; q < mc; q++)  { 
       for j in range(a-d+1):
           ic[q,j] = c[q,j]                              #     for (j = 0; j <= a-d; j++) ictrl[j][q] = ctrl[j][q]; 
       for j in range(b-1,n+1):
           ic[q,j+r+1] = c[q,j]                          #     for (j = b-1; j <= n; j++) ictrl[j+r+1][q] = ctrl[j][q]; 
                                                         #   } 
    
    for j in range(a+1):
        ik[j] = k[j]                                     #   for (j = 0; j <= a; j++)   ik[j] = k[j]; 
    for j in range(b+d,m+1):
        ik[j+r+1] = k[j]                                 #   for (j = b+d; j <= m; j++) ik[j+r+1] = k[j]; 
                                                         # 
    i = b + d - 1                                        #   i = b + d - 1; 
    s = b + d + r                                        #   s = b + d + r; 
     
    for j in range(r,-1,-1):                             #   for (j = r; j >= 0; j--) { 
       while (u[j] <= k[i] and i > a):                   #     while (u[j] <= k[i] && i > a) { 
           for q in range(mc):                           #       for (q = 0; q < mc; q++) 
               ic[q,s-d-1] = c[q,i-d-1]                  #         ictrl[s-d-1][q] = ctrl[i-d-1][q]; 
                                                    
           ik[s] = k[i]                                  #       ik[s] = k[i]; 
           s -= 1                                        #       --s; 
           i -= 1                                        #       --i; 
                                                         #     } 
        
       for q in range(mc):                               #     for (q = 0; q < mc; q++) 
           ic[q,s-d-1] = ic[q,s-d]                       #       ictrl[s-d-1][q] = ictrl[s-d][q]; 
     
       for l in range(1,d+1):                            #     for (l = 1; l <= d; l++)  { 
           ind = s - d + l                               #       ind = s - d + l; 
           alfa = ik[s+l] - u[j]                         #       alfa = ik[s+l] - u[j]; 
           if abs(alfa) == 0:                            #       if (fabs(alfa) == 0.0) 
               for q in range(mc):                       #         for (q = 0; q < mc; q++) 
                   ic[q,ind-1] = ic[q,ind]               #           ictrl[ind-1][q] = ictrl[ind][q]; 
           else:                                         #       else  { 
               alfa = alfa/(ik[s+l] - k[i-d+l])          #         alfa /= (ik[s+l] - k[i-d+l]); 
               for q in range(mc):                       #         for (q = 0; q < mc; q++) 
                   tmp = (1.-alfa)*ic[q,ind] 
                   ic[q,ind-1] = alfa*ic[q,ind-1] + tmp  #           ictrl[ind-1][q] = alfa*ictrl[ind-1][q]+(1.0-alfa)*ictrl[ind][q]; 
                                                         #       } 
                                                         #     } 
       # 
       ik[s] = u[j]                                      #     ik[s] = u[j]; 
       s -= 1                                            #     --s; 
                                                         #   } 
                                                         # 
                                                         #   freevec2mat(ctrl); 
                                                         #   freevec2mat(ictrl); 
                                                         # 
                                                         #   return ierr; 
                                                         # } 
    return ic,ik
#%%

def bincoeff(n,k):
    #  Computes the binomial coefficient. 
    # 
    #      ( n )      n! 
    #      (   ) = -------- 
    #      ( k )   k!(n-k)! 
    # 
    #  b = bincoeff(n,k) 
    # 
    #  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215. 
     
                                                              # double bincoeff(int n, int k) 
                                                              # { 
    b = np.floor(0.5+np.exp(factln(n)-factln(k)-factln(n-k)));      #   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k))); 
     
    return b
     
def factln(n):
    # computes ln(n!) 
    if n <= 1:
        f = 0
        return f
    f = spe.gammaln(n+1) #log(factorial(n));</pre>
    
    return f
#%%


def bspdegelev(d,c,k,t):
    ''' 
    # Function Name:  
    #   bspdegevel - Degree elevate a univariate B-Spline. 
    # Calling Sequence: 
    #   [ic,ik] = bspdegelev(d,c,k,t) 
    # Parameters: 
    #   d	: Degree of the B-Spline. 
    #   c	: Control points, matrix of size (dim,nc). 
    #   k	: Knot sequence, row vector of size nk. 
    #   t	: Raise the B-Spline degree t times. 
    #   ic	: Control points of the new B-Spline. 
    #   ik	: Knot vector of the new B-Spline. 
    # Description: 
    #   Degree elevate a univariate B-Spline. This function provides an 
    #   interface to a toolbox 'C' routine. 
    '''
    mc,nc = c.shape 
                                                          # 
                                                          # int bspdegelev(int d, double *c, int mc, int nc, double *k, int nk, 
                                                          #                int t, int *nh, double *ic, double *ik) 
                                                          # { 
                                                          #   int row,col 
                                                          # 
                                                          #   int ierr = 0; 
                                                          #   int i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, cind, oldr, mul; 
                                                          #   int n, lbz, rbz, save, tr, kj, first, kind, last, bet, ii; 
                                                          #   double inv, ua, ub, numer, den, alf, gam; 
                                                          #   double **bezalfs, **bpts, **ebpts, **Nextbpts, *alfs; 
                                                          # 
    #init ic                                                      #   double **ctrl  = vec2mat(c, mc, nc); 
    ic = np.zeros((mc,nc*(t+1)))                                  #   double **ictrl = vec2mat(ic, mc, nc*(t+1)); 
    ik = np.zeros((t+1)*k.shape[0])
                                                         # 
    n = nc - 1                                               #   n = nc - 1; 
                                                              # 
    bezalfs = np.zeros((d+1,d+t+1))                              #   bezalfs = matrix(d+1,d+t+1); 
    bpts = np.zeros((mc,d+1))                                     #   bpts = matrix(mc,d+1); 
    ebpts = np.zeros((mc,d+t+1))                                 #   ebpts = matrix(mc,d+t+1); 
    Nextbpts = np.zeros((mc,d+1))                                 #   Nextbpts = matrix(mc,d+1); 
    alfs = np.zeros((d,1))                                        #   alfs = (double *) mxMalloc(d*sizeof(double)); 
                                                              # 
    m = n + d + 1                                            #   m = n + d + 1; 
    ph = d + t                                               #   ph = d + t; 
    ph2 = int(ph/2)                                      #   ph2 = ph / 2; 
                                                              # 
                                                              #   // compute bezier degree elevation coefficeients 
    bezalfs[0,0] = 1.                                         #   bezalfs[0][0] = bezalfs[ph][d] = 1.0; 
    bezalfs[d,ph] = 1.                                   # 
     
    for i in np.arange(1,ph2+1):   #1:ph2                                               #   for (i = 1; i <= ph2; i++) { 
        inv = 1/bincoeff(ph,i)                                #     inv = 1.0 / bincoeff(ph,i); 
        mpi = min(d,i)                                        #     mpi = min(d,i); 
                                                              # 
        for j in np.arange(max(0,i-t),mpi+1):   #max(0,i-t):mpi      #     for (j = max(0,i-t); j <= mpi; j++) 
            bezalfs[j,i] = inv*bincoeff(d,j)*bincoeff(t,i-j)  #       bezalfs[i][j] = inv * bincoeff(d,j) * bincoeff(t,i-j); 
                                                              # 
    for i in np.arange(ph2+1,ph):   #ph2+1:ph-1                                          #   for (i = ph2+1; i <= ph-1; i++) { 
        mpi = min(d,i)                                        #     mpi = min(d, i); 
        for j in np.arange(max(0,i-t),mpi+1):   #max(0,i-t):mpi                                   #     for (j = max(0,i-t); j <= mpi; j++) 
            bezalfs[j,i] = bezalfs[d-j,ph-i]         #       bezalfs[i][j] = bezalfs[ph-i][d-j]; 
                                                              # 
    mh = ph                                                  #   mh = ph;       
    kind = ph+1                                              #   kind = ph+1; 
    r = -1                                                   #   r = -1; 
    a = d                                                    #   a = d; 
    b = d+1                                                  #   b = d+1; 
    cind = 1                                                 #   cind = 1; 
    ua = k[0]                                                #   ua = k[0];  
                                                              # 
    for ii in range(mc):  #0:mc-1                                             #   for (ii = 0; ii < mc; ii++) 
        ic[ii,0] = c[ii,0]                                #     ictrl[0][ii] = ctrl[0][ii]; 
    for i in range(ph+1): #0:ph                                                #   for (i = 0; i <= ph; i++) 
        ik[i] = ua                                          #     ik[i] = ua; 
                                                              #   // initialise first bezier seg 
    for i in range(d+1): #0:d                                                 #   for (i = 0; i <= d; i++) 
        for ii in range(mc): #0:mc-1                                          #     for (ii = 0; ii < mc; ii++) 
            bpts[ii,i] = c[ii,i]                       #       bpts[i][ii] = ctrl[i][ii]; 
                                                              #   // big loop thru knot vector 
    while b < m  :                                             #   while (b < m)  { 
        i = b                                                 #     i = b; 
        while b < m and k[b] == k[b+1]  :                      #     while (b < m && k[b] == k[b+1]) 
            b = b + 1                                          #       b++; 
        mul = b - i + 1                                       #     mul = b - i + 1; 
        mh += mul + t                                     #     mh += mul + t; 
        ub = k[b]                                           #     ub = k[b]; 
        oldr = r                                              #     oldr = r; 
        r = d - mul                                           #     r = d - mul; 
                                                              # 
                                                              #     // insert knot u(b) r times 
        if oldr > 0:                                           #     if (oldr > 0) 
            lbz = np.floor((oldr+2)/2)                            #       lbz = (oldr+2) / 2; 
        else :                                                  #     else 
            lbz = 1                                            #       lbz = 1; 
        
        if r > 0 :                                              #     if (r > 0) 
            rbz = ph - np.floor((r+1)/2)                          #       rbz = ph - (r+1)/2; 
        else :                                                  #     else 
            rbz = ph                                           #       rbz = ph; 
        
        if r > 0 :                                             #     if (r > 0) { 
                                                              #       // insert knot to get bezier segment 
            numer = ub - ua                                    #       numer = ub - ua; 
            for q in np.arange(d,mul,-1):    #d:-1:mul+1                                    #       for (q = d; q > mul; q--) 
                alfs[q-mul-1] = numer / (k[a+q]-ua)             #         alfs[q-mul-1] = numer / (k[a+q]-ua); 
           
            for j in np.arange(1,r+1):   #1:r                                           #       for (j = 1; j <= r; j++)  { 
                save = r - j                                    #         save = r - j; 
                s = mul + j                                     #         s = mul + j; 
                                                              # 
                for q in np.arange(d,s-1,-1): #d:-1:s                                     #         for (q = d; q >= s; q--) 
                    for ii in range(mc): #0:mc-1                                 #           for (ii = 0; ii < mc; ii++) 
                        tmp1 = alfs[q-s]*bpts[ii,q]  
                        tmp2 = (1-alfs(q-s))*bpts(ii,q-1)  
                        bpts[ii,q] = tmp1 + tmp2              #             bpts[q][ii] = alfs[q-s]*bpts[q][ii]+(1.0-alfs[q-s])*bpts[q-1][ii]; 
              
                for ii in range(mc): #0:mc-1                                    #         for (ii = 0; ii < mc; ii++) 
                    Nextbpts[ii,save] = bpts[ii,d]       #           Nextbpts[save][ii] = bpts[d][ii]; 
                                                              #     // end of insert knot 
                                                              # 
                                                              #     // degree elevate bezier 
        for i in np.arange(lbz,ph+1):  #lbz:ph                                           #     for (i = lbz; i <= ph; i++)  { 
            for ii in range(mc): #0:mc-1                                       #       for (ii = 0; ii < mc; ii++) 
                ebpts[ii,i] = 0                             #         ebpts[i][ii] = 0.0; 
            mpi = min(d, i)                                    #       mpi = min(d, i); 
            for j in np.arange(max(0,i-t),mpi+1): #max(0,i-t):mpi                                #       for (j = max(0,i-t); j <= mpi; j++) 
                for ii in range(mc): #0:mc-1                                    #         for (ii = 0; ii < mc; ii++) 
                    tmp1 = ebpts[ii,i]  
                    tmp2 = bezalfs[j,i]*bpts[ii,j] 
                    ebpts[ii,i] = tmp1 + tmp2                #           ebpts[i][ii] = ebpts[i][ii] + bezalfs[i][j]*bpts[j][ii]; 
                                                              #     // end of degree elevating bezier 
                                                              # 
        if oldr > 1 :                                           #     if (oldr > 1)  { 
                                                              #       // must remove knot u=k[a] oldr times 
            first = kind - 2                                                    #       first = kind - 2; 
            last = kind                                        #       last = kind; 
            den = ub - ua                                      #       den = ub - ua; 
            bet = np.floor((ub-ik[kind-1]) / den)                   #       bet = (ub-ik[kind-1]) / den; 
                                                              # 
                                                              #       // knot removal loop 
            for tr in np.arange(1,oldr):  #1:oldr-1                                     #       for (tr = 1; tr < oldr; tr++)  { 
                i = first                                       #         i = first; 
                j = last                                        #         j = last; 
                kj = j - kind + 1                               #         kj = j - kind + 1; 
                while j-i > tr :                                  #         while (j - i > tr)  { 
                                                              #           // loop and compute the new control points 
                                                              #           // for one removal step 
                    if i < cind  :                                 #           if (i < cind)  { 
                        alf = (ub-ik[i])/(ua-ik[i])           #             alf = (ub-ik[i])/(ua-ik[i]); 
                        for ii in range(mc):  #0:mc-1                              #             for (ii = 0; ii < mc; ii++) 
                            tmp1 = alf*ic[ii,i] 
                            tmp2 = (1-alf)*ic[ii,i-1]  
                            ic[ii,i] = tmp1 + tmp2             #               ictrl[i][ii] = alf * ictrl[i][ii] + (1.0-alf) * ictrl[i-1][ii]; 
                    if j >= lbz :                                   #           if (j >= lbz)  { 
                        if j-tr <= kind-ph+oldr :                   #             if (j-tr <= kind-ph+oldr) { 
                            gam = (ub-ik[j-tr]) / den            #               gam = (ub-ik[j-tr]) / den; 
                            for ii in range(mc):  #0:mc-1                           #               for (ii = 0; ii < mc; ii++) 
                                tmp1 = gam*ebpts[ii,kj]  
                                tmp2 = (1-gam)*ebpts[ii,kj+1]  
                                ebpts[ii,kj] = tmp1 + tmp2      #                 ebpts[kj][ii] = gam*ebpts[kj][ii] + (1.0-gam)*ebpts[kj+1][ii]; 
                        else :                                      #             else  { 
                            for ii in range(mc):  #0:mc-1                           #               for (ii = 0; ii < mc; ii++) 
                                tmp1 = bet*ebpts[ii,kj]                                      
                                tmp2 = (1-bet)*ebpts[ii,kj+1]                                      
                                ebpts[ii,kj] = tmp1 + tmp2      #                 ebpts[kj][ii] = bet*ebpts[kj][ii] + (1.0-bet)*ebpts[kj+1][ii]; 
                    i += 1                                    #           i++; 
                    j -= 1                                    #           j--; 
                    kj -= 1                                  #           kj--; 
                                                              # 
                first -= 1                               #         first--; 
                last += 1                                 #         last++; 
                                                              #     // end of removing knot n=k[a] 
                                                              # 
                                                              #     // load the knot ua 
        if a != d :                                             #     if (a != d) 
            for i in range(ph-oldr):   #0:ph-oldr-1                                   #       for (i = 0; i < ph-oldr; i++)  { 
                ik[kind] = ua                                 #         ik[kind] = ua; 
                kind += 1                                 #         kind++; 
                                                          # 
                                                          #     // load ctrl pts into ic 
        for j in np.arange(lbz,rbz+1):   #lbz:rbz                                       #     for (j = lbz; j <= rbz; j++)  { 
            for ii in range(mc):   #0:mc-1                                    #       for (ii = 0; ii < mc; ii++) 
                ic[ii,cind] = ebpts[ii,j]            #         ictrl[cind][ii] = ebpts[j][ii]; 
            cind += 1                                 #       cind++; 
                                                          # 
        if b < m :                                           #     if (b < m)  { 
                                                          #       // setup for next pass thru loop 
            for j in range(r): #0:r-1                                      #       for (j = 0; j < r; j++) 
                for ii in range(mc): #0:mc-1                                 #         for (ii = 0; ii < mc; ii++) 
                    bpts[ii,j] = Nextbpts[ii,j]       #           bpts[j][ii] = Nextbpts[j][ii]; 
            for j in np.arange(r,d+1):  #r:d                                        #       for (j = r; j <= d; j++) 
                for ii in range(mc):  #0:mc-1                                 #         for (ii = 0; ii < mc; ii++) 
                    bpts[ii,j] = c[ii,b-d+j]          #           bpts[j][ii] = ctrl[b-d+j][ii]; 
            a = b                                           #       a = b; 
            b += 1                                         #       b++; 
            ua = ub                                         #       ua = ub; 
                                                          #     } 
        else:                                                #     else 
                                                      #       // end knot 
            for i in range(ph+1):   #0:ph                                       #       for (i = 0; i <= ph; i++) 
                ik[kind+i] = ub                            #         ik[kind+i] = ub; 
    # End big while loop                                      #   // end while loop 
                                                              # 
                                                              #   *nh = mh - ph - 1; 
                                                              # 
                                                              #   freevec2mat(ctrl); 
                                                              #   freevec2mat(ictrl); 
                                                              #   freematrix(bezalfs); 
                                                              #   freematrix(bpts); 
                                                              #   freematrix(ebpts); 
                                                              #   freematrix(Nextbpts); 
                                                              #   mxFree(alfs); 
                                                              # 
                                                              #   return(ierr); 
                                                              # } 
                                                              
    # ajout dû au fait qu'on a initialisé trop grand (car difficile d'estimer la taille de ic et ik avant, dépend entre autres de la multiplicité des knots)
    # on enleve les 0 à la fin du knot vector ik
    ik = np.trim_zeros(ik,'b')
    
    # on tronque la matrice des points de contrôle où il faut (revient à enlever les 0, mais si la courbe finit avec un point en (0,0), on n'enlève pas celui-là)
    n = len(ik)-(d+t)-1
    ic = ic[:,0:n]

    return ic,ik 





def derbasisfuns(i,pl,U,nders,u):
    """ 
    i : knot span of u 
    p : degree 
    u : parameter on which we want to evaluate the function 
    nders: number of derivatives 
    U : knot vector 
    
    # i span de u 
    # pl = degrés de la nurbs
    # u = endroit ou l'on veut la fonction
    # nders = numéro de la dérivée désirée
    # U = vecteur de noeud de la fonction """
    
#    import pdb; pdb.set_trace()
    u_knotl=U.copy()
    left = np.zeros((pl+1))
    right = np.zeros((pl+1))
    ndu = np.zeros((pl+1,pl+1))
    ders = np.zeros((nders+1,pl+1))
    ndu[0,0] = 1
    for j in range(pl):   #1:pl
        left[j+1] = u - u_knotl[i-j]   ### rq Ali : i-j au lieu de i-j-1
        right[j+1] = u_knotl[i+j+1] - u ### rq : i+j+1 au lieu de i+j
        saved = 0
        for r in range(j+1):   #0:j-1
            ndu[j+1,r] = right[r+1] + left[j-r+1]
            temp = ndu[r,j]/ndu[j+1,r]
            ndu[r,j+1] = saved + right[r+1]*temp
            saved = left[j-r+1]*temp
        ndu[j+1,j+1] = saved
#    print('checkpoint1 : '+str(ndu))
      
                                    # load basis functions
    for j in range(pl+1):   #0:pl
        ders[0,j] = ndu[j,pl]
#    print('checkpoint2 : '+str(ders))

                                # compute derivatives
    for r in range(pl+1): # 0:pl              # loop over function index
        s1 = 0
        s2 = 1                # alternate rows in array a
        a = np.zeros((nders+1,nders+1))
        a[0,0] = 1
                    # loop to compute kth derivative
        for k in range(nders):   # 1:nders
            d = 0
            rk = r-(k+1)
            pk = pl-(k+1)
            if (r >= (k+1)):
                a[s2,0] = a[s1,0]/ndu[pk+1,rk]
                d = a[s2,0]*ndu[rk,pk]
            if (rk >= -1):
                j1 = 1
            else: 
                j1 = -rk
            if ((r-1) <= pk): 
                j2 = k
            else: 
                j2 = pl-r
            for j in np.arange(j1,j2+0.1):   #j1:j2
                j = int(j)
                a[s2,j] = (a[s1,j] - a[s1,j-1])/ndu[pk+1,rk+j]
                d = d + a[s2,j]*ndu[rk+j,pk]
            if (r <= pk): 
                a[s2,k+1] = -a[s1,k]/ndu[pk+1,r]
                d = d + a[s2,k+1]*ndu[r,pk]
            ders[k+1,r] = d
            j = s1
            s1 = s2
            s2 = j            # switch rows

      
    #     Multiply through by the correct factors
     
    r = pl
    for k in range(nders):   # 1:nders
        for j in range(pl+1):   # 0:pl
            ders[k+1,j] = ders[k+1,j]*r
        r = r*(pl-(k+1))

    return ders

def BasisFunc(i,u,p,U):
    """
    Evaluates the non zero basis functions N_{i,p}(u)
    for u living in the knot span [U_i,U_{i+1}[
    From the Nurbs Book ( Les Piegl, Wayne Tiller)
    """
    N=np.zeros(p+1)
    left=np.zeros(p+1)
    right=np.zeros(p+1)
    N[0]=1.
    for j in range(1,p+1):
        left[j]=u-U[i+1-j]
        right[j]=U[i+j]-u
        saved=0.
        for r in range(j):
            temp=N[r]/(right[r+1]+left[j-r])
            N[r]=saved+right[r+1]*temp
            saved=left[j-r]*temp
        N[j]=saved
    return N

def EvaluateSpline(x,y,Xi,Eta,p,q,c,nx):
    """ Evaluate spline surface """
    # c is the raveled array of control variables 
    se = np.zeros(len(x))
    for k in range(len(x)):
        spanx = findKnotSpan(x[k],Xi,p)    
        spany = findKnotSpan(y[k],Eta,q) 
        Nx = BasisFunc(spanx,x[k],p,Xi)
        Ny = BasisFunc(spany,y[k],q,Eta)
        s = 0 
        for j in range(q+1):
            for i in range(p+1):
                s += Ny[j]*Nx[i]*c[spanx-p+i +(spany-q+j)*nx]
        se[k] = s 
    return se
            

def OneBasisFun(i,u,p,UU):
    """
    Computes the basis function N_{i,p} on the point u 
    Warning: Works only for open knot vectors of type 
    U ={0,...,0,  ....  ,  1,...1} where the first 
    and last knots have a multiplicity equal to p+1 
    From the Nurbs Book ( Les Piegl, Wayne Tiller)
    """
    U=np.copy(UU)
    m=U.shape[0]-1
    
    # Verifying the local support proprety of the Bsplines
    if (u < U[i] or u>U[i+p+1]) :
        Nip = 0.                                                         
        return Nip
    # Special case evaluation on the first and last knots 
    if (i==0 and u==U[0]) or (i==m-p-1 and u==U[m]):
        Nip = 1
        return Nip

    N=np.zeros(p+1)
    # Initialize zeroth-degree functions
    for j in range(p+1):
        if u>=U[i+j] and u<U[i+j+1]:
            N[j]=1.
        else :
            N[j]=0.
    for k in range(1,p+1):
        if N[0]==0.:
            saved=0.
        else:
            saved=((u-U[i])*N[0])/(U[i+k]-U[i])
        for j in range(p-k+1):
            Uleft=U[i+j+1]
            Uright=U[i+j+k+1]
            if N[j+1]==0. :
                N[j]=saved 
                saved=0.
            else:
                temp=N[j+1]/(Uright-Uleft)
                N[j]=saved+(Uright-u)*temp
                saved=(u-Uleft)*temp
    Nip=N[0]
    return Nip 

def derbasisfuncVectorInput(p,U,u,nb_u_values,span,nders):
    """ Returns the basis functions and first derivatives 
    for multiple parameter inputs """ 
    ders_matrix = np.zeros(((nders+1)*(nb_u_values),p+1))
    for i in range(nb_u_values):
        ders = derbasisfuns(span,p,U,nders,u[i])
        ders_matrix[2*i:2*(i+1),:] = ders 
    return ders_matrix[::2,:], ders_matrix[1::2,:] 

def secondOrderDerivBasisFuns(p,U,u):
    nb_u_values = len(u);
    nnz_values = nb_u_values*(p+1)
    indexI      = np.zeros(nnz_values)
    indexJ      = np.zeros(nnz_values)
    ders_values = np.zeros(nnz_values)
    n = len(U) -1 -p   
    k = 0 
    for i in range(nb_u_values):
        span = findspan(n,p,u[i],U)
        ders = derbasisfuns(span,p,U,2,u[i])
        for j in range(p+1):
            indexI[k] = i 
            indexJ[k] = span - p +j 
            ders_values[k] = ders[2][j] 
            k+=1 
    phi_second= sps.csc_matrix(( ders_values, (indexI,indexJ)), shape = (nb_u_values,n) )
    return phi_second   

def global_basisfuns(p,U,u):
    """ Returns the global differential operator N and dNdxi
    N and dNdxi are sparse matrices of size (size(u),nbf)
    where nbf is the toral number of basis functions
    """
    nb_u_values = len(u);
    nnz_values = nb_u_values*(p+1)
    indexI      = np.zeros(nnz_values)
    indexJ      = np.zeros(nnz_values)
    values      = np.zeros(nnz_values)
    ders_values = np.zeros(nnz_values)
    n = len(U) -1 -p   
    k = 0 
    for i in range(nb_u_values):
        span = findspan(n,p,u[i],U)
        ders = derbasisfuns(span,p,U,1,u[i])
        for j in range(p+1):
            indexI[k] = i 
            indexJ[k] = span - p +j 
            values[k] = ders[0][j] 
            ders_values[k] = ders[1][j] 
            k+=1 
    phi  = sps.csc_matrix(( values, (indexI,indexJ)), shape = (nb_u_values,n) )
    dphi = sps.csc_matrix(( ders_values, (indexI,indexJ)), shape = (nb_u_values,n) )
    return phi,dphi


def global_basisfunsWd(p,U,u):
    """ Returns only the function N without derivatives""" 
    nb_u_values = len(u);
    nnz_values = nb_u_values*(p+1)
    indexI      = np.zeros(nnz_values)
    indexJ      = np.zeros(nnz_values)
    values      = np.zeros(nnz_values)
    n = len(U) -1 -p   
    k = 0 
    for i in range(nb_u_values):
        span = findspan(n,p,u[i],U)
        N = BasisFunc(span,u[i],p,U)
        for j in range(p+1):
            indexI[k] = i 
            indexJ[k] = span - p +j 
            values[k] = N[j] 
            k+=1 
    phi  = sps.csc_matrix(( values, (indexI,indexJ)), shape = (nb_u_values,n) )
    return phi 

def Get2dBasisFunctionsAtPts(u,v,U,V,p,q):
    nb_u_values = len(u);
    nnz_values  = nb_u_values*(p+1)*(q+1)
    indexI      = np.zeros(nnz_values)
    indexJ      = np.zeros(nnz_values)
    valuesN     = np.zeros(nnz_values)
    valuesdNdu  = np.zeros(nnz_values)
    valuesdNdv  = np.zeros(nnz_values)

    n = len(U) -1 -p   
    m = len(V) -1 -q 
    nbf = n*m 
    l = 0 
    for k in range(nb_u_values):
        spanu = findspan(n,p,u[k],U)
        spanv = findspan(m,q,v[k],V)
        Nu    = derbasisfuns(spanu,p,U,1,u[k])
        Nv    = derbasisfuns(spanv,q,V,1,v[k])  
        for j in range(q+1):
            for i in range(p+1):   
                valuesN[l]     = Nv[0][j]*Nu[0][i] 
                valuesdNdu[l]  = Nv[0][j]*Nu[1][i] 
                valuesdNdv[l]  = Nv[1][j]*Nu[0][i] 
                indexI[l] = k  
                indexJ[l] = spanu - p + i + (spanv - q + j)*n 
                l = l+1 
    phi     = sps.csc_matrix(( valuesN, (indexI,indexJ)), shape = (nb_u_values,nbf) )
    dphidu  = sps.csc_matrix(( valuesdNdu, (indexI,indexJ)), shape = (nb_u_values,nbf) )
    dphidv  = sps.csc_matrix(( valuesdNdv, (indexI,indexJ)), shape = (nb_u_values,nbf) )
    return phi, dphidu, dphidv

def Get3dBasisFunctionsAtPts(u,v,w,U,V,W,p,q,r):
    nb_u_values = len(u);
    nnz_values = nb_u_values*(p+1)*(q+1)*(r+1)
    indexI      = np.zeros(nnz_values)
    indexJ      = np.zeros(nnz_values)
    valuesN     = np.zeros(nnz_values)
    valuesdNdu  = np.zeros(nnz_values)
    valuesdNdv  = np.zeros(nnz_values)   
    valuesdNdw  = np.zeros(nnz_values)
    
    n = len(U)-1-p
    m = len(V)-1-q 
    l = len(W)-1-r 
    
    nbf = n*m*l
    ll=0
    for up in range(nb_u_values):
        spanu = findspan(n,p,u[up],U)
        spanv = findspan(m,q,v[up],V)
        spanw = findspan(l,r,w[up],W)
        Nu    = derbasisfuns(spanu,p,U,1,u[up])
        Nv    = derbasisfuns(spanv,q,V,1,v[up])
        Nw    = derbasisfuns(spanw,r,W,1,w[up])
        for k in range(r+1):
            for j in range(q+1):
                for i in range(p+1):
                    valuesN[ll]    = Nw[0][k]*Nv[0][j]*Nu[0][i]
                    valuesdNdu[ll] = Nw[0][k]*Nv[0][j]*Nu[1][i]
                    valuesdNdv[ll] = Nw[0][k]*Nv[1][j]*Nu[0][i]
                    valuesdNdw[ll] = Nw[1][k]*Nv[0][j]*Nu[0][i]
                    indexI[ll] = up
                    # Structured grid: i+j*nx+k*nx*ny
                    indexJ[ll] = spanu-p+i +(spanv-q+j)*n +(spanw-r+k)*n*m
                    ll = ll +1
    phi     = sps.csc_matrix(( valuesN,    (indexI,indexJ)), shape = (nb_u_values,nbf) )
    dphidu  = sps.csc_matrix(( valuesdNdu, (indexI,indexJ)), shape = (nb_u_values,nbf) )
    dphidv  = sps.csc_matrix(( valuesdNdv, (indexI,indexJ)), shape = (nb_u_values,nbf) )
    dphidw  = sps.csc_matrix(( valuesdNdw, (indexI,indexJ)), shape = (nb_u_values,nbf) )
    return phi, dphidu, dphidv, dphidw                     


def Get3dBasisFunctionsAtPtsWd(u,v,w,U,V,W,p,q,r):
    nb_u_values = len(u);
    nnz_values = nb_u_values*(p+1)*(q+1)*(r+1)
    indexI      = np.zeros(nnz_values)
    indexJ      = np.zeros(nnz_values)
    valuesN     = np.zeros(nnz_values)
    
    n = len(U)-1-p
    m = len(V)-1-q 
    l = len(W)-1-r 
    
    nbf = n*m*l
    ll=0
    for up in range(nb_u_values):
        spanu = findspan(n,p,u[up],U)
        spanv = findspan(m,q,v[up],V)
        spanw = findspan(l,r,w[up],W)
        Nu    = BasisFunc(spanu,u[up],p,U)
        Nv    = BasisFunc(spanv,v[up],q,V)
        Nw    = BasisFunc(spanw,w[up],r,W)
        for k in range(r+1):
            for j in range(q+1):
                for i in range(p+1):
                    valuesN[ll]    = Nw[k]*Nv[j]*Nu[i]
                    indexI[ll] = up
                    # Structured grid: i+j*nx+k*nx*ny
                    indexJ[ll] = spanu-p+i +(spanv-q+j)*n +(spanw-r+k)*n*m
                    ll = ll +1
    phi     = sps.csc_matrix(( valuesN,    (indexI,indexJ)), shape = (nb_u_values,nbf) )
    return phi         

    
 
def interpolateLinearly(xmin,xmax,lmin,lmax):
    a = (lmax-lmin)/(xmax-xmin)
    b = (lmin*xmax - lmax*xmin)/(xmax-xmin)
    return a,b  


def computeNOELEM3D(nxi, neta, nzeta, p, q , r ):
    """ Creates element connectivity for 3D B-splines """ 
    bf_xi_index   = (np.kron( np.ones(nxi), np.arange(p+1) ) + np.kron(np.arange(nxi), np.ones(p+1) )).reshape((nxi,p+1))
    bf_eta_index  = (np.kron( np.ones(neta), np.arange(q+1) ) + np.kron(np.arange(neta), np.ones(q+1) )).reshape((neta,q+1))
    bf_zeta_index = (np.kron( np.ones(nzeta), np.arange(r+1) ) + np.kron(np.arange(nzeta), np.ones(r+1) )).reshape((nzeta,r+1))
    
    noelem = np.kron(  np.ones_like(bf_zeta_index) ,  np.kron(np.ones_like(bf_eta_index)  ,bf_xi_index) )  + \
          (nxi+p)*np.kron(  np.ones_like(bf_zeta_index), np.kron( bf_eta_index, np.ones_like(bf_xi_index) ) )  + \
          (nxi+p)*(neta+q)*np.kron(  bf_zeta_index, np.kron( np.ones_like(bf_eta_index),np.ones_like(bf_xi_index) ) )  
    
    return noelem.astype('int32') 
    


# nxi = 4 ; neta = 2
# p = 2   ; q= 1 
# bf_xi_index  = (np.kron( np.ones(nxi), np.arange(p+1) ) + np.kron(np.arange(nxi), np.ones(p+1) )).reshape((nxi,p+1))
# bf_eta_index = (np.kron( np.ones(neta), np.arange(q+1) ) + np.kron(np.arange(neta), np.ones(q+1) )).reshape((neta,q+1))
# noelem = np.kron(np.ones_like(bf_eta_index),bf_xi_index) + (nxi+p)*np.kron(bf_eta_index, np.ones_like(bf_xi_index)) 

# nxi = 2 ; neta = 1 ; nzeta = 3 
# p = 2 ; q = 1 ; r = 3 
# bf_xi_index   = (np.kron( np.ones(nxi), np.arange(p+1) ) + np.kron(np.arange(nxi), np.ones(p+1) )).reshape((nxi,p+1))
# bf_eta_index  = (np.kron( np.ones(neta), np.arange(q+1) ) + np.kron(np.arange(neta), np.ones(q+1) )).reshape((neta,q+1))
# bf_zeta_index = (np.kron( np.ones(nzeta), np.arange(r+1) ) + np.kron(np.arange(nzeta), np.ones(r+1) )).reshape((nzeta,r+1))

# noelem = np.kron(  np.ones_like(bf_zeta_index) ,  np.kron(np.ones_like(bf_eta_index)  ,bf_xi_index) )  + \
#          (nxi+p)*np.kron(  np.ones_like(bf_zeta_index), np.kron( bf_eta_index, np.ones_like(bf_xi_index) ) )  + \
#          (nxi+p)*(neta+q)*np.kron(  bf_zeta_index, np.kron( np.ones_like(bf_eta_index),np.ones_like(bf_xi_index) ) )     
         
 
         


# noelemC = np.zeros((nxi*neta*nzeta,(p+1)*(q+1)*(r+1)))
 
# ie=0;  
# for k in range(nzeta):
#     for j in range(neta):
#         for i in range(nxi):
#             t=0
#             for kk in range(r+1):
#                 for jj in range(q+1):
#                     for ii in range(p+1):
#                         noelemC[ie,t] = i+ii+(j+jj)*(nxi+p)+(k+kk)*(nxi+p)*(neta+q) ; 
#                         t+=1
#             ie+=1

# print(np.max(noelem-noelemC))
# print(np.min(noelem-noelemC)) 
    
        
    
    
