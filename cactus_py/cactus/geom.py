"""
This module contains classes for constructing turbine geometry.
It replaces the /CreateGeom folder of *.m files
"""
from __future__ import division, print_function
from numpy import zeros, ones, sqrt, pi, sin, cos
import numpy as np

vawt_defaults = {"REqR" : 1.0,
                 "CR" : 0.07408,
                 "HR" : 2.64,
                 "eta" : 0.42,
                 "BShape" : "parabolic",
                 "CRs" : 0.07408,
                 "TCs" : 0.15}
hawt_defaults = {"RMaxR" : 1.0,
                 "HubRR" : 0.0334,
                 "CR" : [0.0394, 0.1332, 0.1133, 0.0926, 0.0759, 0.0551],
                 "bTwist" : [0.0 ,20.0, 4.715, 0.494, -0.92, -2.5],
                 "bi" : 8.5,
                 "eta" : 0.3,
                 "bCone" : 0.0,
                 "Tilt" : 0.0}

def quatrot(v, theta, nR, origin):
    """Perform rotation of vector(s) v around normal vector nR using the
    quaternion machinery.
    
    v: input vector(s) (can be m x 3 for m vectors to rotate)
    Theta: rotation angle (rad)
    nR: normal vector around which to rotate
    Origin: origin point of rotation
    
    vR: Rotated vector(s) (size m x 3 for m input vectors)"""
    
    # Convert to numpy arrays
    v, nR, origin = np.array(v), np.array(nR), np.array(origin)
    # Force normalize nR
    nR = nR/np.sqrt(np.sum(nR**2))    
    # Quaternion form of v
    # Make origin array with same number of rows as number of input vectors
    oarray = origin*ones(np.shape(v))
    vO = v - oarray
    if len(np.shape(v)) == 1:
        p = np.hstack((0, vO))
    else:
        p = np.hstack((zeros((np.shape(v)[0], 1)), vO))
    # Rotation quaternion and conjugate
    q = np.hstack((cos(theta/2), nR*sin(theta/2)))
    qbar = np.hstack((q[0], -q[1:]))
    QL = np.matrix([[q[0], -q[1], -q[2], -q[3]],
                    [q[1],  q[0], -q[3],  q[2]],
                    [q[2],  q[3],  q[0], -q[1]],
                    [q[3], -q[2],  q[1],  q[0]]])
    QbarR = np.matrix([[qbar[0], -qbar[1], -qbar[2], -qbar[3]],
                       [qbar[1],  qbar[0],  qbar[3], -qbar[2]],
                       [qbar[2], -qbar[3],  qbar[0],  qbar[1]],
                       [qbar[3],  qbar[2], -qbar[1],  qbar[0]]])
    # Rotate p
    p = np.matrix(p)
    pR = p*(QbarR*QL).conj().transpose()
    vR = pR[:, 1:] + oarray
    return vR

class Blade(object):
    """
    Create empty blade with specified number of elements, n_elem. Quarter
    chord line and section normal and tangential vectors defined at zero
    turbine rotation phase.
    
    QC: Quarter chord location (over ref. radius) for each element end 
        (size n_elem + 1).
    n:  Blade section normal vector for each element end (size n_elem + 1). 
        Defines positive angle of attack (AOA positive when relative velocity 
        component is positive in the normal direction.
    t:  Blade section tangent vector (must use rearward chord line) for each 
        element end (size n_elem + 1).
    CtoR: Chord to ref. radius for each element end (size n_elem + 1).
    AreaR: Element area over ref. radius squared for each element (size n_elem).
    iSect: Airfoil section index for each element (size n_elem). Used in
        CACTUS to identify the airfoil data tables (defined in the CACTUS input
        file) to use with that element.
    """
    def __init__(self, n_elem):
        self.n_elem = n_elem
        self.FlipN = 0
        # Element end geometry
        self.QCx = zeros((1, n_elem+1))
        self.QCy = zeros((1, n_elem+1))
        self.QCz = zeros((1, n_elem+1))
        self.tx = zeros((1, n_elem+1))
        self.ty = zeros((1, n_elem+1))
        self.tz = zeros((1, n_elem+1))
        self.CtoR = zeros((1, n_elem+1))
        # Element geometry
        self.PEx = zeros((1, n_elem))
        self.PEy = zeros((1, n_elem))
        self.PEz = zeros((1, n_elem))
        self.tEx = zeros((1, n_elem))
        self.tEy = zeros((1, n_elem))
        self.tEz = zeros((1, n_elem))
        self.nEx = zeros((1, n_elem))
        self.nEy = zeros((1, n_elem))
        self.nEz = zeros((1, n_elem))
        self.sEx = zeros((1, n_elem))
        self.sEy = zeros((1, n_elem))
        self.sEz = zeros((1, n_elem))
        self.ECtoR = zeros((1, n_elem))
        self.EAreaR = zeros((1, n_elem))
        self.iSect = ones((1, n_elem))
        
    def rotate(self, theta, nR, origin):
        """Rotates blade structure around normal vector nR (size 1 x 3) through
        angle theta (rad), using specified origin point (size 1 x 3)."""
        # Rotate QC locations
        QC = np.vstack((self.QCx, self.QCy, self.QCz))
        QCR = quatrot(QC.conj().transpose(), theta, nR, origin)
        self.QCx = QCR[:,0].conj().transpose()
        self.QCy = QCR[:,1].conj().transpose()
        self.QCz = QCR[:,2].conj().transpose()
        
        # Rotate t vector
        t = np.vstack((self.tx, self.ty, self.tz))
        tR = quatrot(t.conj().transpose(), theta, nR, [0,0,0])
        self.tx = tR[:,0].conj().transpose()
        self.ty = tR[:,1].conj().transpose()
        self.tz = tR[:,2].conj().transpose()
        
        # Calculate element geometry
        self.calc_element_geom()
        
    def calc_element_geom(self):
        """Calculates blade element geometry."""
        n_elem = self.n_elem
        FlipN = self.FlipN
        
        for i in range(n_elem):
            PE = np.hstack((self.QCx[0, i+1]+self.QCx[0, i], 
                            self.QCy[0, i+1]+self.QCy[0, i], 
                            self.QCz[0, i+1]+self.QCz[0, i]))/2
            sE = -np.hstack((self.QCx[0, i+1]-self.QCx[0, i], 
                             self.QCy[0, i+1]-self.QCy[0, i], 
                             self.QCz[0, i+1]-self.QCz[0, i])) 
            # nominal element spanwise direction set opposite to QC line in CACTUS
            sEM = sqrt(np.sum(sE**2))
            sE = sE/sEM
            tE = np.hstack((self.tx[0, i+1]+self.tx[0, i], 
                            self.ty[0, i+1]+self.ty[0, i], 
                            self.tz[0, i+1]+self.tz[0, i]))/2
            # Force tE normal to sE
            tE = tE - (tE*sE.conj().transpose())*sE
            tEM = sqrt(np.sum(tE**2))
            if tEM < 1e-10:
                raise RuntimeError('Error: Element t vector must not be \
                                   parallel to quarter chord line.')
            tE = tE/tEM
            # Calc normal vector
            nE = np.cross(sE, tE)
            nE = nE/sqrt(np.sum(nE**2))
            
            # Flip normal direction if requested
            if FlipN == 1:
                nE = -nE
                sE = -sE
            
            self.PEx[0, i] = PE[0]
            self.PEy[0, i] = PE[1]
            self.PEz[0, i] = PE[2]
            self.tEx[0, i] = tE[0]
            self.tEy[0, i] = tE[1]
            self.tEz[0, i] = tE[2]
            self.nEx[0, i] = nE[0]
            self.nEy[0, i] = nE[1]
            self.nEz[0, i] = nE[2]
            self.sEx[0, i] = sE[0]
            self.sEy[0, i] = sE[1]
            self.sEz[0, i] = sE[2]
            
            # Calc element area and chord
            S = np.array([-0.25, 0.75])
            SR = np.array([0.75, -0.25])
            # Calc quad area from two triangular facets
            Px = np.hstack([self.QCx[0, i] + S*self.CtoR[0, i]*self.tx[0, i], 
                            self.QCx[0, i+1] + SR*self.CtoR[0, i+1]*self.tx[0, i+1], 
                            self.QCx[0, i] - 0.25*self.CtoR[0, i]*self.tx[0, i]])
            Py = np.hstack([self.QCy[0, i] + S*self.CtoR[0, i]*self.ty[0, i], 
                            self.QCy[0, i+1] + SR*self.CtoR[0, i+1]*self.ty[0, i+1], 
                            self.QCy[0, i] - 0.25*self.CtoR[0, i]*self.ty[0, i]])
            Pz = np.hstack([self.QCz[0, i] + S*self.CtoR[0, i]*self.tz[0, i], 
                            self.QCz[0, i+1] + SR*self.CtoR[0, i+1]*self.tz[0, i+1], 
                            self.QCz[0, i] - 0.25*self.CtoR[0, i]*self.tz[0, i]])
            V = np.vstack([np.diff(Px), np.diff(Py), np.diff(Pz)])
            A1 = np.cross(V[:,0], V[:,1])/2
            A2 = np.cross(V[:,2], V[:,3])/2
            self.EAreaR[0, i] = sqrt(np.sum(A1**2))+sqrt(np.sum(A2**2))
            # Calc average element chord from area and span
            self.ECtoR[0, i] = self.EAreaR[0, i]/sEM


class Strut(object):
    """
    Create empty strut with specified number of elements, n_elem. Strut element
    locations defined at zero turbine rotation phase.
    
    SE:    Location (over ref. radius) for each element end (size n_elem + 1).
    CtoR:  Chord to ref. radius for each element end (size n_elem + 1).
    AreaR: Element area over ref. radius squared for each element (size n_elem).
    TtoC:  Strut thickness to chord ratio.
    BIndS: Index of the blade to which the first strut element is attached.
    EIndS: Index of the element on blade BIndS where the first strut element 
           is attached.
    BIndE: Index of the blade to which the last strut element is attached.
    EIndE: Index of the element on blade BInd where the last strut element 
           is attached.

    Blade and element indices used in CACTUS to identify the relevant
    blade data to use in calculating strut-blade interference drag.
    For struts that are attached to the rotor shaft at one end (not to
    another blade), set the appropriate BInd and EInd values to zero.
    """
    def __init__(self, n_elem):
        self.n_elem = n_elem
        self.TtoC = 0
        # Element end geometry
        self.MCx = zeros((1, n_elem+1))
        self.MCy = zeros((1, n_elem+1))
        self.MCz = zeros((1, n_elem+1))
        self.CtoR = zeros((1, n_elem+1))
        # Element geometry
        self.PEx = zeros((1, n_elem))
        self.PEy = zeros((1, n_elem))
        self.PEz = zeros((1, n_elem))
        self.sEx = zeros((1, n_elem))
        self.sEy = zeros((1, n_elem))
        self.sEz = zeros((1, n_elem))
        self.ECtoR = zeros((1, n_elem))
        self.EAreaR = zeros((1, n_elem))
        # Blade interference parameters
        self.BIndS = 0
        self.EIndS = 0
        self.BIndE = 1
        self.EIndE = 1
        
    def rotate(self, theta, nR, origin):
        """
        Rotates strut structure around normal vector nR (size 1 x 3) through
        angle Theta (rad), using specified origin point (size 1 x 3).
        """
        # Rotate element locations
        P = np.vstack((self.MCx, self.MCy, self.MCz))
        PR = quatrot(P.conj().transpose(), theta, nR, origin)
        self.MCx = PR[:, 0].conj().transpose()
        self.MCy = PR[:, 1].conj().transpose()
        self.MCz = PR[:, 2].conj().transpose()
        # Calculate element geometry
        self.calc_element_geom()
    
    def calc_element_geom(self):
        """Calculates strut element geometry."""
        for i in range(self.n_elem):
            PE = np.hstack((self.MCx[0, i+1] + self.MCx[0, i],
                            self.MCy[0, i+1] + self.MCy[0, i],
                            self.MCz[0, i+1] + self.MCz[0, i]))/2
            sE = np.hstack((self.MCx[0, i+1] - self.MCx[0, i],
                            self.MCy[0, i+1] - self.MCy[0, i],
                            self.MCz[0, i+1] - self.MCz[0, i]))
            sEM = sqrt(np.sum(sE**2))
            sE = sE/sEM
            self.PEx[0, i] = PE[0]
            self.PEy[0, i] = PE[1]
            self.PEz[0, i] = PE[2]
            self.sEx[0, i] = sE[0]
            self.sEy[0, i] = sE[1]
            self.sEz[0, i] = sE[2]
            # Calc element area and chord
            self.ECtoR[0, i] = (self.CtoR[0, i] + self.CtoR[0, i+1])/2
            self.EAreaR[0, i] = sEM*self.ECtoR[0, i]


class Turbine(object):
    """
    Creates a CACTUS turbine geometry structure.
    
    n_blade : int 
        number of blades.
    n_belem : int
        Number of elements per blade.
    n_strut : int
        Number of struts.
    n_selem : int
        Number of elements per strut.
    ref_r : float    
        Reference radius (ft) defining the scale of the turbine 
        geometry. Will be used to normalize other distances input 
        directly to CACTUS, and also used when outputing dimensional 
        results from CACTUS.
    rot_n : array
        Normal vector (size 1 x 3) of turbine rotation axis. Input value
        used as default when turb_type is empty, but will be overwritten 
        if a turb_type is selected.
    rot_p : array
        Origin point (size 1 x 3) on turbine rotation axis. Input value
        used as default when turb_type is empty, but will be overwritten 
        if a turb_type is selected.
    ref_ar : float
        Reference frontal area scaled by ref_r^2. Used for
        force/torque/power coefficient normalization in CACTUS. Input 
        value used as default when turb_type is empty, but will be 
        overwritten if a turb_type is selected...
    turb_type : string, optional
        Recognized generic turbine type string (see code below). 
        Input empty to just create an empty turbine structure. If input 
        string is recognized, will fill arrays with actual data for given 
        turbine type using additional arguments defined below.
    **kwargs : 
        Additional args used for recognized turbine types (see comments 
        below for definition).
             
    turb_type "VAWT"
    ----------------
    Cross-flow turbine generator for a vertical axis wind turbine 
    (VAWT) with either straight or parabolic blades.
    For n_strut > 0, struts will be evenly distributed amongst blades 
    (n_strut must be a multiple of n_blade) and 
    along rotation axis from the center to the tips.
    
    Additional arguments:
    
    REqR : 
        Equitorial radius to reference radius ratio
    CR :   
        Blade chord to equitorial radius ratio
    HR :   
        Turbine height to equitorial radius ratio
    eta : 
        Blade mount point ratio ((distance behind leading edge of the blade 
        mount point) / (chord))
    BShape : 
        Blade shape; "parabolic" for parabolic blades, None for straight
    CRs :  
        Strut chord to equitorial radius ratio (only used if n_strut > 0)
    TCs :  
        Strut thickness to chord ratio (only used if n_strut > 0)
         
    turb_type "HAWT"
    ----------------
    Axial-flow turbine generator for a horizontal axis wind turbine (HAWT).
    
    Additional arguments:
    
    RMaxR :  
        Turbine radius to reference radius ratio
    HubRR :  
        Hub radius to turbine radius ratio
    CR :     
        Blade chord to turbine radius ratio (n_belem+1 elements ordered root 
        to tip)
    bTwist : 
        Blade planform twist at each element end (deg, w.r.t. blade planform 
        plane (rotor disk plane when bi=0), positive LE into the wind (-x), 
        n_belem + 1 elements ordered root to tip)
    bi :     
        Blade planform incidence (deg, w.r.t. rotor disk plane, positive LE 
        into the wind (-x))
    eta : 
        Blade mount point ratio ((distance behind leading edge of the blade 
        mount point) / (chord))
    bCone :  
        Blade coning angle (deg, positive tip into the wind (-x))
    Tilt :   
        Rotor tilt angle (deg, positive windward axis tilted up)
    """
    def __init__(self, n_blade, n_belem, n_strut, n_selem, ref_r, rot_n,
                 rot_p, ref_ar, turb_type=None, **kwargs):
        
        self.n_blade = n_blade
        self.n_strut = n_strut
        rot_n = np.array(rot_n)
        self.rot_n = rot_n/sqrt(np.sum(rot_n**2)) # force normalize
        self.rot_p = np.array(rot_p)
        self.ref_ar = ref_ar
        self.ref_r = ref_r
        self.turb_type = turb_type
        self.blades = []
        self.struts = []

        # Create blades
        for i in range(n_blade):
            self.blades.append(Blade(n_belem))

        # Create struts
        for i in range(n_strut):
            self.struts.append(Strut(n_selem)) 

        # Fill geometry if turbine type recognized
        if self.turb_type == 'VAWT':
            # Get parameters from kwargs
            if len(kwargs) < 7:
                raise RuntimeError('Not enough inputs for selected turbine type')
            REqR = kwargs["REqR"]
            CR = kwargs["CR"]
            HR = kwargs["HR"]
            eta = kwargs["eta"]
            BShape = kwargs["BShape"]
            CRs = kwargs["CRs"]
            TCs = kwargs["TCs"]
            # Ref to reference radius
            CR = CR*REqR
            HR = HR*REqR
            # Set rotation axis along y
            self.rot_n = [0.0, 1.0, 0.0]
            self.rot_p = [0.0, 0.0, 0.0]
            # Radius ratio function
            yB = np.linspace(0, HR, n_belem+1)
            if BShape == "parabolic":
                # parabolic blades
                rr = REqR*(1.0 - 4.0*(yB/HR - 0.5)**2)
                # Frontal area normalized by ref_r^2
                self.ref_ar = 2*(REqR*HR - 1.0/3.0*HR)
                # Fill element end geometry
                deltac = (eta - 0.25)*CR
                self.blades[0].CtoR = CR*ones((1, n_belem+1))
                self.blades[0].tx = ones((1, n_belem+1))
                self.blades[0].ty = zeros((1, n_belem+1))
                self.blades[0].tz = zeros((1, n_belem+1))
                self.blades[0].QCx = -deltac*ones((1, n_belem+1))
                self.blades[0].QCy = yB*np.ones((1, n_belem+1))
                self.blades[0].QCz = -rr*np.ones((1, n_belem+1))
            else:
                # straight blades
                rr = REqR*ones(np.shape(yB))
                # Frontal area normalized by RefR^2
                self.ref_ar = 2*REqR*HR
                # Fill element end geometry
                deltac = (eta - 0.25)*CR
                self.blades[0].CtoR = CR*ones((1, n_belem+1))
                self.blades[0].tx = ones((1, n_belem+1))
                self.blades[0].ty = zeros((1, n_belem+1))
                self.blades[0].tz = zeros((1, n_belem+1))
                self.blades[0].QCx = -deltac*ones((1, n_belem+1))
                self.blades[0].QCy = yB*np.ones((1, n_belem+1))
                self.blades[0].QCz = -rr*np.ones((1, n_belem+1))
        
            # Calc element geom for first blade
            self.blades[0].calc_element_geom()
            
            # Copy and rotate for other blades
            Phase = np.linspace(0, 2*pi, n_blade+1)
            for i in range(1, n_blade):
                self.blades[i] = self.blades[0]
                self.blades[i].rotate(Phase[i], self.rot_n, self.rot_p)
            # Fill struts on first blade
            if float(n_strut) % n_blade != 0:
                raise RuntimeError('Number of struts must be a multiple of the\
                                    number of blades for the ''VAWT'' input type.')
            NSpB = np.round(n_strut/n_blade)
            yS = np.linspace(0, HR, NSpB+2)
            yS = yS[1:-1]
            rrS = np.interp(yS, yB, rr)
            yC = (yB[1:] + yB[:-1])/2
            
            for i in range(NSpB):
                # Fill element end geometry
                self.struts[i].MCx = zeros((1, n_selem+1))
                self.struts[i].MCy = yS[i]*ones((1, n_selem+1))
                self.struts[i].MCz = -np.linspace(0, rrS[i], n_selem+1)*np.ones((1, n_belem+1))
                self.struts[i].CtoR = CRs*ones((1, n_selem+1))
                self.struts[i].TtoC = TCs
                self.struts[i].BIndS = 0
                self.struts[i].EIndS = 0
                self.struts[i].BIndE = 1
                self.struts[i].EIndE = np.min(np.abs(yC - yS[i]))
                
                # Calc element geom
                self.struts[i].calc_element_geom()
            
            # Copy and rotate for other blades
            for i in range(1, n_blade):
                for j in range(NSpB):
                    SInd = i*NSpB + j
                    self.struts[SInd] = self.struts[j]
                    self.struts[SInd].rotate(Phase[i], self.rot_n, self.rot_p)
                    self.struts[SInd].BIndE = i

        elif self.turb_type=='HAWT':
            # Get vars
            if len(kwargs) < 8:
                raise RuntimeError('Not enough inputs for selected turbine type')
            RMaxR = kwargs["RMaxR"]
            HubRR = kwargs["HubRR"]
            CR = np.array(kwargs["CR"])
            bTwist = np.array(kwargs["bTwist"])
            bi = kwargs["bi"]
            eta = kwargs["eta"]
            bCone = kwargs["bCone"]
            Tilt = kwargs["Tilt"]
            
            # Ref to reference radius
            CR = CR*RMaxR
            HubRR = HubRR*RMaxR
            
            # Set rotation axis along x
            self.rot_n = np.array([1.0, 0.0, 0.0])
            self.rot_p = np.array([0.0, 0.0, 0.0])
            
            # Radius ratio function
            rB = np.linspace(HubRR, RMaxR, n_belem+1)
            # Frontal area normalized by RefR^2
            self.ref_ar = pi*RMaxR**2
        
            # Fill element end data for first blade
            deltac = (eta - 0.25)*CR[0]
            self.blades[0].QCx = zeros((1, n_belem+1))
            self.blades[0].QCy = rB*np.ones((1, n_belem+1))
            self.blades[0].QCz = deltac*ones((1, n_belem+1))
            self.blades[0].CtoR = CR*np.ones((1, n_belem+1))
            sTwist = sin(bTwist/180.0*pi)
            cTwist = cos(bTwist/180.0*pi)
            self.blades[0].tx = sTwist*np.ones((1, n_belem+1))
            self.blades[0].ty = zeros((1, n_belem+1))
            self.blades[0].tz = -cTwist*np.ones((1, n_belem+1))
            # Calc element geom for first blade
            self.blades[0].calc_element_geom()
            # Rotate through incidence and coning angle
            self.blades[0].rotate(bi/180.0*pi, [0,-1,0], [0,0,0])
            self.blades[0].rotate(bCone/180.0*pi, [0,0,1], [0,0,0])
            # Copy and rotate for other blades
            Phase = np.linspace(0, 2*pi, n_blade + 1)
            for i in range(1, n_blade):
                self.blades[i] = self.blades[0]
                self.blades[i].rotate(Phase[i], self.rot_n, self.rot_p)
            # Rotate turbine through tilt angle
            self.rotate(Tilt/180.0*pi, [0,0,-1], [0,0,0])
    
    def rotate(self, theta, nvec, origin):
        """Rotate the turbine."""
        # Rotate turbine rotation axis vector
        self.rot_n = quatrot(self.rot_n, theta, nvec, [0,0,0])
        # Rotate turbine rotation origin point
        self.rot_p = quatrot(self.rot_p, theta, nvec, origin)
        # Rotate blades and struts
        for blade in self.blades:
            blade.rotate(theta, nvec, origin)
        for strut in self.struts:
            strut.rotate(theta, nvec, origin)
            
    def writefile(self, name):
        """Writes the *.geom file."""
        if name[-5:] != ".geom":
            name = name + ".geom"
        with open(name, "w") as f:
            f.write("NBlade:   " + str(self.n_blade) + "\n")
            f.write("NStrut:   " + str(self.n_strut) + "\n")
            try: 
                np.shape(self.rot_n)[1]
                f.write("RotN:   {:e}   {:e}   {:e}\n".format(self.rot_n[0,0],
                        self.rot_n[0,1], self.rot_n[0,2]))
            except IndexError:
                f.write("RotN:   {:e}   {:e}   {:e}\n".format(self.rot_n[0],
                        self.rot_n[1], self.rot_n[2]))
            try: 
                np.shape(self.rot_p)[1]
                f.write("RotP:   {:e}   {:e}   {:e}\n".format(self.rot_p[0,0],
                        self.rot_p[0,1], self.rot_p[0,2]))
            except IndexError:
                f.write("RotP:   {:e}   {:e}   {:e}\n".format(self.rot_p[0],
                        self.rot_p[1], self.rot_p[2]))
            f.write("RefAR:   {:e}\n".format(self.ref_ar))
            f.write("RefR:   {:e}\n".format(self.ref_r))
            if self.turb_type:
                f.write("Type: " + self.turb_type + "\n")
            for n in range(self.n_blade):
                f.write("Blade " + str(n+1) + ":\n")
                f.write("\tNElem:   {:d}\n".format(self.blades[n].n_elem))
                f.write("\tFlipN:   {:d}\n".format(self.blades[n].FlipN))
                f.write("\tQCx:")
                for i in range(np.shape(self.blades[n].QCx)[1]):
                    f.write("  {:e}".format(self.blades[n].QCx[0,i]))
                f.write("\n")
                f.write("\tQCy:")
                for i in range(np.shape(self.blades[n].QCy)[1]):
                    f.write("  {:e}".format(self.blades[n].QCy[0,i]))
                f.write("\n")
                f.write("\tQCz:")
                for i in range(np.shape(self.blades[n].QCz)[1]):
                    f.write("  {:e}".format(self.blades[n].QCz[0,i]))
                f.write("\n")
                f.write("\ttx:")
                for i in range(np.shape(self.blades[n].tx)[1]):
                    f.write("  {:e}".format(self.blades[n].tx[0,i]))
                f.write("\n")
                f.write("\tty:")
                for i in range(np.shape(self.blades[n].ty)[1]):
                    f.write("  {:e}".format(self.blades[n].ty[0,i]))
                f.write("\n")
                f.write("\ttz:")
                for i in range(np.shape(self.blades[n].tz)[1]):
                    f.write("  {:e}".format(self.blades[n].tz[0,i]))
                f.write("\n")
                f.write("\tCtoR:")
                for i in range(np.shape(self.blades[n].CtoR)[1]):
                    f.write("  {:e}".format(self.blades[n].CtoR[0,i]))
                f.write("\n")
                f.write("\tPEx:")
                for i in self.blades[n].PEx[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tPEy:")
                for i in self.blades[n].PEy[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tPEz:")
                for i in self.blades[n].PEz[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\ttEx:")
                for i in self.blades[n].tEx[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\ttEy:")
                for i in self.blades[n].tEy[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\ttEz:")
                for i in self.blades[n].tEz[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tnEx:")
                for i in self.blades[n].nEx[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tnEy:")
                for i in self.blades[n].nEy[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tnEz:")
                for i in self.blades[n].nEz[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tsEx:")
                for i in self.blades[n].sEx[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tsEy:")
                for i in self.blades[n].sEy[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tsEz:")
                for i in self.blades[n].sEz[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tECtoR:")
                for i in self.blades[n].ECtoR[0,:]:
                    f.write("  {:e}".format(i))
                f.write("\n")
                f.write("\tEAreaR:")
                for i in self.blades[n].EAreaR[0,:]:
                    f.write("   {:e}".format(i))
                f.write("\n")
                f.write("\tiSect:")
                for i in self.blades[n].iSect[0,:]:
                    f.write("  {:d}".format(int(i)))
                f.write("\n")
            for n in range(self.n_strut):
                f.write("Strut " + str(n+1) + ":\n")
                f.write("\tNElem:  {:d}\n".format(self.struts[n].n_elem))
                f.write("\tTtoC:  {:e}\n".format(self.struts[n].TtoC))
                f.write("\tMCx:")
                for i in range(np.shape(self.struts[n].MCx)[1]):
                    f.write("  {:e}".format(self.struts[n].MCx[0,i]))
                f.write("\n")
                f.write("\tMCy:")
                for i in range(np.shape(self.struts[n].MCy)[1]):
                    f.write("  {:e}".format(self.struts[n].MCy[0,i]))
                f.write("\n")
                f.write("\tMCz:")
                for i in range(np.shape(self.struts[n].MCz)[1]):
                    f.write("  {:e}".format(self.struts[n].MCz[0,i]))
                f.write("\n")
                f.write("\tCtoR:")
                for i in range(np.shape(self.struts[n].CtoR)[1]):
                    f.write("  {:e}".format(self.struts[n].CtoR[0,i]))
                f.write("\n")
                f.write("\tPEx:")
                for i in range(np.shape(self.struts[n].PEx)[1]):
                    f.write("  {:e}".format(self.struts[n].PEx[0,i]))
                f.write("\n")
                f.write("\tPEy:")
                for i in range(np.shape(self.struts[n].PEy)[1]):
                    f.write("  {:e}".format(self.struts[n].PEy[0,i]))
                f.write("\n")
                f.write("\tPEz:")
                for i in range(np.shape(self.struts[n].PEz)[1]):
                    f.write("  {:e}".format(self.struts[n].PEz[0,i]))
                f.write("\n")
                f.write("\tsEx:")
                for i in range(np.shape(self.struts[n].sEx)[1]):
                    f.write("  {:e}".format(self.struts[n].sEx[0,i]))
                f.write("\n")
                f.write("\tsEy:")
                for i in range(np.shape(self.struts[n].sEy)[1]):
                    f.write("  {:e}".format(self.struts[n].sEy[0,i]))
                f.write("\n")
                f.write("\tsEz:")
                for i in range(np.shape(self.struts[n].sEz)[1]):
                    f.write("  {:e}".format(self.struts[n].sEz[0,i]))
                f.write("\n")
                f.write("\tECtoR:")
                for i in range(np.shape(self.struts[n].ECtoR)[1]):
                    f.write("  {:e}".format(self.struts[n].ECtoR[0,i]))
                f.write("\n")
                f.write("\tEAreaR:")
                for i in range(np.shape(self.struts[n].EAreaR)[1]):
                    f.write("  {:e}".format(self.struts[n].EAreaR[0,i]))
                f.write("\n")
                f.write("\tBIndS:  {:d}\n".format(self.struts[n].BIndS))
                f.write("\tEIndS:  {:d}\n".format(self.struts[n].EIndS))
                f.write("\tBIndE:  {:d}\n".format(self.struts[n].BIndE))
                f.write("\tEIndE:  {:d}\n".format(self.struts[n].BIndE))
    def plot(self, options=None):
        """Plot turbine geometry. This function still needs to be written."""
        pass
    
    
if __name__ == "__main__":
    n_blade = 2
    n_belem = 5
    n_strut = 0
    n_selem = 0
    ref_r = 18.15
    rot_n = [0, 0, 1]
    rot_p = [0, 0, 0]
    ref_ar = 1.0
    turb_opts = {"REqR" : 0.5,
                 "CR" : 0.1,
                 "HR" : 1.0,
                 "eta" : 0.5,
                 "BShape" : None,
                 "CRs" : 0.1,
                 "TCs" : 0.12}
    turb_opts = hawt_defaults
    turbine = Turbine(n_blade, n_belem, n_strut, n_selem, ref_r, rot_n,
                      rot_p, ref_ar, turb_type="HAWT", **turb_opts)
    print(turbine.blades[0].tz)
#    turbine.writefile("test")
