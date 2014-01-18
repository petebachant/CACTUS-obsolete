"""
This module contains classes for constructing turbine geometry
It replaces CreateBlade.m and CreateStrut.m

"""

from numpy import zeros, ones, sqrt, pi, sin, cos
import numpy as np

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
#    O = origin[ones(np.size(v,1),1), :] 
#    O = origin
    # Why are we using "O" as a variable name?
    vO = v - O
    p = np.array([np.zeros(v.shape[0], 1), vO]).reshape(1, -1)   
#    p = [zeros(np.size(v,1),1),vO]
    
    # Rotation quaternion and conjugate
#    q = np.matrix([cos(theta/2), [nR*sin(theta/2)]])
    q = np.array([cos(theta/2), nR*sin(theta/2)]).reshape(1, -1)
    qbar = np.array([q[0], - q[(2 -1):4]]).reshape(1, -1)
#    qbar=[q[0], -q[1:3]]
    
    # These should be matrices
#    QL = np.matrix([[q[0], -q[1], -q[2], -q[3]],
#                    [q[1],  q[0], -q[3],  q[2]],
#                    [q[2],  q[3],  q[0], -q[1]],
#                    [q[3], -q[2],  q[1],  q[0]]])
#       
#    QbarR = np.matrix([[qbar[0], -qbar[1], -qbar[2], -qbar[3]],
#                       [qbar[1],  qbar[0],  qbar[3], -qbar[2]],
#                       [qbar[2], -qbar[3],  qbar[0],  qbar[1]],
#                       [qbar[3],  qbar[2], -qbar[1],  qbar[0]]])
    QL = np.array([q[0], -q[1], -q[2], -q[3], 
                   q[1], q[0], -q[3], q[2], 
                   q[2], q[3], q[0], - q[1], 
                   q[3], - q[2], q[1], q[0]]).reshape(1, -1)
    QbarR = np.array([qbar[0], -qbar[1], -qbar[2], -qbar[3], 
                      qbar[1], qbar[0], qbar[3], -qbar[2], 
                      qbar[2], -qbar[3], qbar[0], qbar[1], 
                      qbar[3], qbar[2], -qbar[1], qbar[0]]).reshape(1, -1)
    # Rotate p
    pR = p*(QbarR*QL).transpose()
#    vR = pR[:,2:4] + O
    vR = pR[:, (2 -1):4] + O
    return vR

class Blade(object):
    """
    Create empty blade with specified number of elements, n_elem. Quarter
    chord line and section normal and tangential vectors defined at zero
    turbine rotation phase.
    
    QC: Quarter chord location (over ref. radius) for each element end (size n_elem + 1).\n
    n: Blade section normal vector for each element end (size n_elem + 1). Defines
    positive angle of attack (AOA positive when relative velocity component
    is positive in the normal direction.\n
    t: Blade section tangent vector (must use rearward chord line) for each element end
    (size n_elem + 1).\n
    CtoR: Chord to ref. radius for each element end (size n_elem + 1).\n
    AreaR: Element area over ref. radius squared for each element (size n_elem).\n
    iSect: Airfoil section index for each element (size n_elem). Used in
    CACTUS to identify the airfoil data tables (defined in the CACTUS input
    file) to use with that element.\n
    
    """
    def __init__(self, n_elem):
        self.n_elem = n_elem
        self.QCx = zeros(1, self.n_elem+1)
        self.QCy = zeros(1, self.n_elem+1)
        self.QCz = zeros(1, self.n_elem+1)
        self.nx = zeros(1, self.n_elem+1)
        self.ny = zeros(1, self.n_elem+1)
        self.nz = zeros(1, self.n_elem+1)
        self.tx = zeros(1, self.n_elem+1)
        self.ty = zeros(1, self.n_elem+1)
        self.tz = zeros(1, self.n_elem+1)
        self.CtoR = zeros(1, self.n_elem+1)
        self.AreaR = zeros(1, self.n_elem)
        self.iSect = ones(1, self.n_elem)
        
    def rotate(self, theta, nR, origin):
        """Rotates blade structure around normal vector nR (size 1 x 3) through
        angle theta (rad), using specified origin point (size 1 x 3)."""
        # Rotate QC locations
        QC = [self.QCx, self.QCy, self.QCz]
        QCR = quatrot(QC.transpose(), theta, nR, origin)
        self.QCx = QCR[:,0].transpose()
        self.QCy = QCR[:,1].transpose()
        self.QCz = QCR[:,2].transpose()
        
        # Rotate n and t vectors
        t = [self.tx, self.ty, self.tz]
        tR = quatrot(t.transpose(),theta, nR, [0,0,0])
        self.tx = tR[:,0].transpose()
        self.ty = tR[:,1].transpose()
        self.tz = tR[:,2].transpose()
        
        n = [self.nx, self.ny, self.nz]
        nR = quatrot(n.transpose(), theta, nR, [0,0,0])
        self.nx = nR[:,0].transpose()
        self.ny = nR[:,1].transpose()
        self.nz = nR[:,2].transpose()
        
    def update_area(self):
        pass
        # Put the code from UpdateBElemArea.m here


class Strut(object):
    """
    Create empty strut with specified number of elements, n_elem. Strut element
    locations defined at zero turbine rotation phase.
    
    SE: Location (over ref. radius) for each element end (size n_elem + 1).\n
    CtoR: Chord to ref. radius for each element end (size n_elem + 1).\n
    AreaR: Element area over ref. radius squared for each element (size n_elem).\n
    TtoC: Strut thickness to chord ratio.\n
    BIndS: Index of the blade to which the first strut element is attached.\n
    EIndS: Index of the element on blade BIndS where the first strut element is attached.\n
    BIndE: Index of the blade to which the last strut element is attached.\n
    EIndE: Index of the element on blade BInd where the last strut element is attached.\n
       Blade and element indicies used in CACTUS to identify the relevant
       blade data to use in calculating strut-blade interference drag. \n
       For struts that are attached to the rotor shaft at one end (not to
       another blade), set the appropriate BInd and EInd values to zero. \n
   """
    def __init__(self, n_elem):
        self.n_elem = n_elem
        self.SEx = zeros(1,n_elem+1)
        self.SEy = zeros(1,n_elem+1)
        self.SEz = zeros(1,n_elem+1)
        self.CtoR = zeros(1,n_elem+1)
        self.AreaR = zeros(1,n_elem)
        self.TtoC = 0
        self.BIndS = 0
        self.EIndS = 0
        self.BIndE = 1
        self.EIndE = 1
        
    def rotate(self, theta, nR, Origin):
        """
        Rotates strut structure around normal vector nR (size 1 x 3) through
        angle Theta (rad), using specified origin point (size 1 x 3).
        """
        # Rotate element locations
        P = [self.SEx, self.SEy, self.SEz]
        PR = quatrot(P.conj().transpose(), Theta, nR, Origin)
        self.SEx = PR[:, 0].conj().transpose()
        self.SEy = PR[:, 1].conj().transpose()
        self.SEz = PR[:, 2].conj().transpose()
#        SR = S
        
    def update_area(self):
        pass


class Turbine(object):
    """
    Creates a CACTUS turbine geometry structure.
    
    n_blade: Number of blades. \n
    n_belem: Number of elements per blade. \n
    n_strut: Number of struts.
    n_selem: Number of elements per strut.
    ref_r: Reference radius (ft) defining the scale of the turbine geometry.
    Will be used to normalize other distances input directly to CACTUS, and
    also used when outputing dimensional results from CACTUS.
    rot_n: Normal vector (size 1 x 3) of turbine rotation axis. Input value
    used as default when turb_type is empty, but will be overwritten if a turb_type 
    is selected...
    rot_p: Origin point (size 1 x 3) on turbine rotation axis. Input value
    used as default when turb_type is empty, but will be overwritten if a turb_type 
    is selected...
    ref_ar: Reference frontal area scaled by ref_r^2. Used for
    force/torque/power coefficient normalization in CACTUS. Input value
    used as default when turb_type is empty, but will be overwritten if a turb_type 
    is selected...
    turb_type: Recognized generic turbine type string (see code below). Input empty to 
    just create an empty turbine structure. If input string is
    recognized, will fill arrays with actual data for given turbine type
    using additional arguments defined below.
    varargin: Additional args (comma separated) used for recognized turbine
    types (see comments in code below for definition).
    """
    def __init__(self, n_blade, n_belem, n_strut, n_selem, ref_r, rot_n,
                 rot_p, ref_ar, turb_type=None, varargin=None):
        
        self.n_blade = n_blade
        self.n_strut=n_strut
        self.rot_n = rot_n/sqrt(np.sum(rot_n**2)) # force normalize
        self.rot_p = rot_p
        self.ref_ar = ref_ar
        self.ref_r = ref_r
        self.turb_type = turb_type
        self.blades = []
        self.struts = []

        # Create blades
        for i in xrange(n_blade):
            self.blades.append(Blade(n_belem))

        # Create struts
        for i in xrange(n_strut):
            self.struts.append(Strut(n_selem)) 


        # Fill geometry if turbine type recognized
        if turb_type == 'VAWT':
            """
            Cross-flow turbine generator for a vertical axis wind turbine (VAWT) with either straight or parabolic blades.
            For n_strut > 0, struts will be evenly distributed amongst blades (n_strut must be a multiple of n_blade) and 
            along rotation axis from the center to the tips.
            Additional arguments:
                REqR: Equitorial radius to reference radius ratio
                CR: Blade chord to equitorial radius ratio
                HR: Turbine height to equitorial radius ratio
                eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
                BShape: 0 for straight blades, 1 for parabolic blades
                CRs: Strut chord to equitorial radius ratio (only used if n_strut > 0)
                TCs: Strut thickness to chord ratio (only used if n_strut > 0)
            """
        # Get vars
            
#        if len(varargin) < 7:
#            error('Not enough inputs for selected turbine type')
#        end 
#        REqR=varargin{1};
#        CR=varargin{2};
#        HR=varargin{3};
#        eta=varargin{4};
#        BShape=varargin{5};
#        CRs=varargin{6};
#        TCs=varargin{7};
#        
#        % Ref to reference radius
#        CR=CR*REqR;
#        HR=HR*REqR;
#        
#        % Set rotation axis along y
#        T.RotN=[0,1,0];
#        T.RotP=[0,0,0];
#        
#        % Radius ratio function
#        yB=linspace(0,HR,n_belem+1);
#        if BShape
#            % parabolic blades
#            rr=REqR*(1-4*(yB/HR-.5).^2);
#            % Frontal area normalized by ref_r^2
#            T.ref_ar=2*(REqR*HR-1/3*HR);
#            
#            % Fill first blade
#            deltac=(eta-.25)*CR;
#            T.B(1).CtoR=CR*ones(1,n_belem+1);
#            T.B(1).tx=ones(1,n_belem+1);
#            T.B(1).ty=zeros(1,n_belem+1);
#            T.B(1).tz=zeros(1,n_belem+1);
#            % normal vector (machine inward)
#            drdy=-8/HR*(yB/HR-1/2);
#            n=[zeros(1,n_belem+1);drdy;ones(1,n_belem+1)];
#            nmag=sqrt(sum(n.^2));
#            n=n./nmag(ones(1,3),:);
#            Ind=find(n(3,:)<0);
#            n(:,Ind)=-n(:,Ind);
#            T.B(1).nx=n(1,:);
#            T.B(1).ny=n(2,:);
#            T.B(1).nz=n(3,:);
#            T.B(1).QCx=-deltac*ones(1,n_belem+1);
#            T.B(1).QCy=yB;
#            T.B(1).QCz=-rr;
#        else
#            % straight blades
#            rr=REqR*ones(size(yB));
#            % Frontal area normalized by ref_r^2
#            T.ref_ar=2*REqR*HR;
#            
#            % Fill first blade
#            deltac=(eta-.25)*CR;
#            T.B(1).CtoR=CR*ones(1,n_belem+1);
#            T.B(1).tx=ones(1,n_belem+1);
#            T.B(1).ty=zeros(1,n_belem+1);
#            T.B(1).tz=zeros(1,n_belem+1);
#            % normal vector (machine inward)
#            T.B(1).nx=zeros(1,n_belem+1);
#            T.B(1).ny=zeros(1,n_belem+1);
#            T.B(1).nz=ones(1,n_belem+1);
#            T.B(1).QCx=-deltac*ones(1,n_belem+1);
#            T.B(1).QCy=yB;
#            T.B(1).QCz=-rr;
#        end
#    
#        % Calc element area
#        T.B(1)=UpdateBElemArea(T.B(1));
#        
#        % Copy and rotate for other blades
#        Phase=linspace(0,2*pi,n_blade+1);
#        for i=2:n_blade
#            T.B(i)=RotateBlade(T.B(1),Phase(i),T.RotN,T.RotP);
#        end
#        
#        % Fill struts on first blade
#        if mod(n_strut,n_blade)~=0
#            error('Number of struts must be a multiple of the number of blades for the ''VAWT'' input type.');
#        end
#        NSpB=round(n_strut/n_blade);
#        yS=linspace(0,HR,NSpB+2);
#        yS=yS(2:end-1);
#        rrS=interp1(yB,rr,yS);
#        yC=(yB(2:end)+yB(1:end-1))/2;
#        
#        for i=1:NSpB
#            T.S(i).SEx=zeros(1,n_selem+1);
#            T.S(i).SEy=yS(i)*ones(1,n_selem+1);
#            T.S(i).SEz=-linspace(0,rrS(i),n_selem+1);
#            T.S(i).CtoR=CRs*ones(1,n_selem+1);
#            T.S(i).TtoC=TCs;
#            T.S(i).BIndS=0;
#            T.S(i).EIndS=0;
#            T.S(i).BIndE=1;
#            [m,T.S(i).EIndE]=min(abs(yC-yS(i)));
#            
#            % Calc element area
#            T.S(i)=UpdateSElemArea(T.S(i));
#        end
#        
#        % Copy and rotate for other blades
#        for i=2:n_blade
#            for j=1:NSpB
#                SInd=(i-1)*NSpB+j;
#                T.S(SInd)=RotateStrut(T.S(j),Phase(i),T.RotN,T.RotP);
#                T.S(SInd).BInd=i;
#            end
#        end
#    
#        
#    elseif strcmp(turb_type,'HAWT')==1
#        
#        % Axial-flow turbine generator for a horizontal axis wind turbine (HAWT).
#    	% Additional arguments:
#        % RMaxR: Turbine radius to reference radius ratio
#        % HubRR: Hub radius to turbine radius ratio
#        % CR: Blade chord to turbine radius ratio (n_belem+1 elements ordered root to tip)
#        % bTwist: Blade planform twist at each element end (deg, w.r.t. rotor disk plane, positive LE into the wind (-x), n_belem+1 elements ordered root to tip)
#        % bi: Blade planform incidence (deg, w.r.t. rotor disk plane, positive LE into the wind (-x))
#        % eta: Blade mount point ratio ((distance behind leading edge of the blade mount point) / (chord))
#        % bCone: Blade coning angle (deg, positive tip into the wind (-x))
#        % Tilt: Rotor tilt angle (deg, positive windward axis tilted up)
#    
#        % Get vars
#        if length(varargin)<8
#            error('Not enough inputs for selected turbine type');
#        end 
#        RMaxR=varargin{1};
#        HubRR=varargin{2};
#        CR=varargin{3};
#        bTwist=varargin{4};
#        bi=varargin{5};
#        eta=varargin{6};
#        bCone=varargin{7};
#        Tilt=varargin{8};
#        
#        % Ref to reference radius
#        CR=CR*RMaxR;
#        HubRR=HubRR*RMaxR;
#        
#        % Set rotation axis along x
#        T.RotN=[1,0,0];
#        T.RotP=[0,0,0];
#        
#        % Radius ratio function
#        rB=linspace(HubRR,RMaxR,n_belem+1);
#        % Frontal area normalized by ref_r^2
#        T.ref_ar=pi*RMaxR^2;
#    
#        % Fill first blade
#        deltac=(eta-.25)*CR(1);
#        T.B(1).QCx=zeros(1,n_belem+1);
#        T.B(1).QCy=rB;
#        T.B(1).QCz=deltac*ones(1,n_belem+1);
#        T.B(1).CtoR=CR;
#        sTwist=sin(bTwist/180*pi);
#        cTwist=cos(bTwist/180*pi);
#        T.B(1).tx=sTwist;
#        T.B(1).ty=zeros(1,n_belem+1);
#        T.B(1).tz=-cTwist;
#        % normal vector (machine rearward (x))
#        T.B(1).nx=cTwist;
#        T.B(1).ny=zeros(1,n_belem+1);
#        T.B(1).nz=sTwist;
#    
#        % Calc element area
#        T.B(1)=UpdateBElemArea(T.B(1));
#        
#        % Rotate through incidence and coning angle
#        T.B(1)=RotateBlade(T.B(1),bi/180*pi,[0,-1,0],[0,0,0]);
#        T.B(1)=RotateBlade(T.B(1),bCone/180*pi,[0,0,1],[0,0,0]);
#        
#        % Copy and rotate for other blades
#        Phase=linspace(0,2*pi,n_blade+1);
#        for i=2:n_blade
#            T.B(i)=RotateBlade(T.B(1),Phase(i),T.RotN,T.RotP);
#        end
#        
#        % Rotate turbine through tilt angle
#        T=RotateTurbine(T,Tilt/180*pi,[0,0,-1],[0,0,0]);
#        
#    end
    
    def rotate(self, theta, nvec, origin):
        """Rotate the turbine."""
        
        # Rotate turbine rotation axis vector
        self.rotn = quatrot(self.rotn, theta, nvec, (0,0,0))
        
        # Rotate turbine rotation origin point
        self.rotp = quatrot(self.rotp, theta, nvec, origin)
        
        for blade in self.blades:
            blade.rotate(theta, nvec, origin)
            
        for strut in self.struts:
            strut.rotate(theta, nvec, origin)
            
    def writefile(self):
        """Writes the *.geom file."""
        pass
    
    def plot(self, options=None):
        """Plot turbine geometry."""
        pass
    
    
if __name__ == "__main__":
    vector = (1, 0, 0)
    nvector = (0, 1, 0)
    origin = (0, 0, 0)
    theta = pi/2
    print quatrot(vector, theta, nvector, origin)